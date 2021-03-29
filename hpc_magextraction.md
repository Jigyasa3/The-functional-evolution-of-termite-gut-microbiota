## Extracting MAGs from metagenome contigs

```
#!/bin/bash
#SBATCH --job-name=metawrap
#SBATCH --partition=compute
#SBATCH --time=4-0
#SBATCH --mem=100G
#SBATCH --ntasks=3
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=megan_nt_%A-%a.out
#SBATCH -o out_megan.%j
#SBATCH -e err_megan.%j

#SBATCH --array 0-52
##num=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

# command to run
module load python/2.7.18

#call the folders-
IN_DIR="/bucket/BourguignonU/Jigs_backup/working_files/tagswitching/assembly/272_nov2"
READS_DIR="/bucket/BourguignonU/Jigs_backup/working_files/tagswitching/binning/reads/272_nov2"

cd ${IN_DIR}/
files=(*.fasta)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/_smallerkmer.fasta/}  #remove.fasta extension

#check if reads are unziped and renamed as READS_1.fastq and READS_2.fastq
#mv .local .local_backup #docker tries to read the python from .local instead of its own as /home/ directory is mapped to the docker.

#create a sample specific output directory-
OUT_DIR="<toadd>"
mkdir ${OUT_DIR}/${file2}

metawrap binning -o ${OUT_DIR}/${file2} -t 96 -a ${IN_DIR}/${file1} --concoct ${READS_DIR}/${file2}*fastq

#to transfer to /bucket/
#ls -d *_* > filename.txt
#while read line;do mkdir ${line};done<filename.txt ##in the bucket folder.
#while read line;do cp -r ${line}/concoct_bins/ ${OUT_DIR}/${line};done <filename.txt
#while read line;do ls ${line}/concoct_bins/*fa | wc -l;done < filename.txt ##to check how many bins formed

```


## CheckM software-completeness, contamination, 43-single copy marker genes

```
#!/bin/bash
#SBATCH --job-name=checkmtaxa
#SBATCH --partition=largemem
#SBATCH --time=4-0
#SBATCH --mem=120G
#SBATCH --ntasks=3

#SBATCH --output=check_nt_%A-%a.out
#SBATCH -o out_check.%j
#SBATCH -e err_check.%j

#SBATCH --array 0-26

#load modules-

module load Prodigal/2.6.2
module load python/3.7.3
module load pplacer/v1.1
module load hmmer/3.1b2

#num=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

IN_DIR="/bucket/BourguignonU/Jigs_backup/working_files/tagswitching/binning/bins/301_may"
cd ${IN_DIR}/

files=(*) #read the sample folders
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/.fasta/}


TAXA_DIR="/flash/BourguignonU/Jigs/binning/301_may"
#mkdir ${TAXA_DIR}/${file1}

checkm lineage_wf -x fa -u 5 -m 5 --tab_table ${IN_DIR}/${file1}/concoct_bins ${TAXA_DIR}/${file1} -t 16
#-m, --multi MULTI     maximum number of multi-copy phylogenetic markers before defaulting to domain-level marker set (default: 10)
#-u, --unique UNIQUE   minimum number of unique phylogenetic markers required to use lineage-specific marker set (default: 10)

```

## extract the 43-single copy marker genes 

```
#run "checkm_taxonomy.sh" script to get checkM marker genes out.

#get the marker genes from /checkM/storage/tree <WORKS!>
>for i in *-*;do grep ">" ${i}/storage/tree/*faa > markergenes_headers/${i}-markergenes-fastaheader.txt;done
>for i in *markergenes-fastaheader.txt;do awk -F":>" '{print $2}' ${i} > extracted-${i};done

#extract filename and contig name from the file- <FOR ONE FILE>
> while read line; do samplename=`echo ${line}| awk -F"&&" '{print $1}'`; contigname=`echo ${line} | awk -F"&&" '{print $2}'`; echo $samplenam
e;echo $contigname;done <extracted-301-90-markergenes-fastaheader.txt
#create a filename for each bin in the file-
>while read line; do samplename=`echo ${line}| awk -F"&&" '{print $1}'`; contigname=`echo ${line} | awk -F"&&" '{print $2}'`; touch ${samplena
me} ;done <extracted-301-90-markergenes-fastaheader.txt
#create a tab seperated file-
>awk -F"&&" '{print $1"\t"$2}' extracted-301-90-markergenes-fastaheader.txt > tabseparated-extracted-301-90-markergenes-fastaheader.txt
#extract the contig names from the file into samplename file-
>for i in 301-90.bin.*;do grep "$i" tabseparated-extracted-301-90-markergenes-fastaheader.txt > ${i};done
#seqtk on the contig name for each sample name-

#extract filename and contig name from the file- <FOR ALL THE FILES>
>for i in extracted*-markergenes-fastaheader.txt;do while read line; do samplename=`echo ${line}| awk -F"&&" '{print $1}'`; contigname=`echo ${line} | awk -F"&&" '{print $2}'`; touch ${samplename};done <${i};done
#check if all the bin files have been created using "touch" command-
>ls *bin* > test.txt
>sort all-bins-taxonomy-names.txt test.txt | uniq -u > unique-bins.txt # to compare "touch" files to all-bins file!
>for i in extracted-*;do awk -F"&&" '{print $1"\t"$2}' ${i} > tabseparated-${i};done
#get the bins with phyla level taxonomy info-
>for i in *.bin.*;do grep "$i" tabseparated-extracted-* > ${i};done
>for i in *.bin.*;do  sed -i 's/.*\t//g' ${i};done #keep only the contig info
#seqtk-
>for i in *.bin.*;do seqtk subseq ${IN_DIR}/${i}.fa-nucleotide.fna ${i} > markersequences-${i}.fna;done

##extract the bin contig names-
>for i in named-*;do file1=`echo ${i} | sed 's/named-//g'`;file2=`echo ${file1}| sed 's/.bin.*//g'`; cat named*${file2}* > contignames-${file2}.txt;done

##to get one fasta header per file in a filelist, but separate fasta files-
while read line; do samplename=`echo ${line}| awk '{print $1}'`; proteiname=`echo ${line}| awk '{print $2}'` ;for i in ${samplename}.fna_cds_from_genomic.fna;do filterbyname.sh in=${i} out=PF00164.20-${i} include=t names=${proteiname};done ;done < tabseparated-extracted-all-ncbi-PF00164.20.txt
#----------------------------------------------------

##OR-
#for each marker gene-
for i in *-*;do grep ">" ${i}/storage/tree/PF00164.20.masked.faa > markergenes_headers/PF00164.20-${i}-header.txt;done

#add sample name to the fasta header- <in fna files>
for i in *fna;do awk '/>/{sub(">","&"FILENAME"__");sub(/\.fasta-nucleotidesequence.fna/,x)}1' ${i} | sed 's/ #.*$//g' > fna_files/renamed-${i};done

#add sample name to header file-
for i in *txt;do awk -F"&&" '{print $1"\t"$2}' ${i} | sed 's/>//g' | sed 's/\t/__/g' > 2-${i};done

#run seqtk for each marker gene-
for i in renamed-*.fna ;do file1=`echo ${i}| sed 's/renamed-//g'`; file2=`echo ${file1} | sed 's/.bin.*$//g'`; seqtk subseq ${i} ${IN_DIR}/renamed-PF00164.20-${file2}-header.txt > ${OUT_DIR}/PF00164.20-${i};done

```


## checkm -qa to get completeness and contamination

```
#!/bin/bash
#SBATCH --job-name=checkmqa
#SBATCH --partition=largemem
#SBATCH --time=4-0
#SBATCH --mem=120G
#SBATCH --ntasks=3

#SBATCH --output=check_nt_%A-%a.out
#SBATCH -o out_check.%j
#SBATCH -e err_check.%j

#SBATCH --array 0-26

#load modules-

module load Prodigal/2.6.2
module load python/3.7.3
module load pplacer/v1.1
module load hmmer/3.1b2

#num=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

IN_DIR="/flash/BourguignonU/Jigs/binning/301_may2" #read checkmout folders
cd ${IN_DIR}/
mkdir checkm_qa

files=(*) #read the sample folders
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/.fasta/}


OUT_DIR="/flash/BourguignonU/Jigs/binning/301_may2/checkm_qa"
touch ${OUT_DIR}/${file1}_taxonomy.txt

checkm qa ${IN_DIR}/${file1}/lineage.ms ${IN_DIR}/${file1} -o 3 -f ${OUT_DIR}/${file1}_taxonomy.txt
#output format 3- summary of bin quality for increasingly basal lineage-specific marker sets

```


## GTDB taxonomy annotation

```
#!/bin/bash
#SBATCH --job-name=gtdb2
#SBATCH --partition=largemem
#SBATCH --time=3-0
#SBATCH --mem=150G
##SBATCH --ntasks=3

#SBATCH --output=gtdb2_nt_%A-%a.out
#SBATCH -o out_gtdb2.%j
#SBATCH -e err_gtdb2.%j
##SBATCH --array 0-358

CONTIG_DIR="/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/bin_taxonomy/all_my_bins/fastafile_withheaders/onefile"
#cd ${CONTIG_DIR}/

#files=(*fasta) #read the sample folders
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
#file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/.fasta/}

#load module-
module load python/3.7.3
module load hmmer/3.1b2
module load Prodigal/2.6.2
#other dependancies-FASTANI,MASH,PPLACER IN .bashrc file

OUT_DIR1="/flash/BourguignonU/Jigs/gtdb/gtdb_mags/identify"
OUT_DIR2="/flash/BourguignonU/Jigs/gtdb/gtdb_mags/align"
OUT_DIR3="/flash/BourguignonU/Jigs/gtdb/gtdb_mags/classify"

/home/j/jigyasa-arora/.local/bin/gtdbtk identify --genome_dir ${CONTIG_DIR} --out_dir ${OUT_DIR1} --extension fasta --cpus 2

/home/j/jigyasa-arora/.local/bin/gtdbtk align --identify_dir ${OUT_DIR1} --out_dir ${OUT_DIR2} --cpus 2

/home/j/jigyasa-arora/.local/bin/gtdbtk classify --genome_dir ${CONTIG_DIR} --align_dir ${OUT_DIR2} --out_dir ${OUT_DIR3} -x fasta --cpus 2


```

## Extract out information from MAGs according to Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) (Bowers et al 2017)
```
## extract rRNA sequences from each bin-
module load hmmer/3.1b2
module load ncbi-blast/2.7.1+
module load mafft/7.305

metaxa2 -i bin.1.fasta -o ${OUT_DIR}/metaxa2_bin.1 -f fasta --mode m --plus T
#--plus T : to specify that blast+ is used


#count of rRNA sequences per bin-
awk -F, '{print FILENAME","$0}' metaxa2_bin.1.summary.txt > 2-metaxa2_bin.1.summary.txt #add filename to output file with sufix ".summary.txt"
cat 2-metaxa2_*.summary.txt > all-metaxa2-bins.txt #concatenate them all

#------------------------------------------------------------------------------------------------------------------------------------------
##extract tRNA sequences from each bin-
module load python/3.7.3
file1="bin.1.fasta"
tRNAscan-SE -B -o ${OUT_DIR}/trna-scan-${file1} -m ${OUT_DIR}/${file1}.stats ${file1}

#-A : search for archaeal tRNAs
#-B : bacterial tRNAs only
#-G: use general tRNA model (cytoslic tRNAs from all 3 domains included)
#-m : save statistics summary

##NOTE- run tRNAscan with -A option for MAGs annotated as archaea.

grep "Total tRNAs:" bin.1.fasta.stats.txt > output-bin.1.fasta.stats.txt #get the line containing the final count of tRNA in the MAG.
cat output*.stats.txt > all-tRNAcounts-mags.txt #concatenate all the stats together in a file

#--------------------------------------------------------------------------------------------------------------------------------------------
##extract (a) total no. of contigs present in each MAG, (b) mean contig length, (c) maximum contig length for each MAG-

(a)
grep -c ">" named-301-92.bin.97.fasta > ../contig_stats/totalcount-named-301-92.bin.97.fasta.txt
awk -F, '{print FILENAME","$0}' totalcount-named-301-92.bin.97.fasta.txt > 2-totalcount-named-301-92.bin.97.fasta.txt #add filename to the output file
cat 2-totalcount-named-* > all-bins-totalcontigcount.txt


for (b) and (c)
module load emboss/6.6.0
infoseq -auto -only -name -length named-301-92.bin.97.fasta > info-named-301-92.bin.97.fasta.txt

(b)
awk '{x+=$2; next} END{print x/NR}' info-named-301-92.bin.97.fasta.txt > mean-info-named-301-92.bin.97.fasta.txt #get mean value from column2.
awk -F, '{print FILENAME","$0}' mean-info-named-301-92.bin.97.fasta.txt > 2-mean-info-named-301-92.bin.97.fasta.txt #add filename to the output file
cat 2-mean-info-named-* > all-bins-meancontiglength.txt

(c)
awk 'NR==1{max = $2 + 0; next} {if ($2 > max) max = $2;} END {print max}' info-named-301-92.bin.97.fasta.txt > max-info-named-301-92.bin.97.fasta.txt #get max. value from column2.
awk -F, '{print FILENAME","$0}' max-info-named-301-92.bin.97.fasta.txt > 2-max-info-named-301-92.bin.97.fasta.txt #add filename to the output file
cat 2-max-info-named-* > all-bins-maxcontiglength.txt

#---------------------------------------------------------------------------------------------------------------------------------------------
##coverage of contigs assembled into MAGs-


```
