## DIAMOND BLAST against the nr database to get prokaryotic and non-prokaryotic contigs out. LCA annotation of each contig using MEGAN6 command line version.

```
#!/bin/bash
#SBATCH --job-name=megan_diamondx_230
#SBATCH --partition=largemem
#SBATCH --time=14-0
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --mail-user=jigyasa-arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
 
#SBATCH --output=megan_nt_%A-%a.out
#SBATCH -o out_megan.%j
#SBATCH -e err_megan.%j
 
#SBATCH --array 1-26
num=$(printf "%02d" $SLURM_ARRAY_TASK_ID) #add two zeros infront

# command to run
export PATH=$PATH:/home/j/jigyasa-arora/local
module load ncbi-blast/2.6.0+
 
 
IN_DIR="/230-run"
OUT_DIR="/230-run"
DB_DIR="/work/student/jigyasa-arora/BourguignonU_data/main_ID230+ID231_2/ID230/joinedfiles/MEGAN/db"
TAX_DIR="/work/student/jigyasa-arora"

 
diamond blastx --query ${IN_DIR}/230-${num}_newname_scaffolds.fasta --db ${DB_DIR}/nr --outfmt 5 --out ${IN_DIR}/230-${num}_newname_scaffolds.fasta.xml --threads 15

#convert to txt format file-

xvfb-run /work/student/jigyasa-arora/tools/blast2lca -i ${IN_DIR}/230-${num}_newname_scaffolds.fasta.xml -f BlastXML -o megan-matches-${num}.txt -a2t /work/student/jigyasa-arora/prot_acc2tax-Nov2018X1.abin


#megan run on raw joined reads, no cut-offs
#outfmt 5 = xml format
#daa format is directly used in megan
#-alg -> to use lca algorithm


```

## Extract out prokaryotic contigs out
```

grep "d__Bacteria" megan-matches-229-01.txt > bacterial-megan-matches-229-01.txt
grep "d__Archaea" megan-matches-229-01.txt > archaeal-megan-matches-229-01.txt

cat bacterial-megan-matches-229-01.txt archaeal-megan-matches-229-01.txt > names-bacterial-2-newname-megan-matches-229-01.txt

seqtk subseq allcontigs-229-01.fasta names-bacterial-2-newname-megan-matches-229-01.txt > bacterial-contigs-229-01.fasta #extract prokaryotic contigs from nr
```

## GTDB database to annotate the prokaryotic contigs

```
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 3-0
#SBATCH --mem=150G
#SBATCH -e err_wood.%j

#SBATCH --array 0-226
#load the libraries-
module load ncbi-blast/2.10.0+
module load edirect/13.8
eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib` #to run!
echo 'eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`' >> ~/.profile #to run!
echo 'export MANPATH=$HOME/perl5/man:$MANPATH' >> ~/.profile #to run!

#call the folders-
IN_DIR="/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/222samples_taxonomy/gtdb/microbial_contigsequences"
DB_DIR="/flash/BourguignonU/DB/GTDB_Kraken/gtdb_r89_54k"
OUT_DIR="/flash/BourguignonU/Jigs/gtdb/kraken_all"

cd ${IN_DIR}/
files=(*fasta)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/filename-/}
#file3=${file2/.fasta/}

#seqtk subseq ${IN_DIR1}/${file1} ${IN_DIR2}/names-bacterial-2-newname-megan-matches-${file3}.txt > ${OUT_DIR}/bacterial-contigs-${file3}.fasta #extract prokaryotic contigs from nr
# Start 'myprog' with input from bucket, and output in /flash/

kraken2 --db ${DB_DIR} ${IN_DIR}/${file1} --output ${OUT_DIR}/kraken-output-${file1}.txt --report ${OUT_DIR}/kraken-report-${file1}.txt --report-zero-counts --use-mpa-style --use-names --threads 10

```

## generating the final taxonomy file, with annotation per contig-

```
files=(2-kraken-output-allprokaryotic_contigs-*.txt) #read the kraken output
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/2-kraken-output-allprokaryotic_contigs-/}

module load R/3.6.1
Rscript krakenoutputs.R ${CONTIG_DIR}/${file1} ${CONTIG_DIR}/kraken-report-allprokaryotic_contigs-${file2} ${OUT_DIR}/final-krakenoutput-${file2}
```

```
#save as "krakenoutputs.R"
#USAGE- Rscript [kraken] [taxonomy] [output_file]

args <- commandArgs(TRUE)
kraken <- args[1] #2-kraken file
taxonomy <- args[2] #report file
output_file <- args[3]


##getting a final kraken file-
#awk -F"\t" '{print $2"\t"$3}' kraken-output-301-92.txt > 2-kraken-output-301-92.txt

##in R-
library(dplyr)

kraken<-read.csv(file=kraken,sep="\t",header=FALSE)
kraken$V2<-gsub(" \\(taxid.*$","",kraken$V2)

taxonomy<-read.csv(file=taxonomy,header=FALSE,sep="\t")
taxonomy$V1<-as.character(taxonomy$V1)
taxonomy$lca<-gsub("^.*\\|","",taxonomy$V1)
taxonomy$lca<-gsub("^.*__","",taxonomy$lca)

#remove duplicates per lca-
taxonomy2<-taxonomy
levels<-taxonomy2%>%select(lca)%>%unique
levels<-as.vector(levels$lca)

taxonomy2$lca <- factor(taxonomy2$lca, levels = levels)
taxonomy2 <- taxonomy2[order(taxonomy2$lca), ]
taxonomy3<-taxonomy2[!duplicated(taxonomy2[,c('lca')]),]


#note- the two joining columns need to be character vectors-
kraken2<-merge(kraken,taxonomy3,by.x="V2",by.y="lca")
kraken3<-kraken2%>%filter(V2!="unclassified")
kraken3$V2.y<-NULL
kraken3$V2<-NULL

write.csv(kraken3,file=output_file)
```

## COG marker gene taxonomy

```
#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --partition=compute
#SBATCH --time=1-0
#SBATCH --mem=180G
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=iqtree_nt_%A-%a.out
#SBATCH -o out_iqtree.%j
#SBATCH -e err_iqtree.%j
#SBATCH --array 0-225

#load modules-
DNA_DIR="/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/all_prokka_outputs/all_fna_sequences_Dec2019"
cd ${DNA_DIR}/

files=(2-renamed-*fasta)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/2-renamed-/}
file3=${file2/-prokka.fasta/}

PROTEIN_DIR="/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/all_prokka_outputs/all_protein_sequences_Dec2019/with_filenames"
OUT_DIR="/flash/BourguignonU/Jigs/markers/allmarkers_re"

fetchMGs.pl -m extraction -x /home/j/jigyasa-arora/local/fetchMGs/bin -d ${DNA_DIR}/${file1} ${PROTEIN_DIR}/filename-${file3}-prokka.faa-smallername.faa -o ${OUT_DIR}/${file3}-fetchmoutput -t 16 ##bin directory where you cloned the repo.

```

## group data by COG markers. Taxonomic annotation of markers (.fna) files same way as contig taxonomic annotation done above.

```
while read line;do while read cogs;do cp ${line}/${cogs}*faa allfetchm_nucoutput/${line}-${cogs}.faa;done < allcogs.txt ;done <filesnames.txt

while read line;do while read cogs;do cp ${line}/${cogs}*fna allfetchm_nucoutput/${line}-${cogs}.fna;done < allcogs.txt ;done <filesnames.txt

##perform taxonomic annotation of COGs as done above using KRAKEN2
```

## check1- Does each contig has one annotation?

```
module load R/3.6.1
# in R-
krakenoutput<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/222samples_taxonomy/gtdb/krakenoutput/2-kraken-output-allprokaryotic_contigs-229-01.fasta.txt",header=FALSE,sep="\t")

library(dplyr)

contigs<-krakenoutput%>%select(V1)%>%unique()
nrow(contigs)
[1] 716080
nrow(krakenoutput)
[1] 716080
```

## check2-how many prokaryotic contigs were annotated-

```
#total no. of assembled contigs-
grep -c ">" all-samples-contigs.fasta
[1] 480971195

#nr microbial contigs-
wc -l allsamples-archaeaonly-newname-megan-matches.txt
[1] 1483566
wc -l allsamples-bacteriaonly-newname-megan-matches.txt
[1] 146547721 

#GTDB microbial contigs-
cat final-krakenoutput-* > allsamples-krakentaxonomy-feb2021.txt
sort contignames.txt | uniq -u > uniq-contignames.txt
wc -l uniq-contignames.txt
[1] 40254274

```

## check3-Any difference between GTDB contig annotation and GTDB marker gene annotation?

```
#in R
module load R/3.6.1

gtdb_taxonomy<-read.csv("allsamples-krakentaxonomy-feb2021.txt")
gtdb_taxonomy$V1.y<-gsub("\\|","_",gtdb_taxonomy$V1.y)
gtdb_taxonomy$phyla<-gsub("_c__.*$","",gtdb_taxonomy$V1.y)
gtdb_taxonomy$class<-gsub("_o__.*$","",gtdb_taxonomy$V1.y)
gtdb_taxonomy$order<-gsub("_f__.*$","",gtdb_taxonomy$V1.y)
gtdb_taxonomy$family<-gsub("_g__.*$","",gtdb_taxonomy$V1.y)

cog1<-read.csv("tpm_cogs_allsamples_feb2021_1000bpscontigs.csv") #the marker genes present in contigs >1000bps.
cog1$markertaxonomy<-gsub("\\|","_",cog1$markertaxonomy)                      
cog1$phyla<-gsub("_c__.*$","",cog1$markertaxonomy)                            
cog1$class<-gsub("_o__.*$","",cog1$markertaxonomy)
cog1$order<-gsub("_f__.*$","",cog1$markertaxonomy)
cog1$family<-gsub("_g__.*$","",cog1$markertaxonomy)

#joining at phyla level
gtdb_cog1_phyla<-merge(cog1,gtdb_taxonomy,by.x=c("fullnames","phyla"),by.y=c("V1.x","phyla"))
nrow(gtdb_cog1_phyla)
[1] 65818
nrow(cog1)
[1] 84463

#joining at class level
gtdb_cog1_class<-merge(cog1,gtdb_taxonomy,by.x=c("fullnames","class"),by.y=c("V1.x","class"))
nrow(gtdb_cog1_class)
[1] 64575
nrow(cog1)
[1] 84463

#joining at order level
gtdb_cog1_order<-merge(cog1,gtdb_taxonomy,by.x=c("fullnames","order"),by.y=c("V1.x","order"))
nrow(gtdb_cog1_order)
[1] 57970
nrow(cog1)
[1] 84463

#joining at family level
gtdb_cog1_family<-merge(cog1,gtdb_taxonomy,by.x=c("fullnames","family"),by.y=c("V1.x","family"))
nrow(gtdb_cog1_family)
[1] 48008
nrow(cog1)
[1] 84463

#---------------------------------------------------------------------------------------------------------------------------------------------------
cog2<-read.csv("tpm_cogs_allsamples_feb2021_1000bpscontigs_100counts1tpm.csv") #the marker genes with the same filters as functional genes.
cog2$markertaxonomy<-gsub("\\|","_",cog2$markertaxonomy)                      
cog2$phyla<-gsub("_c__.*$","",cog2$markertaxonomy)                            
cog2$class<-gsub("_o__.*$","",cog2$markertaxonomy)
cog2$order<-gsub("_f__.*$","",cog2$markertaxonomy)
cog2$family<-gsub("_g__.*$","",cog2$markertaxonomy)

#at family level-
gtdb_cog2_family<-merge(cog2,gtdb_taxonomy,by.x=c("fullnames","family"),by.y=c("V1.x","family"))
nrow(gtdb_cog2_family)
[1] 17563
nrow(cog2)
[1] 31928

write.csv(gtdb_cog2_family,file="tpm_cogs_allsamples_feb2021_1000bpscontigs_100counts1pm_similartocontigs.csv")
#-------------------------------------------------------------------------------------------------------------------------------------------------------
##how many marker genes are present on single contigs that donot match GTDB contig annotation?
cog2_notgtdb<-anti_join(cog2, gtdb_cog2_family, by=c("fullnames","family")) #markers on same contigs
#write.csv(cog2_notgtdb,file="cog2_notgtdbcontigs.csv")
nrow(cog2_notgtdb)
[1] 14365
nrow(cog2)
[1] 31928
14365/31928*100
[1] 44.99186

samecontigs<-merge(cog2_notgtdb,gtdb_cog2_family,by="fullnames")
> nrow(samecontigs)
[1] 13438
13438/14365*100
[1] 93.54682

##NOTE- At family level taxonomic annotation, 55% of marker gene taxonomy matches the GTDB contig taxonomy. Out of 44% that donot match, 99% of this 44% are present on the same contig as the ones that do match, making their presence redundant i.e. only 1% of data is lost.
Use "tpm_cogs_allsamples_feb2021_1000bpscontigs_100counts1pm_similartocontigs.csv" file for all statistical analysis.

```

