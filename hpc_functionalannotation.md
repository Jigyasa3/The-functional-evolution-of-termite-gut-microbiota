## CAZY annotation

### a)map the protein sequences to CAZY database
```
#!/bin/bash
#SBATCH --job-name=hummrCAZY250above
#SBATCH --partition=compute
#SBATCH --time=7-0
#SBATCH --mem=80G
#SBATCH --ntasks=1
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=hmmr250_nt_%A-%a.out
#SBATCH -o out_hmmr250.%j
#SBATCH -e err_hmmr250.%j

##SBATCH --array 1-100
##num=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

#load modules-
module load hmmer/3.1b2

cd /work/BourguignonU/jigyasa/microcerotermes/functional_annotation/all_prokka_outputs

#files=(*.faa)
#file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/-prokka.faa/}
#file3=${file2/.txt/}


DB_DIR="/work/student/jigyasa-arora/BourguignonU_data/main_ID230+ID231_2/ID230/joinedfiles/K21-kate_assembly/prodigal/500above/cazy_output/hummer_cazy"


#for metagenomes use e-5 (dbCAN db recommendation, Ye et al, 2012)
#non-stringent hmmscan-
hmmscan -E 1.0e-5 -o 229-01-HMMoutputfile --tblout 229-01-HMM-persequence-output --domtblout 229-01-HMM-perdomain-output ${DB_DIR}/dbCAN-fam-HMMs.txt 229-01-prokka.faa

#NOTE- Even though a non-stringent evlaue cutoff used in the initial analysis, the final version has evalue <1e-30 (see below).
```

### b) get the hmmscan output file to a text delimited file

`sh hmmscan-parser.sh 229-01-HMM-perdomain-output > dbcan-229-01-HMM-perdomain-output`

### c) use a more stringent evalue cutoff for final CAZYme annotation. Get best hit per protein sequence

```
awk -F"," '$7<1e-30 {print $2"\t"$3"\t"$4"\t"$7"\t"$12}' dbcan-229-01-HMM-perdomain-output > evalue30-229-01-HMM-perdomain-output

export LC_ALL=C LC_LANG=C; sort -k2,2 -k5,5gr -k4,4g evalue30-cazy-output-229-01.txt > sorted-evalue30-cazy-output-229-01.txt #sort by proteinID, bitscore, evalue

for next in $(cut -f2 sorted-evalue30-cazy-output-229-01.txt | sort | uniq -u); do grep -w -m 1 "$next" sorted-evalue30-cazy-output-229-01.txt; done > top-sorted-evalue30-cazy-229-01.txt
# cut and sort proteinID column, grep fullwords (-w) and stop after first match (-m 1)

cat top-sorted-evalue30-* > allsamples-evalue30-cazymes.txt

```

### d) check if each protein is annotated only once, if not then get the best match out-

```
awk -F"\t" '{print $6"_"$2}' allsamples-evalue30-cazymes.txt | sed 's/"//g' > colnames.txt #get the proteinIDs per samples

sort colnames.txt | uniq -u > 2-colnames.txt #sort and get uniq IDs out
wc -l *colnames*
[1] 168640 2-colnames.txt
[2] 168640 colnames.txt

```



## KEGG annotation-Using default parameters. 

```
#!/bin/bash
#SBATCH --job-name=bedtoolsmap
#SBATCH --partition=largemem
#SBATCH --time=7-0
#SBATCH --mem=120G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=bedtoolsmap_nt_%A-%a.out
#SBATCH -o out_bedtoolsmap.%j
#SBATCH -e err_bedtoolsmap.%j

#SBATCH --array 0-225
#num=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

files=(*.faa)
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/-prokka.faa/}
#file3=${file2/.txt/}

#module load hmmer/3.1b2

/work/student/jigyasa-arora/kofamscan/bin/kofamscan-1.1.0/exec_annotation -o ${file1}-kopfamdetail ${file1} -p /work/student/jigyasa-arora/kofamscan/db/profiles -k /work/student/jigyasa-arora/kofamscan/db/ko_list --cpu 10 -f mapper --tmp-dir ${file1}-tmp
#-f mapper -one can be converted to heatmap
#-f detail -for tpm analysis
```

### a)if multiple KEGG IDs found per proteinID-
```
cat annotated-*kofam > all-metagenomes-kofam-annotation.txt

#in R-
kofam<-read.csv("all-metagenomes-kofam-annotation.txt",header=FALSE,sep="\t")

library(data.table)
setDT(kofam)
kofam2<-as.data.table(kofam)[, toString(V3), by = list(V1,V2)]
write.csv(kofam2,file="all-metagenomes-multiplekofamids-annotation.txt")

colnames(kofam2)<-c("samples","proteins","kofamid")
kofam2$sums<-count.fields(textConnection(kofam2$kofamid), sep = ",") #get the sum of kofamid column
library(dplyr)
kofam3<-kofam2%>%filer(sums>1)

library(stringr)
kofam3%>%filter(str_detect(kofamid, "K00297")) #check for each protein annotated as lignocellulose kegg IDs if they are annotated with other ids too.
write.csv(kofam3,file="final-all-metagenomes-multiplekofamids-annotation.txt") #after removing all the duplicate annotations. (see below)

##NOTE- proteins annotated with lignocellulose degradating pathway KEGG ids (in Suppmentary table S11) have single KEGG ID. Those that don't are listed below-
K00395 (aprB) is also annotated as K05337 (ferrodoxin) as aprB contains 4FE-4S cluster. Not removed.
K00266 (gltD) is also annotated as K05796 (electron transport protein HydN). They were removed.
K00297 (metF) is annotated as K03521 and K00548. They were removed.
K00925 (ack) is annotated as K00932. They were removed.
K00926 (arcC) is annotated as K07491. They were removed.
K01491 (folD) is annotated as K00949. They were removed.
K01491 (folD) is annotated as K01945. They were removed. 
K01915 (glnA) is annotated as K01778. They were removed. 
K01915 (glnA) is annotated as K09470. They were removed. 
K01938 (fhs) is annotated as K00288. They were removed. 
K00402(mcrG) is annotated as K03421. They were removed. 
K02591 (nifK) is annotated as  K02592. They were removed.
K01428 (ureB) is annotated as K01427. They were removed. 
K01429 (ureC) is annotated as K01438. They were removed.
K03320 (amtB) is annotated as K04751. They were removed.

```

## Pfam annotation for hydrogenases catalytic subunit classification. Proteins were annotated with PfamA database, and those annotated as PF00374 (NiFe) or PF02906 (Fe-Fe) were further classified via https://services.birc.au.dk/hyddb/ Hydrogenases database.

```
#!/bin/bash
#SBATCH --job-name=hmm
#SBATCH --partition=compute
#SBATCH --time=7-0
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=hmm_nt_%A-%a.out
#SBATCH -o out_hmm.%j
#SBATCH -e err_hmm.%j

#SBATCH --array 0-224
##num=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

#module load-
module load hmmer/3.1b2
module load  singularity/3.0.3

cd /work/BourguignonU/jigyasa/microcerotermes/functional_annotation/funcitonal_outputs/faa_files

files=(filtered-*faa) #protein files
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/filtered-/}
#file3=${file2/.txt/}

##USING FILTERED FILES-
#to remove asterisk containing fasta headers-
#grep -B1 "*" 230-59-150above-fraggene.fasta.faa > headers_with_asterisk.txt #-B1 : print a line before "*"
#grep "^>"  headers_with_asterisk.txt |sed 's/>//g' > headers_with_asterisk2.txt #to get the header info
#extract out the fasta headers-
#python3 remove_asterisk.py ${IN_DIR}/230-59-150above-fraggene.fasta.faa ${IN_DIR}/headers_with_asterisk2.txt > ${IN_DIR}/230-59-150above-fraggene.filtered.faa
#or
#filterbyname.sh in=229-${num}-prokka.faa out=remaining/filtered-229-${num}-prokka.faa include=f names=2-lines-with-astrisk-229-${num}-prokka.faa.txt substring=t ignorejunk


DB_DIR="/work/student/jigyasa-arora/BourguignonU_data/main_ID230+ID231_2/ID230/230run_nomismatch_joinedfiles/protein/pfam_hmm"

singularity exec /work/student/jigyasa-arora/pfa.sif pfam_scan.pl -fasta ${file1} -dir ${DB_DIR}/ -outfile pfam-output-${file2} -e_seq 1e-4 -cpu 10

#e-value cut-off used is 1e-4
```

### a) Pfam annotated proteins were further filtered by a stringent evalue cutoff.

`awk -F"\t" '$13 <1e-30 {print $0}' pfam-output-229-01.txt > evalue30-pfam-output-229-01.txt`

### b) Extract NiFe and FeFe catalytic subunits out-
```
grep "PF00374" evalue30-pfam-output-229-01.txt < nife-evalue30-pfam-output-229-01.txt
grep "PF02906" evalue30-pfam-output-229-01.txt < fefe-evalue30-pfam-output-229-01.txt
```

### c) Generate a final pfam annotated file for all samples-

```
cat evalue_30_pfam-output-* > allsamples-evalue30-pfamoutput.txt
awk '{print $1","$6","$16}' allsamples-evalue30-pfamoutput.txt > 3columns-allsamples-evalue30-pfamoutput.txt #extract proteinID, pfamID, sampleID out

#in R-
pfam<-read.csv("pfam<-read.csv("3columns-allsamples-evalue30-pfamoutput.txt",header=FALSE)

library(data.table)
setDT(pfam)
pfam2<-as.data.table(pfam)[, toString(V2), by = list(V1,V3)] #convert pfam annotation column toString using proteinID and sampleID columns as groups.
write.csv(pfam2,file="3columns-allsamples-multiplepfamids-evalue30-pfamoutput.txt")

##NOTE- for pfam annotation, one annotation per proteinID was not taken into consideration as only hydrogenases were examined from pfam annotation. Hydrogenase annotation via HydDB requires examination of protein sequences, removing the problem of double annotation of proteinID. Those protein sequences annotated as Hydrogenases were kept for statistical analysis in R. 
```
