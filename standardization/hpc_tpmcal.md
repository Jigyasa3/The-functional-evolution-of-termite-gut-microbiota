Get Transcripts per million (TPM) per protein annotated as prokaryote based on GTDB taxonomy (hpc_taxonomy.md). This version of TPM calculation doesnot account for GC bias in data.

## TPM calculations for metagenome contigs

```
#!/bin/bash
#SBATCH --job-name=samtools-markdup
#SBATCH --partition=compute
#SBATCH --time=4-0
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=samtools-markdup_%A-%a.out
#SBATCH -o out_samtools-markdup.%j
#SBATCH -e err_samtools-markdup.%j

#SBATCH --array 0-26
#num=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
IN_DIR="/flash/BourguignonU/Jigs/tpm_2021/229samples"

cd ${IN_DIR}/
files=(*fasta)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/.fasta/}
file3=${file2/prok-filename-/}

#load modules-
module load bowtie2/2.2.6
module load samtools/1.9 #use the newest version of samtools

OUT_DIR="/flash/BourguignonU/Jigs/tpm_2021/229samples"
READS_DIR="/bucket/BourguignonU/Jigs_backup/working_files/230_run/trimmed_output/trimmed_reads/flash/BourguignonU/Jigs/zips/230_run/bucket/BourguignonU/Jigs_backup/working_files/230_run/trimmed_files/trimmed_reads"

#generate a sam file-
bowtie2-build ${IN_DIR}/${file1} ${IN_DIR}/${file2}
bowtie2 -x ${IN_DIR}/${file2} -1 ${READS_DIR}/${file3}_R1_paired.fq.gz -2 ${READS_DIR}/${file3}_R2_paired.fq.gz -S ${OUT_DIR}/${file2}-gtdb.sam

#generate a bam file-
samtools sort -n -o ${OUT_DIR}/${file2}-gtdb-sorted.bam ${OUT_DIR}/${file2}-gtdb.sam ##these are all the reads (to get all counts stats)

samtools view -b -F 4 ${OUT_DIR}/${file2}-gtdb-sorted.bam > ${OUT_DIR}/${file2}-prok.bam ##theses are prokaryotic reads  (to get all mapped to prok stats)

#samtools fixmate -m ${OUT_DIR}/${file2}-gtdb-sorted.bam ${OUT_DIR}/${file2}-gtdb-fixmate.bam
##-m tag is used by markdup

# Markdup needs position order
#samtools sort -o ${OUT_DIR}/${file2}-gtdb-positionsort.bam ${OUT_DIR}/${file2}-gtdb-fixmate.bam

# Finally mark duplicates
#samtools markdup -rsS ${OUT_DIR}/${file2}-gtdb-positionsort.bam ${OUT_DIR}/${file2}-gtdb-markdup.bam


##r- remove duplicates, s-print basic stats, S-remove the supplementary reads of duplcates as duplicates
```

## generate bed file from gtf file, get genelengths, get multicov file-

```
cd ${OUT_DIR}/
files=(*-gtdb-markdup.bam)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf$
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/gtdb-/}
file3=${file2/-gtdb-positionsort.bam/}

module load bowtie2/2.2.6
module load samtools/1.3.1
module load bedtools/v2.25.0


##convert the gtf to bed file- <perform this step interactively>
gff2bed convert2bed < seed_final2.gtf > seed_final2.bed
sort -k 1,1 -k2,2n seed_final2.bed > sorted_seed_final2.bed #create_bed_file.sh

#index the bam file and generate the bedtools map file-
samtools index ${OUT_DIR}/${file1}
bedtools multicov -bams ${OUT_DIR}/${file1} -bed ${IN_DIR}/sorted_${file3}_prokka.gff.gtf.bed > ${OUT_DIR}/multicov_output_${file3}.txt

#get genelengths- <perform this step interactively>
awk -F"\t" '{print $10"\t"$3-$2}' 229-01-prokka.map.gtf.bed | sed 's/gene_id //g' > genelength-229-01-prokka.map.gtf.bed.txt

#manipulate multicov file- <perform this step interactively>
cat ${IN_DIR}/multicov_output_seed-${num} | sed 's/gene_id//g' | awk '{print $1,$10,$11}'| tr ' ' '\t' > ${OUT_DIR}/multicov_output_seed_count-${num}.txt

#Reason why multicov can be used-The multicov tool excludes duplicate alignments unless you use the -D option and it excludes failed QC reads unless you use the -F option.

```

## generate TPM for each protein- Save the R script as "tpm_cal.R". To run the R script on multiple files on linux use the following code-

`module load R/3.6.1`

`Rscript tpm_cal.R genelength-229-01-prokka.map.gtf.bed.txt multicov_output_seed_count-229-01.txt TPM-229-01.txt`

``` 
#{save as tpm_cal.R}
library(dplyr)
library(plyr)
library(tidyverse)

#USAGE- Rscript [gene_length] [multicov_output] [output_file]

#args <- commandArgs(trailingOnly = T) #to run Rscript from command line
args <- commandArgs(TRUE)
gene_length <- args[1]
multicov_output <- args[2]
output_file <- args[3]



length<-read.table(gene_length)
colnames(length)<-c("gene_name", "gene_length")
count<-read.table(multicov_output)
count<-count%>%select(V1,V11,V12)
colnames(count)<-c("contig_name","gene_name","gene_count")

#merge the two files-
length_count<-merge(length, count, by.x="gene_name", by.y="gene_name")

length_count2<-length_count%>%mutate(countbylength=(length_count$gene_count*250)/length_count$gene_length)
sum_countbylength<-sum(length_count2$countbylength)
length_count2<-length_count2%>%mutate(TPM=(length_count2$countbylength*1000000)/(sum_countbylength))

write.table(length_count2, file=output_file, col.names=TRUE, sep="\t")

#formula of TPM used-> (read_count of each gene * library read_length * 10^6)/(gene_length of each gene * sum[read_count*library_read_length/gene_length])
```

## to get mapped read counts and TPM for reads mapped to prokaryotic contigs only-

`Rscript prokTPM_cal.R TPM-229-01.txt prokcontignames-229-01.txt prokTPM-229-01.txt`

```
#{save as prokTPM_cal.R}
library(dplyr)
library(stringr)

#USAGE- Rscript [orginal_tpm] [prokcontigs] [output_file]

#args <- commandArgs(trailingOnly = T) #to run Rscript from command line
args <- commandArgs(TRUE)
original_tpm <- args[1]
prokcontigs <- args[2]
output_file <- args[3]

original_tpm<-read.csv(file=original_tpm,sep="\t")
original_tpm$fullnames<-paste(original_tpm$file_name,original_tpm$contig_name,sep="_")

prokcontigs<-read.csv(file=prokcontigs,header=FALSE)

gtdb_tpm<-original_tpm%>%filter(fullnames %in% prokcontigs$V1) #get only prokcontig proteins out

gtdb_tpm2<-gtdb_tpm%>%filter(!str_detect(gene_name, "_gene"))
gtdb_tpm3<-gtdb_tpm2%>%filter(!str_detect(gene_name, "_mRNA"))
gtdb_tpm3$TPM<-NULL #remove the previous TPM values

sum_countbylength<-sum(gtdb_tpm3$countbylength)
gtdb_tpm3<-gtdb_tpm3%>%mutate(prokTPM=(gtdb_tpm3$countbylength*1000000)/(sum_countbylength))

write.table(gtdb_tpm3, file=output_file, col.names=TRUE, sep="\t")

#formula of TPM used-> (read_count of each gene * library read_length * 10^6)/(gene_length of each gene * sum[read_count*library_read_length/gene_length])
```

## join protein annotations (hpc_protanno.md) with TPM files to get the final files for statistical analysis-

`module load R/3.6.1`


```
 #in R

tpm<-read.csv("all-samples-prokTPM.txt",sep="\t") #cat TPM*.txt > all-samples-prokTPM.txt
markers<-read.csv("all-40markers.txt")
colnames(markers)<-c("cogmarkers","X","samplename","genename","fullnames")
tpm2<-merge(tpm,markers,by.x=c("file_name","gene_name"),by.y=c("samplename","genename"),all.x=TRUE) #merge with COG markers annotation file

cazy<-read.csv("all-cazymes.txt",header=FALSE)
colnames(cazy)<-c("file_name","cazy_anno","gene_name")
tpm3<-merge(tpm2,cazy,by=c("file_name","gene_name"),all.x=TRUE) #merge with cazyme annotation file



```


