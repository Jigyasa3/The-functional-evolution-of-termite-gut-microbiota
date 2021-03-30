## TPM calculation of relative abundance accounting for GC bias also

### Use https://github.com/Jigyasa3/termite_146guts_microbes_function/blob/main/hpc_taxonomy.md to get the taxonomy of each contig. Extract the microbial contigs and the protein coding nucleotide sequences in microbial contigs-
```
##get microbial contigs out-
grep "d__" taxa-gtdb-lca-method2-matches-filename-230-13.fasta.txt > microbesonly-taxa-gtdb-lca-method2-matches-filename-230-13.fasta.txt
##get microbial fna sequences out-
awk -F"," '{print $3}' microbesonly-taxa-gtdb-lca-method2-matches-filename-230-13.fasta.txt | sed 's/"//g' > microbesonly-contigheaders-230-13.txt
#-----------------------------------------------------------------------------------------------------------------------
##run microbes_fna.R to get the proteinIDs corresponding to contignames

module load R/3.6.1
R
#load files-
gtf<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/all_functions_all_taxonomy/gtf_files_Dec2019/named-gtffiles/filename-230-13-prokka.map.gtf",header=FALSE,sep="\t")
gtf$V1<-gsub("-prokka.map.gtf","",gtf$V1)
gtf$names<-paste(gtf$V1,gtf$V2,sep="_")
gtf$V10<-gsub("gene_id ","",gtf$V10)
gtf$proteinnames<-paste(gtf$V1,gtf$V10,sep="_")


microbescontigs<-read.csv("microbesonly-contigheaders-230-13.txt",header=FALSE)

library(dplyr)
extracted<-gtf%>%filter(names %in% microbescontigs$V1)
write.csv(extracted,file="gtdbmicrobial-230-13-prokka.map.gtf")

#-----------------------------------------------------------------------------------------------------------------------
##get protein coding nucleotide sequences from each contig-
awk -F"," '{print $13}' gtdbmicrobial-230-13-prokka.map.gtf | sed 's/"//g' > proteinnames-230-13.txt
seqtk subseq 2-renamed-230-13-prokka.fasta ${IN_DIR}/proteinnames-230-13.txt > microbial_fna/gtdbmicrobial-230-13-prokka.fasta


```


### Run Salmon to get GC corrected TPM values for each protein coding nucleotide sequences-
```
module load salmon/1.4.0
salmon index -p 30 -t ${FNA_DIR}/gtdbmicrobial-230-13-prokka.fasta -i ${OUT_DIR}/gtdbmicrobial-230-13-prokka_index

salmon quant -i ${OUT_DIR}/gtdbmicrobial-230-13-prokka_index --libType IU -1 ${READS_DIR}/230-13_R1_paired.fq.gz -2 ${READS_DIR}/230-13_R2_paired.fq.gz -o ${OUT_DIR}/230-13.quant --meta -p 30
#the salmon scripts and the parameters are taken from https://github.com/bxlab/metaWRAP/blob/master/bin/metawrap-modules/quant_bins.sh


```
