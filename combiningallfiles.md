## All files (from hpc_functionalannotation.md , hpc_taxonomy.md and hpc_tpmcal.md) were joined in R-

### a) Generating the final COG taxonomy file-
```
## in R
module load R/3.6.1

tpm<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/all-samples-prokTPM.txt",header=TRUE,sep="\t")
tpm$fullproteinnames<-paste(tpm$file_name,tpm$gene_name,sep="_")

cogs<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/allcogs-allsamples-finalkrakenoutput.csv")
colnames(cogs)<-c("markers","X","fullproteinnames","markertaxonomy")

tpm_cogs<-merge(tpm,cogs,by="fullproteinnames")
write.csv(tpm_cogs,file="tpm_cogs_allsamples.csv")

```

### b) Generating the final functional annotation file-

```
library(dplyr)

tpm<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/all-samples-prokTPM.txt",header=TRUE,sep="\t")
tpm$fullproteinnames<-paste(tpm$file_name,tpm$gene_name,sep="_")

cazymes<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/allsamples-evalue30-cazymes.txt",,header=FALSE,sep="\t")
colnames(cazymes)<-c("contig_name","proteins","annotation","evalue","bitscore","samples")
cazymes$fullproteinnames<-paste(cazymes$samples,cazymes$proteins,sep="_")
cazymes<-cazymes%>%select(samples,proteins,annotation,fullproteinnames)

kofam<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/final-all-metagenomes-multiplekofamids-annotation.txt")
colnames(kofam)<-c("X","samples","proteins","annotation")
kofam$fullproteinnames<-paste(kofam$samples,kofam$proteins,sep="_")
kofam<-kofam%>%select(samples,proteins,annotation,fullproteinnames)

pfam<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/3columns-allsamples-multiplepfamids-evalue30-pfamoutput.txt")
colnames(pfam)<-c("X","proteins","samples","annotation")
pfam<-pfam%>%select(samples,proteins,annotation)
pfam$fullproteinnames<-paste(pfam$samples,pfam$proteins,sep="_")


allannotations<-rbind(cazymes,kofam,pfam)%>%as.data.frame()
allannotations_tpm<-merge(allannotations,tpm,by.x=c("fullproteinnames","proteins","samples"),by.y=c("fullproteinnames","gene_name","file_name"))
write.csv(allannotations_tpm,file="allannotations_tpm.csv")

taxonomy<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/allsamples-krakentaxonomy-feb2021.txt")

allannotations_tpm_taxonomy<-merge(allannotations_tpm,taxonomy,by.x="fullnames",by.y="V1.x") #merge by contigIDs and samplenames
write.csv(allannotations_tpm_taxonomy,file="allannotations_tpm_taxonomy.csv")

```

## c) Generating files for statistical analysis

```


```
