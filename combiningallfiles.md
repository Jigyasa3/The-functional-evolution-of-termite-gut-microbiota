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

contigs<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_2020/all-contigs-readcount.csv")
contigs<-contigs%>%select(contig_names,sample_names,contig_length)%>%unique()

cogs_contigs<-merge(cogs,contigs,by.x=c("file_name","contig_names"),by.y=c("sample_names","contig_names"))
cogs_contigs$contig_length<-as.numeric(as.character(cogs_contigs$contig_length))
cogs_contigs_1000above<-cogs_contigs%>%filter(contig_length>=1000) #get cogs present in contigs >1000bps.
write.csv(cogs_contigs_1000above,file="tpm_cogs_allsamples_feb2021_1000bpscontigs.csv")

cogs_contigs_1000above2<-cogs_contigs%>%filter(contig_length>=1000 & read_count >=100 & TPM>=1) #get cogs present in contigs >1000bps and all the same filters applied to functionally annotated genes.
write.csv(cogs_contigs_1000above2,file="tpm_cogs_allsamples_feb2021_1000bpscontigs_100counts1tpm.csv")

##NOTE-Use "tpm_cogs_allsamples_feb2021_1000bpscontigs_100counts1pm_similartocontigs.csv" file for all statistical analysis. See "hpc_taxonomy.md" file for details.
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
allannotations_tpm_taxonomy<-read.csv("allannotations_tpm_taxonomy.csv")
allannotations_tpm_taxonomy$gene_count<-as.numeric(as.character(allannotations_tpm_taxonomy$gene_count))
allannotations_tpm_taxonomy$prokTPM<-as.numeric(as.character(allannotations_tpm_taxonomy$prokTPM))

contig_length<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_2020/all-contigs-readcount.csv")

allannotations_tpm_taxonomy_contiglen<-merge(allannotations_tpm_taxonomy,contig_length,by.x=c("samples","contig_name"),by.y=c("sample_names","contig_names"))
write.csv(allannotations_tpm_taxonomy_contiglen,file="allannotations_tpm_taxonomy_contiglen.csv")

joined1<-allannotations_tpm_taxonomy_contiglen%>%filter(contig_length>=1000 & read_count >=100 & TPM>=1)
write.csv(joined1,file="all-metagenome-annotation-taxonomy-1000bps.txt")


## NOTE- for statistical analysis ensure that there are no duplicates per proteinID otherwise they will counted in summing.
setDT(joined1)
joined1_genes<-joined1[, .(TPM = sum(prokTPM)), by = .(samplename,annotation)]
write.csv(joined1_genes,file="1000bps_genes.txt")

joined3<-allannotations_tpm_taxonomy_contiglen%>%filter(contig_length>=5000 & read_count >=100 & TPM>=1)
write.csv(joined3,file="all-metagenome-annotation-taxonomy-5000bps.txt")

setDT(joined3)
joined3_genes<-joined3[, .(TPM = sum(prokTPM)), by = .(samplename,annotation)]
write.csv(joined3_genes,file="5000bps_genes.txt")
```

## c.1) Adding prokTPM to NiFe and FeFe hydrogenases catalytic subunits annotated by Hyddb-

```
joined1<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/all-metagenome-annotation-taxonomy-1000bps_feb2021.txt")

#FeFe hydrogenases
fefe<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/bin_taxonomy/all_my_bins/prokka_output/all_faa/pathways/hydrogenases/metagenome_contigs/evalue30_cutoff/2-fefe-hyddb-results.csv",header=FALSE,sep=",")

allmeta_fefe<-merge(fefe,allmeta,by.x="V1",by.y="fullproteinnames")
allmeta_fefe<-allmeta_fefe%>%select(V1,V2,samples,contig_name,proteins,gene_length,gene_count,prokTPM,full_contig_name,contig_length)%>%unique()
write.csv(allmeta_fefe,file="fefe_hyddb_prokTPM.csv")

library(data.table)
setDT(allmeta_fefe)
allmeta_fefe2<-allmeta_fefe[, .(TPM = sum(prokTPM)), by = .(samples,V2)] 
write.csv(allmeta_fefe2,file="fefe_hydd_groups_prokTPM.csv")

#NiFe hydrogenases
nife<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/bin_taxonomy/all_my_bins/prokka_output/all_faa/pathways/hydrogenases/metagenome_contigs/evalue30_cutoff/2-nife-hyddb-results.csv",header=FALSE,sep=",")

allmeta_nife<-merge(nife,allmeta,by.x="V1",by.y="fullproteinnames")
allmeta_nife<-allmeta_nife%>%select(V1,V2,samples,contig_name,proteins,gene_length,gene_count,prokTPM,full_contig_name,contig_length)%>%unique()
write.csv(allmeta_nife,file="nife_hyddb_prokTPM.csv")

library(data.table)
setDT(allmeta_nife)
allmeta_nife2<-allmeta_nife[, .(TPM = sum(prokTPM)), by = .(samples,V2)] 
write.csv(allmeta_nife2,file="nife_hydd_groups_prokTPM.csv")

```
