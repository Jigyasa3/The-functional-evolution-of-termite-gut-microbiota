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

tpm<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_march2021/all-samples-salmon-tpm.txt",header=TRUE,sep="\t") #from salmon analysis
colnames(tpm)[which(names(tpm) == "Name")] <- "fullproteinnames"

cazymes<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/allsamples-evalue30-cazymes.txt",header=FALSE,sep="\t")
colnames(cazymes)<-c("contig_name","proteins","annotation","evalue","bitscore","samples")
cazymes$fullproteinnames<-paste(cazymes$samples,cazymes$proteins,sep="_")
cazymes<-cazymes%>%select(samples,proteins,annotation,fullproteinnames)
write.csv(cazymes,file="2-allsamples-evalue30-cazymes.txt")

#kofam<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/final-all-metagenomes-multiplekofamids-annotation.txt") #for the previous version
kofam<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_april2021/final-all-metagenomes-individualkofamids-annotation.txt") #the newest version
colnames(kofam)<-c("X","samples","proteins","annotation")
kofam$fullproteinnames<-paste(kofam$samples,kofam$proteins,sep="_")
kofam<-kofam%>%select(samples,proteins,annotation,fullproteinnames)
write.csv(kofam,file="2-final-all-metagenomes-multiplekofamids-annotation.txt")

pfam<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/3columns-allsamples-multiplepfamids-evalue30-pfamoutput.txt")
colnames(pfam)<-c("X","proteins","samples","annotation")
pfam<-pfam%>%select(samples,proteins,annotation)
pfam$fullproteinnames<-paste(pfam$samples,pfam$proteins,sep="_")
write.csv(pfam,file="2-3columns-allsamples-multiplepfamids-evalue30-pfamoutput.txt")

allannotations<-rbind(cazymes,kofam,pfam)%>%as.data.frame()
allannotations_tpm<-merge(allannotations,tpm,by.x=c("fullproteinnames"),by.y=c("fullproteinnames"))
write.csv(allannotations_tpm,file="allannotations_tpm.csv")

contig_length<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_march2021/all-contigs-readcount-jan2021.csv")

allannotations_tpm_contiglength<-merge(allannotations_tpm,contig_length,by.x=c("samples","proteins"),by.y=c("sample_names","gene_names"))
allannotations_tpm_contiglength$full_contig_name2<-paste(allannotations_tpm_contiglength$samples,allannotations_tpm_contiglength$contig_names,sep="_")
write.csv(allannotations_tpm_contiglength,file="allannotations_tpm_contiglength.csv")

taxonomy<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_march2021/allsamples-taxa-gtdb-lca.txt",header=FALSE) #from GTDB lca
library(forcats) #https://stackoverflow.com/questions/39126537/replace-na-in-a-factor-column
taxonomy$V7<-fct_explicit_na(taxonomy$V7, "0")

taxonomy2<-taxonomy%>%mutate(taxa=ifelse(V7==0,as.character(V5),as.character(V7)))

allannotations_tpm_contiglength_taxonomy<-merge(allannotations_tpm_contiglength,taxonomy2,by.x="full_contig_name2",by.y="V3") #merge by samplenames_contigIDs
write.csv(allannotations_tpm_contiglength_taxonomy,file="allannotations_tpm_contiglength_taxonomy.csv")

```

## c) Generating files for statistical analysis

```
allannotations_tpm_contiglength_taxonomy<-read.csv("allannotations_tpm_contiglength_taxonomy.csv")
allannotations_tpm_contiglength_taxonomy$contig_length<-as.numeric(as.character(allannotations_tpm_contiglength_taxonomy$contig_length))
allannotations_tpm_contiglength_taxonomy$TPM<-as.numeric(as.character(allannotations_tpm_contiglength_taxonomy$TPM))


joined1<-allannotations_tpm_contiglength_taxonomy%>%filter(contig_length>=1000 & TPM>=1)
write.csv(joined1,file="all-metagenome-annotation-taxonomy-1000bps.txt")


## NOTE- for statistical analysis ensure that there are no duplicates per proteinID otherwise they will counted in summing.
library(data.table)
setDT(joined1)
joined1_genes<-joined1[, .(TPM = sum(TPM)), by = .(samples,annotation)]
write.csv(joined1_genes,file="1000bps_genes.txt")

joined3<-allannotations_tpm_contiglength_taxonomy%>%filter(contig_length>=5000 & TPM>=1)
write.csv(joined3,file="all-metagenome-annotation-taxonomy-5000bps.txt")

setDT(joined3)
joined3_genes<-joined3[, .(TPM = sum(TPM)), by = .(samples,annotation)]
write.csv(joined3_genes,file="5000bps_genes.txt")
```

### If there is "\<NA>" in a column with factor variables in R-
```
library(forcats) #https://stackoverflow.com/questions/39126537/replace-na-in-a-factor-column
test$V7<-fct_explicit_na(test$V7, "0")

test2<-test%>%mutate(taxa=ifelse(V7==0,as.character(V5),as.character(V7)))
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
