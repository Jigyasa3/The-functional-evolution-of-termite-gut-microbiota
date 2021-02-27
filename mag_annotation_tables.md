## extracting lignocellulose degrading genes (KEGG annotation) from MAGs

```
module load R/3.6.1

head(reductiveacetogenesis_genes.txt) #this is a comma separated file with KEGG IDs for a pathways
K00198,acsA_ec1.2.7.4
K14138,acsB_ec2.3.1.169
K00194,acsD_ec2.1.1.245
K00197,acsC_ec2.1.1.245
K15023,acsE_ec2.1.1.258
K01938,fhs_ec6.3.4.3
K01491,folD_ec1.5.1.5_ec3.5.4.9
K13990,ftcD_ec4.3.1.4
K00297,metF_ec1.5.1.20
K00625,pta_ec2.3.1.8
K00925,ack_ec2.7.2.1
K00122,fdhF1_ec1.17.1.9

Rscript pathway_extract_kegg.R reductiveacetogenesis_genes.txt allkofam_bins.txt reductiveacetogenesis_allbins.txt
```

```
#save as "pathway_extract_kegg.R"
#USAGE- Rscript [methano] [kofam] [output_file1]

args <- commandArgs(TRUE)
methano <- args[1] #comma separated file
kofam <- args[2]
output_file1 <- args[3]



methano<-read.csv(file=methano,header=FALSE)
file<-read.csv(file=kofam,header=TRUE,sep="\t")

joined<-merge(file,methano,by.x="koid",by.y="V1")
write.csv(joined,file=output_file1)
```
## converting gene presence to binary (0/1)

```
head(659bins-gtdb-taxonomy_fullnames.txt) #taxonomy file generated from gtdb-tk taxonomic annotation of bins
"","user_genome","classification","completeness","contamination","gc","N50","size"
"1","229-01.bin.107","Archaea_Thermoplasmatota_Thermoplasmata_Methanomassiliicoccales_Methanomethylophilaceae_Methanoplasma_s__",71.43,1.612,0.502,2180,1032627
"2","229-01.bin.118","Bacteria_Bacteroidota_Bacteroidia_Bacteroidales_P3__s__",44.41,0.268,0.389,1764,916589
"3","229-01.bin.138","Bacteria_Actinobacteriota_Actinomycetia_Actinomycetales_Bifidobacteriaceae_Ancillula_s__",30.95,1.333,0.37,1371,544050
"4","229-01.bin.23","Bacteria_Proteobacteria_Gammaproteobacteria_Burkholderiales_Rhodocyclaceae_Azovibrio_s__",54.96,1.664,0.582,1628,1247006

Rscript bins_spread_binary.R reductiveacetogenesis_allbins.txt 659bins-gtdb-taxonomy_fullnames.txt reductiveacetogenesis_allbins_spreadbinary.txt

```

```
#save as "bins_spread_binary.R"
#USAGE- Rscript [methano] [taxonomy] [output_file]

args <- commandArgs(TRUE)
methano <- args[1] #comma separated output file from "pathway_extract_kegg.R"
taxonomy <- args[2]
output_file <- args[3]

library(dplyr)
library(tidyr)

#read the output from "pathway_extract_kegg.R" or "pathway_extract_pfam.R"
methano<-read.csv(file=methano)

#get taxonomy-
taxonomy<-read.csv(file=taxonomy)


#convert koid to presence/absence + spread-
joined_test<-methano%>%select(file_name,koid) ##imp. to extract only one distinct column to spread against.
methano_spread<-joined_test%>% mutate(yesno = 1) %>% distinct %>% spread(koid, yesno, fill = 0)

#get the taxonomy column back into final file-
methano_spread2<-merge(methano_spread,taxonomy,by.x="file_name",by.y="user_genome")

write.csv(methano_spread2,file=output_file)

```


## BLASTp analysis against ANNOTREE database

```

```


## Generating the final BLAST results file

```
Rscript blast_outputs.R 
```


```
#save as "blast_outputs.R"
#USAGE- Rscript [blast] [outputfile1] [outputfile2]

args <- commandArgs(TRUE)
blast <- args[1] #comma separated file
outputfile1 <- args[2] #above60perc100aa-blastresult.txt
outputfile2 <- args[3] #notkofam.txt


library(dplyr)
library(tidyr)

##perform the following for each file separately-

blast<-read.csv(file=blast,header=FALSE,sep="\t")
blast$file_name<-gsub("_.*$","",blast$V1)
blast$protein1<-unlist(lapply(strsplit(as.character(blast$V1),split="_"),"[",2))
blast$protein2<-unlist(lapply(strsplit(as.character(blast$V1),split="_"),"[",3))
blast$protein_name<-paste(blast$protein1,blast$protein2,sep="_")

blast$fullname<-paste(blast$file_name,blast$protein_name,sep="__")

###fullname column "229-01.bin.107__GLHJBGCB_00001" is binname and protein-name

#get the best match for each protein sequence with percent identity > 60% and sequence length > 100bps
blast2<-blast%>%group_by(fullname)%>%filter(V3 >=60 & V4 >=100)%>%as.data.frame()
blast3<-blast2%>%group_by(fullname)%>%summarize(max_V4=max(V4,na.rm=TRUE)) #get the best match


write.csv(blast2,file=outputfile1)

```
## Adding BLAST results to KEGG results to generate the final table

```
##bacterial genes-
kofam<-read.csv("wood-ljungdahl-ikedaohtsubo-bacterial-kofamoutput-spreadbinary.txt")
gtdb<-read.csv("notkofam-above60-bacterialwood-spreadbinary.txt",sep="\t") #presence shown as 1*


joined2<-merge(kofam,gtdb,by=c("file_name","classification","completeness","contamination","gc","N50","size"),all=TRUE)

joined2[is.na(joined2)] <- 0

#create a final column per KOID with presence of gene using kofamdb indicated by "1", blast by "1*" nad absence by "0"-

joined2<-joined2%>%mutate(K00122=ifelse(K00122.x=="1" & K00122.y=="1*","1",ifelse(K00122.x=="0" & K00122.y=="1*","1*",ifelse(K00122.x=="1" & K00122.y=="0","1","0"))))
joined2$K00122.x<-NULL
joined2$K00122.y<-NULL

joined2<-joined2%>%mutate(K00194=ifelse(K00194.x=="1" & K00194.y=="1*","1",ifelse(K00194.x=="0" & K00194.y=="1*","1*",ifelse(K00194.x=="1" & K00194.y=="0","1","0"))))
joined2$K00194.x<-NULL
joined2$K00194.y<-NULL

joined2<-joined2%>%mutate(K00197=ifelse(K00197.x=="1" & K00197.y=="1*","1",ifelse(K00197.x=="0" & K00197.y=="1*","1*",ifelse(K00197.x=="1" & K00197.y=="0","1","0"))))
joined2$K00197.x<-NULL
joined2$K00197.y<-NULL

joined2<-joined2%>%mutate(K00198=ifelse(K00198.x=="1" & K00198.y=="1*","1",ifelse(K00198.x=="0" & K00198.y=="1*","1*",ifelse(K00198.x=="1" & K00198.y=="0","1","0"))))
joined2$K00198.x<-NULL
joined2$K00198.y<-NULL

joined2<-joined2%>%mutate(K00297=ifelse(K00297.x=="1" & K00297.y=="1*","1",ifelse(K00297.x=="0" & K00297.y=="1*","1*",ifelse(K00297.x=="1" & K00297.y=="0","1","0"))))
joined2$K00297.x<-NULL
joined2$K00297.y<-NULL

joined2<-joined2%>%mutate(K00625=ifelse(K00625.x=="1" & K00625.y=="1*","1",ifelse(K00625.x=="0" & K00625.y=="1*","1*",ifelse(K00625.x=="1" & K00625.y=="0","1","0"))))
joined2$K00625.x<-NULL
joined2$K00625.y<-NULL

joined2<-joined2%>%mutate(K00925=ifelse(K00925.x=="1" & K00925.y=="1*","1",ifelse(K00925.x=="0" & K00925.y=="1*","1*",ifelse(K00925.x=="1" & K00925.y=="0","1","0"))))
joined2$K00925.x<-NULL
joined2$K00925.y<-NULL

joined2<-joined2%>%mutate(K01491=ifelse(K01491.x=="1" & K01491.y=="1*","1",ifelse(K01491.x=="0" & K01491.y=="1*","1*",ifelse(K01491.x=="1" & K01491.y=="0","1","0"))))
joined2$K01491.x<-NULL
joined2$K01491.y<-NULL

joined2<-joined2%>%mutate(K01938=ifelse(K01938.x=="1" & K01938.y=="1*","1",ifelse(K01938.x=="0" & K01938.y=="1*","1*",ifelse(K01938.x=="1" & K01938.y=="0","1","0"))))
joined2$K01938.x<-NULL
joined2$K01938.y<-NULL

joined2<-joined2%>%mutate(K13990=ifelse(K13990.x=="1" & K13990.y=="1*","1",ifelse(K13990.x=="0" & K13990.y=="1*","1*",ifelse(K13990.x=="1" & K13990.y=="0","1","0"))))
joined2$K13990.x<-NULL
joined2$K13990.y<-NULL

joined2<-joined2%>%mutate(K14138=ifelse(K14138.x=="1" & K14138.y=="1*","1",ifelse(K14138.x=="0" & K14138.y=="1*","1*",ifelse(K14138.x=="1" & K14138.y=="0","1","0"))))
joined2$K14138.x<-NULL
joined2$K14138.y<-NULL

joined2<-joined2%>%mutate(K15023=ifelse(K15023.x=="1" & K15023.y=="1*","1",ifelse(K15023.x=="0" & K15023.y=="1*","1*",ifelse(K15023.x=="1" & K15023.y=="0","1","0"))))
joined2$K15023.x<-NULL
joined2$K15023.y<-NULL

joined2$domain<-gsub("_.*$","",joined2$classification)
joined2<-joined2%>%filter(domain != "Archaea") #removing archaeal MAGs, as these are bacterial genes.

write.csv(joined2,file="kofam-gtdb-ikedaohtsubo-bacterial-spreadbinary.txt")

#------------------------------------------------------------------------------------------
##CHANGING THE ORDER OF GENE COLUMNS in MAGs-

bacteria_reductive_mag<-read.csv("kofam-gtdb-ikedaohtsubo-bacterial-spreadbinary.txt",row.names = 1)
rownames(bacteria_reductive_mag)<-bacteria_reductive_mag$file_name
bacteria_reductive_mag2<-bacteria_reductive_mag%>%select(-c(K00193,classification,completeness,contamination,gc,N50,size,domain))

genenames_bacteria<-c("K01938","K01491","K13990","K00297","K00198","K14138","K00194","K00197","K15023","K00625","K00925","K00122") #geneorder
genenames_bacteria[!(genenames_bacteria) %in% colnames(bacteria_reductive_mag2)]
bacteria_reductive_mag_ordered<-bacteria_reductive_mag2[,genenames_bacteria] #order the metagenome contigs data
bacteria_reductive_mag_ordered$file_name<-rownames(bacteria_reductive_mag_ordered)
mag_data<-bacteria_reductive_mag%>%select(file_name,classification,completeness,contamination,gc,N50,size,domain)

bacteria_reductive_mag_ordered2<-merge(bacteria_reductive_mag_ordered,mag_data,by="file_name")
write.csv(bacteria_reductive_mag_ordered2,file="kofam-gtdb-ikedaohtsubo-bacterial-spreadbinary-ordered.csv")

```
