## In this version, kofamscan was run with output format "mapper" which didnot give the threshold value or evalue per KEGG annotation. So, all the geneIDs with duplicates were removed.

```

#extract proteinIDs annotated with KEGGid-
for i in *-prokka.faa-kopfam;do awk '$2 ~ /^K/ { print $0 }' ${i} > annotated-${i};done #extract rows that starts with "K" in the second column


#remove duplicates-
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
