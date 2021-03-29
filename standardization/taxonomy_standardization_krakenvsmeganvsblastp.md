## Compare different methods of marker-gene taxonomy annotation- KRAKEN2, manual Blastp, MEGAN using GTDB database and DIAMOND lca (outmft=102)

### A) MEGAN analysis- Which filters to use for contig annotation (blastx based)
```
#the filters used in blast2rma can give different results to the same contigs-


filter1<-read.csv("meganoutput-nr-matches-230-12-contigs.txt",header=FALSE,sep="\t")
##filter used: --maxExpected 1e-15 --minPercentIdentity 50 --lcaCoveragePercent 60

filter2<-read.csv("meganoutput-newfilters-nr-matches-230-12-contigs.txt",header=FALSE,sep="\t")
##filters used: --maxExpected 1e-15

joined1<-merge(filter1,filter2,by="V1")
write.csv(joined1,"meganoutput-230-12-twodiff_filters.csv")

filter3<-read.csv("meganoutput-newfilters2-nr-matches-230-12-contigs.txt",header=FALSE,sep="\t")
##filters used: --minScore 140 --minReadLength 50 --maxExpected 1e-25 --minPercentIdentity 65

joined2<-merge(joined1,filter3,by="V1",all.x=TRUE)
write.csv(joined2,file="meganoutput-230-12-threediff_filters.csv")

#Q/ how many contig level annotation match between filter1 and filter2?
9 contig annotations dont match out of 36 contigs. Out of these 9, 6 annotations vary at phyla level.

#Q/ how many contig level annotation match between filter1, filter2 and filter3?
11 contig annotations donot match out of 36 contigs, the mismatch is at phyla level or there is no-match for that contig in  filter3 file. 5 contig annotations out of 11 vary between filter1 and filter2 also i.e. not supported in any filtering criteria.
8 contig annotations vary at taxonomic positions-at family level.

NOTE-Filter 3 will be used as it is the most stringent of all three thereby removing any false positive results.
```

### B) Does MEGAN blastp of each protein match with MEGAN blastx of all proteins in a contig and manual blastp results?
```
##protein MEGAN script-
diamond blastp --db ${DB_DIR}/nr.dmnd --query ${IN_DIR}/${file1} --outfmt 100 --out ${OUT_DIR}/nr-nolongreads-matches-${file1} --threads 15
#convert daa to rma file-
/home/j/jigyasa-arora/local/megan/tools/blast2rma --in ${file1} --format DAA --blastMode BlastP --out ${file1}.rma --minScore 140 --minReadLength 50 --maxExpected 1e-25 --minPercentIdentity 65 --lcaAlgorithm naive  --mapDB ${DB_DIR}/megan-map-Jan2021.db --threads 16 --verbose
#get the final results in tab delimited output-
/home/j/jigyasa-arora/local/megan/tools/rma2info --in ${file1} --read2class GTDB --paths --out meganoutput-${file1}.txt

##contig MEGAN script-
diamond blastx --db ${DB_DIR}/nr.dmnd --query ${IN_DIR}/${file1} --outfmt 100 --out ${OUT_DIR}/nr-matches-${file1} --threads 15 --long-reads
#convert daa to rma file-
/home/j/jigyasa-arora/local/megan/tools/blast2rma --in ${OUT_DIR}/${file1} --format DAA --blastMode BlastX --out ${OUT_DIR}/${file1}.rma --longReads --minScore 140 --minReadLength 50 --maxExpected 1e-25 --minPercentIdentity 65 --lcaAlgorithm longReads  --mapDB ${DB_DIR}/megan-map-Jan2021.db --threads 16 --verbose
#get the final results in tab delimited output-
/home/j/jigyasa-arora/local/megan/tools/rma2info --in ${file1} --read2class GTDB --paths --out meganoutput-${file1}.txt

##blastp script-
diamond blastp --db ${DB_DIR}/gtdb_decipher_r89.faa.dmnd --query ${IN_DIR}/${file1} --outfmt 6 --out ${OUT_DIR}/gtdb-ver89-matches-${file1}.txt --threads 15
#apply filter-
awk -F"\t" '$4/($8-$7)>0.75 && $4/($10-$9) > 0.75 && $11<1e-15 {print $0}' gtdb-ver89-matches-230-12-Microcerotermes-1000bpscontigs_100counts1tpm.txt > 75cov-evalue-15-gtdb-ver89-matches-230-12-Microcerotermes-1000bpscontigs_100counts1tpm.txt
#get the top match per protein ID-
for next in $(cut -f1 75cov-evalue-15-gtdb-ver89-matches-230-12-Microcerotermes-1000bpscontigs_100counts1tpm.txt | sort -u); do grep -w -m 1 "$next" 75cov-evalue-15-gtdb-ver89-matches-230-12-Microcerotermes-1000bpscontigs_100counts1tpm.txt; done > top-sorted-filteredvalue-markes-230-12.txt


##comparisons-
CONTIG2-
manual blastp of proteins in a contig show-
Bacteria;Spirochaetota;Spirochaetia;Treponematales;Treponemataceae_B;Treponema_E;Treponema_Eprimitia for 5 out of 5 proteins in contig2.
     
MEGAN annotation of each protein shows-
Bacteria;Spirochaetota;Spirochaetia for 2 out of 5 proteins
Bacteria;Spirochaetota;Spirochaetia;Treponematales;Treponemataceae_B for 3 out of 5 proteins.

MEGAN annotation of the whole contig2 (blastx) shows-
Bacteria;Spirochaetota;Spirochaetia;Treponematales;Treponemataceae_B;Treponema_E


CONTIG3-
manual blastp-
Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Anaerotignaceae for 3 out of 6 proteins.

MEGAN annotation of each protein-
Bacteria;Firmicutes_A;Clostridia;Lachnospirales (order level) for all 6 proteins

MEGAN annotation of the whole contig3-
Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospiraceae


CONTIG1-
manual blastp-
Bacteria;Fibrobacterota;Chitinivibrionia;Chitinivibrionales;Chitinispirillaceae;GUT77 for 8 out of 9 proteins.
The one protein that is annotated differently, is annotated as 
Bacteria;Fibrobacterota;Chitinivibrionia;Chitinivibrionales;Chitinispirillaceae;Chitinispirillum;Chitinispirillum alkaliphilum

MEGAN annotation of each protein-
Bacteria;Fibrobacterota;Chitinivibrionia;Chitinivibrionales;Chitinispirillaceae;Chitinispirillum;Chitinispirillum alkaliphilum;Chitinispirillum alkaliphilum ACht6-1 for 8 out of 9 proteins.
The one protein that is annotated differently, is not the same protein as manual blastp results and is just annotated as Bacteria (domain level).

MEGAN annotation of the whole contig1-
Bacteria;Fibrobacterota;Chitinivibrionia;Chitinivibrionales;Chitinispirillaceae;Chitinispirillum;Chitinispirillum alkaliphilum;Chitinispirillum alkaliphilum ACht6-1

#NOTE- MEGAN annotation of contig and proteins match each other, even though it is at a higher taxonomic level.
```

### C) KRAKEN results of each protein vs whole contig vs blastp results-
```
kraken<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/tpm_cogs_allsamples_feb2021_1000bpscontigs_100counts1tpm.csv")
kraken$phyla<-gsub("_c__.*$","",kraken$markertaxonomy)
kraken$class<-gsub("_o__.*$","",kraken$markertaxonomy)
kraken$order<-gsub("_f__.*$","",kraken$markertaxonomy)
kraken$family<-gsub("_g__.*$","",kraken$markertaxonomy)


blastp<-read.csv("top-sorted-filteredvalue-markes-230-12.txt",header=FALSE,sep="\t")
blastp$V2<-gsub(";","_",blastp$V2)
blastp$phyla<-gsub("_c__.*$","",blastp$V2)
blastp$class<-gsub("_o__.*$","",blastp$V2)
blastp$order<-gsub("_f__.*$","",blastp$V2)
blastp$family<-gsub("_g__.*$","",blastp$V2)

blastp$phyla<-gsub("^.*p__","p__",blastp$phyla)
blastp$class<-gsub("^.*p__","p__",blastp$class)
blastp$order<-gsub("^.*p__","p__",blastp$order)
blastp$family<-gsub("^.*p__","p__",blastp$family)

#join by protein-ids-
joined_all<-merge(kraken,blastp,by.x=c("fullproteinnames"),by.y=c("V1"))
write.csv(joined_all,file="kraken_blastp_230-12_cogs.csv")

nrow(kraken)
[1] 177

#join by taxonomy-
joined_phyla<-merge(kraken,blastp,by.x=c("fullproteinnames","phyla"),by.y=c("V1","phyla"))
nrow(joined_phyla)
[1] 95

joined_class<-merge(kraken,blastp,by.x=c("fullproteinnames","class"),by.y=c("V1","class"))
nrow(joined_class)
[1] 94
joined_order<-merge(kraken,blastp,by.x=c("fullproteinnames","order"),by.y=c("V1","order"))
nrow(joined_order)
[1] 83
joined_class<-merge(kraken,blastp,by.x=c("fullproteinnames","class"),by.y=c("V1","class"))
nrow(joined_class)
[1] 94
joined_family<-merge(kraken,blastp,by.x=c("fullproteinnames","family"),by.y=c("V1","family"))
nrow(joined_family)
[1] 38

#NOTE- ~50% of taxonomy match between manual Blastp and KRAKEN protein annotation.
#----------------------------------------------------------------------------------------------------------------------------
krakencontig<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/paper1/markergenes/markers-rpkm/individualanalysis_feb2021/230-12-Microcerotermes-krakencontigtaxonomy.txt",header=FALSE)
#generated from "tpm_cogs_allsamples_feb2021_1000bpscontigs.csv"

krakencontig$V16<-gsub("\\|","_",krakencontig$V16)
krakencontig$phyla<-gsub("_c__.*$","",krakencontig$V16)
krakencontig$class<-gsub("_o__.*$","",krakencontig$V16)
krakencontig$order<-gsub("_f__.*$","",krakencontig$V16)
krakencontig$family<-gsub("_g__.*$","",krakencontig$V16)

joined_contig_phyla<-merge(krakencontig,blastp,by.x=c("V5","phyla"),by.y=c("V1","phyla"))
nrow(joined_contig_phyla)
[1] 95

#NOTE- there is a difference in blastp annotation and KRAKEN2 annotation for marker genes and  contig annotation. Only 53% match with manual blastp annotation.

```

### d)DIAMOND lca 

```
##standardize the parameters for DIAMOND lca-

###CONTIGS-
#method1: --sensitive --query-cover 60 --id 65 --evalue 1e-24
awk '$2!=0 {print $0}' gtdb-lca-method1-matches-prot-all-COG0087.fasta.txt > output-gtdb-lca-method1-matches-prot-all-COG0087.fasta.txt

#method2: --query-cover 60
awk '$2!=0 {print $0}' gtdb-lca-method2-matches-prot-all-COG0087.fasta.txt > output-gtdb-lca-method2-matches-prot-all-COG0087.fasta.txt

#------------------------------------------------------------------------------------------------------------------------------------------------------
##creating the accession_to_taxid_taxonomy file in R-
#NOTE- The final output is "gtdb_ver95_alllca_taxid.csv" which is used in all downstream analysis. This code doesn't have to be run

module load R/3.6.1
R
alltaxonomy<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/work_backup_scripts_databases/230run_nomismatch_joinedfiles/230run_nomismatch_joinedfiles/taxonomy/gtdb_to_diamond/alltaxonomy_r95.tsv",header=FALSE,sep="\t")
taxids<-read.csv("/apps/unit/BioinfoUgrp/DB/diamondDB/GTDB/r95/taxID_info.tsv",header=TRUE,sep="\t")

#create taxonomy lca paths for each phylogenetic level-
alltaxonomy$domain<-unlist(lapply(strsplit(as.character(alltaxonomy$V2),split=";"),"[",1))
alltaxonomy$phyla<-unlist(lapply(strsplit(as.character(alltaxonomy$V2),split=";"),"[",2))
alltaxonomy$class<-unlist(lapply(strsplit(as.character(alltaxonomy$V2),split=";"),"[",3))
alltaxonomy$order<-unlist(lapply(strsplit(as.character(alltaxonomy$V2),split=";"),"[",4))
alltaxonomy$family<-unlist(lapply(strsplit(as.character(alltaxonomy$V2),split=";"),"[",5))
alltaxonomy$genus<-unlist(lapply(strsplit(as.character(alltaxonomy$V2),split=";"),"[",6))
alltaxonomy$species<-unlist(lapply(strsplit(as.character(alltaxonomy$V2),split=";"),"[",7))

#create lca-
alltaxonomy$phylalca<-paste(alltaxonomy$domain,alltaxonomy$phyla,sep="_")
alltaxonomy$orderlca<-paste(alltaxonomy$domain,alltaxonomy$phyla,alltaxonomy$class,alltaxonomy$order,sep="_")
alltaxonomy$familylca<-paste(alltaxonomy$domain,alltaxonomy$phyla,alltaxonomy$class,alltaxonomy$order,alltaxonomy$family,sep="_")
alltaxonomy$genuslca<-paste(alltaxonomy$domain,alltaxonomy$phyla,alltaxonomy$class,alltaxonomy$order,alltaxonomy$family,alltaxonomy$genus,sep="_")
alltaxonomy$specieslca<-paste(alltaxonomy$domain,alltaxonomy$phyla,alltaxonomy$class,alltaxonomy$order,alltaxonomy$family,alltaxonomy$genus,alltaxonomy$species,sep="_")
alltaxonomy$subsepecieslca<-gsub(";","_",alltaxonomy$V2)

#get each path out as separate files-
phylalca<-alltaxonomy%>%dplyr::select(phyla,phylalca)%>%unique()
colnames(phylalca)<-c("name","lca")
write.csv(phylalca,file="gtdb_ver95_phylalca.csv")

classlca<-alltaxonomy%>%dplyr::select(class,classlca)%>%unique()
colnames(classlca)<-c("name","lca")
write.csv(classlca,file="gtdb_ver95_classlca.csv")

orderlca<-alltaxonomy%>%dplyr::select(order,orderlca)%>%unique()
colnames(orderlca)<-c("name","lca")
write.csv(orderlca,file="gtdb_ver95_orderlca.csv")

familylca<-alltaxonomy%>%dplyr::select(family,familylca)%>%unique()
colnames(familylca)<-c("name","lca")
write.csv(familylca,file="gtdb_ver95_familylca.csv")

genuslca<-alltaxonomy%>%dplyr::select(genus,genuslca)%>%unique()
colnames(genuslca)<-c("name","lca")
write.csv(genuslca,file="gtdb_ver95_genuslca.csv")

specieslca<-alltaxonomy%>%dplyr::select(species,specieslca)%>%unique()
colnames(specieslca)<-c("name","lca")
write.csv(specieslca,file="gtdb_ver95_specieslca.csv")

subspecieslca<-alltaxonomy%>%dplyr::select(V1,subsepecieslca)%>%unique()
colnames(subspecieslca)<-c("name","lca")
write.csv(subspecieslca,file="gtdb_ver95_subspecieslca.csv")

#get the final lca file-
alllca<-rbind(phylalca,classlca,orderlca,familylca,genuslca,specieslca,subspecieslca)%>%as.data.frame()
write.csv(alllca,file="gtdb_ver95_alllca.csv")

#merge to taxid file-
taxid<-read.csv("/apps/unit/BioinfoUgrp/DB/diamondDB/GTDB/r95/taxID_info.tsv",header=TRUE,sep="\t")
taxid_taxa<-merge(taxid,taxa,by="name",all.x=TRUE)
write.csv(taxid_taxa,file="gtdb_ver95_alllca_taxid.csv")

##NOTE- nrow(alllca) == nrow(taxID_info.tsv) i.e. no. of elements completely match.

#---------------------------------------------------------------------------
##get the taxonomy per marker gene from DIAMOND blast output- 

#to run-
module load R/3.6.1
Rscript gtdb_diamondlca.R gtdb_ver95_alllca_taxid.csv output-gtdb-lca-method1-matches-prot-all-COG0087.fasta.txt taxa-output-gtdb-lca-method1-matches-prot-all-COG0087.fasta.txt

#---------------------------------------------------------------------------
<called as "gtdb_diamondlca.R">
#USAGE- Rscript [taxaid] [blast] [output_file]

args <- commandArgs(TRUE)
taxaid <- args[1] #gtdb_ver95_alllca_taxid.csv file
blast <- args[2]
output_file <- args[3]

taxa<-read.csv(file=taxaid,header=TRUE)
taxa$X<-NULL

blast<-read.csv(file=blast,header=FALSE,sep="\t")
library(dplyr)
blast2<-blast%>%filter(V3<=1e-24)

#merge-
taxa_blast<-merge(blast2,taxa,by.x="V2",by.y="taxID")
write.csv(taxa_blast,file=output_file)



```


```
##comparison with other taxonomic methods-

###CONTIGS-
#NOTE- for the 3 contigs of sample 230-12, the taxonomic annotation of diamond lca (102 outfmt) is at a higher taxonomic level than MEGAN lca.
#eg-
#MEGAN output-
230-12_NODE_61number61  GTDB;Bacteria;Firmicutes_A;Clostridia;Lachnospirales;Lachnospiraceae;
230-12_NODE_59number59  GTDB;Bacteria;Spirochaetota;Spirochaetia;Treponematales;Treponemataceae_B;Treponema_E;
230-12_NODE_189number189 GTDB;Bacteria;Fibrobacterota;Chitinivibrionia;Chitinivibrionales;Chitinispirillaceae;Chitinispirillum;Chitinispirillum alkaliphilum;Chitinispirillum alkaliphilum ACht6-1;

#DIAMOND lca output-
"230-12_NODE_61number61","d__Bacteria_p__Firmicutes_A_c__Clostridia"
"230-12_NODE_59number59","d__Bacteria_p__Spirochaetota_c__Spirochaetia_o__Treponematales"
"230-12_NODE_189number189","d__Bacteria_p__Fibrobacterota_c__Chitinivibrionia_o__Chitinivibrionales_f__Chitinispirillaceae"

#BLASTp output and Blobtools-
"230-12_NODE_61number61" has 3/6 proteins annotated at family level in BLASTp. So, class level LCA in diamond lca might be very unspecific(?!) BLASTp, MEGAN and Blobtools disagree at genus/family level for this contig. So maybe order level classification should have been used by diamond lca but instead class level is used.
"230-12_NODE_59number59" has 5/5 proteins annotated at genus level in BLASTp. Blobtools annotates this contig at species level. MEGAN annotates the contig at genus level. So, order level LCA in diamond lca is not specific.
"230-12_NODE_189number189" has 8/9 proteins annotated at genus level in BLASTp. MEGAN and Blobtools annotation are at species level. But BLASTp annotation does not match the other two methods (g_GUT77 vs Chitinispirillum alkaliphilum ACht6-1). So, family level LCA classification is good.

##NOTE- DIAMOND lca is more stringent in its taxonomic classification thereby giving outputs at higher taxonomic levels (class, order and family). But the four methods cannot be dirrectly compared as MEGAN and Blobtools first blast the contigs to nr database and then map to GTDB database, while BLASTp is against proteins rather than contigs so it is using a different blast algorithm. 

#four different diamond blast parameter combinations examined for contig data-
--sensitive --query-cover 60 --id 65 --evalue 1e-24 : no result
--query-cover 30 :no result
--max-hsps 0 --evalue 1e-24 : gives results
--max-hsps 0 --query-cover 30 : no results

#"--max-hsps 0 --evalue 1e-24" parameter was chosen.

```
