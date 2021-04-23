## Psuedo genes are rampant in endosymbionts, and they have been examined in endosymbionts of flagellates in lower termites (Hongoh-san lab). 
## Ways to detect Pseudo genes in protein sequences?

### a)

### b) Using no. of asterisk inside the translated protein sequences as a proxy of pseudo genes (genetic code 11)
```
#Q/ How many bacterial proteins have pseudo genes? <using gtdb blastx for taxonomic annotation>

#get microbial protein sequences out- (using the "proteinnames-301-20.txt" for example)
for i in filename*.faa;do filename=`echo ${i}| sed 's/filename-//g'`; filename2=`echo ${filename}| sed 's/-prokka.faa-smallername.faa//g'` ; seqtk subseq ${i} ${IN_DIR}/proteinnames-${filename2}.txt > ../gtdbmicrobialproteins_march2021/microbial-${filename2}.faa ;done

#get proteins with asterisk in them-
for i in microbial-*.faa ;do awk '/^>/ {printf("\n%s\t",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${i} | grep "*" > asterisk-${i};done

#count the no. of proteins in original file and asterisk containing file-
for i in *.faa;do grep -c ">" ${i} > count-${i}.txt;done

#add filename to the count file-
for i in count-*.txt; do awk '{print FILENAME"\t"$0}' ${i} > 2-${i};done

#in R-
> asterisk<-read.csv("all-asterisk-proteins.txt",header=FALSE,sep="\t")
> asterisk$V1<-gsub("count-asterisk-microbial-","",asterisk$V1)
> asterisk$V1<-gsub(".faa.txt","",asterisk$V1)
> colnames(asterisk)<-c("samples","asterisk")

> microbial<-read.csv("all-proteins.txt",header=FALSE,sep="\t")
> microbial$V1<-gsub("count-microbial-","",microbial$V1)
> microbial$V1<-gsub(".faa.txt","",microbial$V1)
> colnames(microbial)<-c("samples","all")
> joined<-merge(asterisk,microbial,by="samples")

> joined$perc<-(joined$asterisk/joined$all)*100
> mean(joined$perc)
[1] 1.608938
> max(joined$perc)
[1] 8.818765
> min(joined$perc)
[1] 0

##NOTE- Overall there are 1.6% microbial proteins on average with internal stop codon.
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Q/ How many bacterial proteins functionally annotated for the paper's stats have pseudo genes?

#extract proteinids from "allannotations_tpm_contiginfo_taxonomy2.csv" (this file is used for all downstream stats)
awk -F"," '{print $4}' allannotations_tpm_contiginfo_taxonomy2.csv > allannotated-proteins.txt

#get annotated microbial proteins out-
cat microbial-*.faa > allmicrobialproteins.faa
seqtk subseq allmicrobialproteins.faa allannotated-proteins.txt > allmicrobialannotatedproteins.faa

#get proteins with asterisk-
awk '/^>/ {printf("\n%s\t",$0);next; } { printf("%s",$0);}  END {printf("\n");}' allmicrobialannotatedproteins.faa | grep "*" > asterisk-allmicrobialannotatedproteins.faa

grep -c ">" asterisk-allmicrobialannotatedproteins.faa
[1] 10302
grep -c ">" allmicrobialannotatedproteins.faa
[1] 4543277

##NOTE- Using genetic code 11, there are 0.2% proteins in the annotated microbial sequences that have an internal stop codon.
```
