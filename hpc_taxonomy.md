## nr database to get prokaryotic and non-prokaryotic contigs out

## GTDB database to annotate the prokaryotic contigs

```
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 3-0
#SBATCH --mem=150G
#SBATCH -e err_wood.%j

#SBATCH --array 0-226
#load the libraries-
module load ncbi-blast/2.10.0+
module load edirect/13.8
eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib` #to run!
echo 'eval `perl -I ~/perl5/lib/perl5 -Mlocal::lib`' >> ~/.profile #to run!
echo 'export MANPATH=$HOME/perl5/man:$MANPATH' >> ~/.profile #to run!

#call the folders-
IN_DIR="/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/222samples_taxonomy/gtdb/microbial_contigsequences"
DB_DIR="/flash/BourguignonU/DB/GTDB_Kraken/gtdb_r89_54k"
OUT_DIR="/flash/BourguignonU/Jigs/gtdb/kraken_all"

cd ${IN_DIR}/
files=(*fasta)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/filename-/}
#file3=${file2/.fasta/}

#seqtk subseq ${IN_DIR1}/${file1} ${IN_DIR2}/names-bacterial-2-newname-megan-matches-${file3}.txt > ${OUT_DIR}/bacterial-contigs-${file3}.fasta #extract prokaryotic contigs from nr
# Start 'myprog' with input from bucket, and output in /flash/

kraken2 --db ${DB_DIR} ${IN_DIR}/${file1} --output ${OUT_DIR}/kraken-output-${file1}.txt --report ${OUT_DIR}/kraken-report-${file1}.txt --report-zero-counts --use-mpa-style --use-names --threads 10

```

## generating the final taxonomy file, with annotation per contig-

`

```
#save as "krakenoutputs.R"



```

## check1- Does each contig has one annotation?

```
module load R/3.6.1
# in R-
krakenoutput<-read.csv("/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/222samples_taxonomy/gtdb/krakenoutput/2-kraken-output-allprokaryotic_contigs-229-01.fasta.txt",header=FALSE,sep="\t")

library(dplyr)

contigs<-krakenoutput%>%select(V1)%>%unique()
nrow(contigs)
[1] 716080
nrow(krakenoutput)
[1] 716080
```

## check2-how many prokaryotic contigs (from nr) were annotated by GTDB-
