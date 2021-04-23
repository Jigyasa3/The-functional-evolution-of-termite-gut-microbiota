## taxonomic annotation of contigs

## taxonomic annotation of marker-genes

### A) extract marker-genes
```
#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --partition=compute
#SBATCH --time=1-0
#SBATCH --mem=180G
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=iqtree_nt_%A-%a.out
#SBATCH -o out_iqtree.%j
#SBATCH -e err_iqtree.%j
#SBATCH --array 0-225

#load modules-
DNA_DIR="/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/all_prokka_outputs/all_fna_sequences_Dec2019"
cd ${DNA_DIR}/

files=(2-renamed-*fasta)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/2-renamed-/}
file3=${file2/-prokka.fasta/}

PROTEIN_DIR="/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/all_prokka_outputs/all_protein_sequences_Dec2019/with_filenames"
OUT_DIR="/flash/BourguignonU/Jigs/markers/allmarkers_re"

fetchMGs.pl -m extraction -x /home/j/jigyasa-arora/local/fetchMGs/bin -d ${DNA_DIR}/${file1} ${PROTEIN_DIR}/filename-${file3}-prokka.faa-smallername.faa -o ${OUT_DIR}/${file3}-fetchmoutput -t 16 ##bin directory where you cloned the repo.
```

### B) Group data by each COG-i.e. 40COG files containing all samples' marker genes
```
while read line;do while read cogs;do cp ${line}/${cogs}*faa allfetchm_nucoutput/${line}-${cogs}.faa;done < allcogs.txt ;done <filesnames.txt

while read line;do while read cogs;do cp ${line}/${cogs}*fna allfetchm_nucoutput/${line}-${cogs}.fna;done < allcogs.txt ;done <filesnames.txt

```

### C) Based on https://github.com/Jigyasa3/termite_146guts_microbes_function/blob/main/standardization/taxonomy_standardization_krakenvsmeganvsblastp.md, use DIAMOND lca for taxonomic analysis-

```
#!/bin/bash
#SBATCH -p compute
#SBATCH -t 3-0
#SBATCH --mem=350G
#SBATCH -e err_wood.%j

#SBATCH --array 0
#load the libraries-
module load bioinfo-ugrp-modules
module load DB/diamondDB/GTDB/r95
module load DIAMOND/2.0.4.142

#call the folders-
IN_DIR="/bucket/BourguignonU/Jigs_backup/working_files/AIMS/AIM2/tpm_functional_annotation/functional_annotation/all_prokka_outputs/all_contigs_Dec2019"
OUT_DIR="/flash/BourguignonU/Jigs/gtdb/gtdb_allcontigs"

cd ${IN_DIR}/
files=(filename-272-15.fasta)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/filename-/}
#file3=${file2/.fasta/}

diamond blastx --db ${DIAMONDDB}/gtdb.dmnd --query ${IN_DIR}/${file1} --outfmt 102 --out ${OUT_DIR}/gtdb-lca-method2-matches-${file1}.txt --max-hsps 0 --evalue 1e-24 --threads 15
#--max-hsps 0 : to get all best hits per query
#--evalue 1e-24
#--query-cover 60 was tested, but it is too stringent for termite microbiome that doesnt have a match in the database. the results were empty.

#-------------------------------------------------------------------------------------------------------------------------------------
#get the annotations that are not zero from the above script-

awk '$2!=0 {print $0}' gtdb-lca-method2-matches-prot-all-COG0087.fasta.txt > output-gtdb-lca-method2-matches-prot-all-COG0087.fasta.txt

#-----------------------------------------------------------------------------------------------------------------------------------------
##get the taxonomy per marker gene from DIAMOND blast output- 
module load R/3.6.1
Rscript gtdb_diamondlca.R gtdb_ver95_alllca_taxid.csv output-gtdb-lca-method1-matches-prot-all-COG0087.fasta.txt taxa-output-gtdb-lca-method1-matches-prot-all-COG0087.fasta.txt

#------------------------------------------------------------------------------------------------------------------------------------------
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
