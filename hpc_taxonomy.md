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

### C) Based on https://github.com/Jigyasa3/termite_146guts_microbes_function/blob/main/standardization/taxonomy_standardization_krakenvsmeganvsblastp.md, use MEGAN for taxonomic analysis-

```
#a)
diamond blastp --db ${DB_DIR}/nr.dmnd --query ${IN_DIR}/${file1} --outfmt 100 --out ${OUT_DIR}/nr-matches-${file1} --threads 15
#--outfmt 100- DAA file output

#b)
/home/j/jigyasa-arora/local/megan/tools/blast2rma --in ${file1} --format DAA --blastMode BlastP --out ${file1}.rma --minScore 140 --minReadLength 50 --maxExpected 1e-25 --minPercentIdentity 65 --lcaAlgorithm naive  --mapDB ${DB_DIR}/megan-map-Jan2021.db --threads 16 --verbose
#--minScore 100
#filters mainly used from http://megan.informatik.uni-tuebingen.de/t/missing-family-rank-in-blast2lca-output/577/5
#--minReadLength 50 == trimmomatic results

#c)
/home/j/jigyasa-arora/local/megan/tools/rma2info --in ${file1} --read2class GTDB --paths --out meganoutput-${file1}.txt
#based on http://megan.informatik.uni-tuebingen.de/t/missing-family-rank-in-blast2lca-output/577/5
#--read2class gives GTDB taxonomy to each protein.


```
