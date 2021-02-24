## CAZY annotation

### a)map the protein sequences to CAZY database
```
#!/bin/bash
#SBATCH --job-name=hummrCAZY250above
#SBATCH --partition=compute
#SBATCH --time=7-0
#SBATCH --mem=80G
#SBATCH --ntasks=1
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=hmmr250_nt_%A-%a.out
#SBATCH -o out_hmmr250.%j
#SBATCH -e err_hmmr250.%j

##SBATCH --array 1-100
##num=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

#load modules-
module load hmmer/3.1b2

cd /work/BourguignonU/jigyasa/microcerotermes/functional_annotation/all_prokka_outputs

#files=(*.faa)
#file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/-prokka.faa/}
#file3=${file2/.txt/}


DB_DIR="/work/student/jigyasa-arora/BourguignonU_data/main_ID230+ID231_2/ID230/joinedfiles/K21-kate_assembly/prodigal/500above/cazy_output/hummer_cazy"


#for metagenomes use e-5 (dbCAN db recommendation, Ye et al, 2012)
#non-stringent hmmscan-
hmmscan -E 1.0e-5 -o 229-01-HMMoutputfile --tblout 229-01-HMM-persequence-output --domtblout 229-01-HMM-perdomain-output ${DB_DIR}/dbCAN-fam-HMMs.txt 229-01-prokka.faa

#NOTE- Even though a non-stringent evlaue cutoff used in the initial analysis, the final version has evalue <1e-30 (see below).
```

### b)get the hmmscan output file to a text delimited file

`sh hmmscan-parser.sh 229-01-HMM-perdomain-output > dbcan-229-01-HMM-perdomain-output`

## c) use a more stringent evalue cutoff for final CAZYme annotation

`awk -F"," '$7<1e-30 {print $2","$3","$4}' dbcan-229-01-HMM-perdomain-output > evalue30-229-01-HMM-perdomain-output`



## KEGG annotation

## Pfam annotation
