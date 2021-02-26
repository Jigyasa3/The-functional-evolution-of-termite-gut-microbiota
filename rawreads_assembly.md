## get read quality information

```
#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=largemem
#SBATCH --time=7:00:00
#SBATCH --mem=124G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END


#SBATCH --output=slurm-%j.base
#SBATCH -o out_nfastqc.%A_%a
#SBATCH -e err_nfastqc.%A_%a

# command to run
module load fastqc/0.11.5

for i in /work/student/jigyasa-arora/BourguignonU_data/main_ID230+ID231_2/ID230/230run_nomismatch/230-*.fastq.gz
do
        fastqc -o /work/student/jigyasa-arora/BourguignonU_data/main_ID230+ID231_2/ID230/230run_nomismatch/fastqc/ $i
done


for i in /work/student/jigyasa-arora/BourguignonU_data/main_ID230+ID231_2/ID230/230run_nomismatch/trimmed_output/230-*.fq.gz
do
        fastqc -o /work/student/jigyasa-arora/BourguignonU_data/main_ID230+ID231_2/ID230/230run_nomismatch/trimmed_output/fastqc/ $i
done

```

## trim raw reads

```
#!/bin/bash
#SBATCH -p largemem
#SBATCH -t 10-0
#SBATCH --mem=80G
#SBATCH -e err_mt_mafft_protein.%j

#SBATCH --array 0-50
##num=$(printf "%04d" $SLURM_ARRAY_TASK_ID) #add two zeros infront
#call the folders-
IN_DIR="/bucket/BourguignonU/Jigs_backup/301run/301_nov1"
PAIRED_DIR="/flash/BourguignonU/Jigs/tag_switching301/paired"
UNPAIRED_DIR="/flash/BourguignonU/Jigs/tag_switching301/unpaired"
#FINAL_DIR="/bucket/BourguignonU/Jigs_backup/working_files/tagswitching/trimmomatic/301_may"

cd ${IN_DIR}/
files=(*_R1_001.fastq.gz)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/R1/R2}
file3=${file1/_R1_001.fastq.gz/}

#load modules from amd modules
module load amd-modules
module load Trimmomatic/0.33


# Start 'myprog' with input from bucket, and output in /flash/

java -jar /apps/free72/Trimmomatic/0.33/lib/trimmomatic-0.33.jar PE ${IN_DIR}/${file1} ${IN_DIR}/${file2} ${PAIRED_DIR}/${file3}_R1_paired.fq.gz ${UNPAIRED_DIR}/${file3}_R1_unpaired.fq.gz ${PAIRED_DIR}/${file3}_R2_paired.fq.gz ${UNPAIRED_DIR}/${file3}_R2_unpaired.fq.gz SLIDINGWINDOW:4:30 MINLEN:50 HEADCROP:16 LEADING:3 TRAILING:3

##explaination-
#LEADING-3 -> to remove leading bad reads below quality 3
#TRAILING-3 -> to remove trailing bad reads below quality 3
#SLIDINGWINDOW:4:30 -> 30 is the quality of base. 30 is good PHRED quality score
#MINLEN=50 ->as default 36 is too low and when 100 was used, it retained reads below 100 bps
#HEADCROP=16 ->removing the contaminating barcode reads from the 5' end
```

## contig assembly

```
#!/bin/bash
#SBATCH -p largemem
#SBATCH -t 14-0
#SBATCH --mem=180G
#SBATCH -e err_spades.%j

#SBATCH --array 0-53
##num=$(printf "%04d" $SLURM_ARRAY_TASK_ID) #add two zeros infront
#call the folders-
IN_DIR="/flash/BourguignonU/Jigs/spades/272_nov1"

cd ${IN_DIR}/
files=(*R1_paired.fq.gz)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
file2=${file1/R1/R2}
file3=${file1/_R1_paired.fq.gz/}

#load modules from amd modules
module load amd-modules
module load python/2.7.18
module load SPAdes/3.11.1

# Start 'myprog' with input from bucket, and output in /flash/
spades.py --meta --only-assembler -k 21,31,41,51,71 --pe1-1 ${IN_DIR}/${file1} --pe1-2 ${IN_DIR}/${file2} -o ${file3}_smallerkmer.fasta


```
