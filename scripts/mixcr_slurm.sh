#!/bin/bash
#SBATCH --job-name=mixcr_ht        # Job name
#SBATCH --cpus-per-task=110         # Run on a single CPU
#SBATCH --mem=100gb                 # Job memory request
#SBATCH --time=10:00:00           # Time limit hrs:min:sec
#SBATCH --output=JobName.%j.log   # Standard output and error log
#SBATCH --mail-user=your_amazing_email@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=medium


# script to extract repertoires from healthy tissues
source ~/.bashrc
conda activate mixcr

echo "Start analysis" > run_log.txt

for file in *_1.fastq.gz; do
    SAMPLE=${file%_1.fastq.gz}
    
    echo "$(date): Processing $SAMPLE" >> run_log.txt

    mixcr analyze rna-seq --species hsa --threads 110 \
        --no-json-reports \
        "${SAMPLE}_1.fastq.gz" "${SAMPLE}_2.fastq.gz" "../clontypes/$SAMPLE"
    
        
    echo "$(date): $SAMPLE analysis done" >> run_log.txt
done

echo "ALL DONE" >> run_log.txt

conda deactivate
