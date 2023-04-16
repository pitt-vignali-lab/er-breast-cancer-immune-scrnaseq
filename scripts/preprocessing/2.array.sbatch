#!/bin/bash
#SBATCH -p htc
#SBATCH --nodes 1
#SBATCH --job-name count+velco
#SBATCH -t 2-00:00:00
#SBATCH --cpus-per-task=10
#SBATCH --output=Logs/%A_%a.out
#SBATCH --mem=150g

module load cellranger/6.1.1
module load velocyto/0.17
module load gcc/10.2.0
module load samtools/1.9

cd ./Outcomes

FASTQS=($((cd ../ && ls -d *GEX) | head -$SLURM_ARRAY_TASK_ID | tail -1))

mv ../${FASTQS}/*R1*gz ../${FASTQS}/${FASTQS}_S1_L001_R1_001.fastq.gz
mv ../${FASTQS}/*R2*gz ../${FASTQS}/${FASTQS}_S1_L001_R2_001.fastq.gz

cellranger count --id ${FASTQS} \
                 --fastqs=../${FASTQS} \
                 --transcriptome=/bgfs/genomics/refs/CellRanger/refdata-gex-GRCh38-2020-A \
                 --localcores=10 \
                 --localmem=150

#### velocyto #####

samtools sort -@ 10  -t CB -O BAM -o \
./${FASTQS}/outs/cellsorted_possorted_genome_bam.bam \
./${FASTQS}/outs/possorted_genome_bam.bam

rmsk_gtf=/bgfs/dvignali/cuij/hg38RepeatRmsk.gtf 

cellranger_gtf=/bgfs/genomics/refs/CellRanger/refdata-gex-GRCh38-2020-A/genes/genes.gtf

cellranger_outDir=./${FASTQS}

velocyto run10x -m $rmsk_gtf  $cellranger_outDir $cellranger_gtf

mv ./${FASTQS}/velocyto/${FASTQS}.loom ./