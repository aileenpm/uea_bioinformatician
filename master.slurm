#!/usr/bin/env bash

#SBATCH --partition jic-long,jic-medium           # using SSD and only jic* nodes have them (spec. j512* nodes)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time 1-0:00                             # time (D-H:m)
#SBATCH --output ./slurm_output/align.%A_%a.out   # STDOUT
#SBATCH --mem=40G                                 # memory pool for ALL cores
#SBATCH --array=1-10                              # default - override in 'sbatch' command
#SBATCH --gres=ssd
#SBATCH --localscratch=ssd:60                     # request 60GB SSD - enough?

# set up functions from Hugh
function info() {
  echo "--> $(date "+%F %T") INFO : $@" >&2
}

cpus=${SLURM_CPUS_PER_TASK}
ssd=${SLURM_LOCAL_SCRATCH}

# enter bash "strict" mode
set -euo pipefail

# activate conda env
info "activating conda environment..."

conda init bash
. "/hpc-home/magilin/miniconda3/etc/profile.d/conda.sh"
conda activate /jic/scratch/projects/bravo/magilin/resources/mirna-alignment/conda-envs/mirna-identification-env

# list directory names
info "creating directories..."
d=$(pwd -P)

qc_dir="${d}/qc"

qc_raw_dir="${qc_dir}/raw"
qc_trim_dir="${qc_dir}/trim"
bowtielog_dir="${d}/bowtielog"

intermediate_dir="${d}/intermediate_data"
filtered_fa="${intermediate_dir}/filtered_fa"

final_dir="${d}/final_data"
mirna_bam_dir="${final_dir}/bam"
mirna_counts_dir="${final_dir}/counts"


mkdir -p "${qc_dir}"
mkdir -p "${qc_raw_dir}"
mkdir -p "${qc_trim_dir}"
mkdir -p "${bowtielog_dir}"
mkdir -p "${intermediate_dir}"
mkdir -p "${filtered_fa}"
mkdir -p "${final_dir}"
mkdir -p "${mirna_bam_dir}"
mkdir -p "${mirna_counts_dir}"

# get the fastq file/dir for this array job
fq_files=$1                     # the file containing a list of fastq files/dirs
raw_fq=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" ${fq_files})

info "processing $fq_files (line ${SLURM_ARRAY_TASK_ID}) : $raw_fq"

# move to ssd
info "entering ssd directory ${ssd} ..."
cd $ssd

# copy bowtie indexes (quicker than building from scratch)
info "copying bowtie indexes ..."
srun cp /hpc-home/magilin/magilin-s/resources/mirna-alignment/bowtie-indexes/*.ebwt .


# 1. qc raw data
info "running fastqc - initial raw read QC..."

srun fastqc \
  ${raw_fq} \
  -t ${cpus} \
  -o "${qc_raw_dir}"

# 2. trim adapters
info "trimming adapters..."

sample_id=$(basename ${raw_fq} | sed 's/.fastq.gz//')
trim_fq="${sample_id}_trim.fq"         
  
srun cutadapt \
  -j ${cpus} \
  -q 15 \
  -u 1 \
  -a TGGAATTCTCGGGTGCCAAGG \
  -m 15 \
  -M 30 \
  -o $trim_fq \
  $raw_fq

# 3. qc trimmed data
info "running fastqc - trimmed reads QC..."

srun fastqc \
  ${trim_fq} \
  -t ${cpus} \
  -o "${qc_trim_dir}"
  
# 4. convert fq to fa for bowtie
info "convert .fq to .fa for bowtie..."

trim_fa="${sample_id}_trim.fa"

cat $trim_fq \
| paste - - - - \
| sed 's/^@/>/g' \
| cut -f 1,2 \
| tr '\t' '\n' > $trim_fa

# 6. Align reads to genome
info "step 01: aligning clean reads to genome..."

s01_genome_aligned_fa="${sample_id}_s01_genome_aligned.fa"
s01_genome_bowtielog_file="${bowtielog_dir}/${sample_id}_s01_genome_aligned_bowtielog.txt"

#s01_genome_aligned_sam="${sample_id}_s01_genome_aligned.sam"                             #uncomment for .sam file
#s01_genome_aligned_bam_file="${mirna_bam_dir}/${sample_id}_s01_genome_aligned.bam"       #uncomment for .bam file

srun bowtie ro18-genome \
  --threads ${cpus} \
  -n 0 \
  -f $trim_fa \
  --quiet \
  --al $s01_genome_aligned_fa \
  2> $s01_genome_bowtielog_file
  #  -S $s01_genome_aligned_sam \          #uncomment to get .sam file

#srun samtools view -bS $s01_genome_aligned_sam | samtools sort -o $s01_genome_aligned_bam_file -O bam -      #uncomment to get .bam file


# 7. Filter genome-aligned reads against ncrna
info "step 02: filter out ncrna-aligned reads..."

s02_ncrna_filtered_fa="${sample_id}_s02_ncrna_filtered.fa"
s02_ncrna_bowtielog_file="${bowtielog_dir}/${sample_id}_s02_ncrna_filtered_bowtielog.txt"

srun bowtie ncrna \
  --threads ${cpus} \
  -f $s01_genome_aligned_fa \
  --quiet \
  --un $s02_ncrna_filtered_fa \
  2> $s02_ncrna_bowtielog_file

# 8. Filter out repeats
info "step 03: filter out repeat-aligned reads..."

s03_repeats_filtered_fa="${sample_id}_s03_repeats_filtered.fa"
s03_repeats_bowtielog_file="${bowtielog_dir}/${sample_id}_s03_repeats_filtered_bowtielog.txt"

srun bowtie repeats \
  --threads ${cpus} \
  -f $s02_ncrna_filtered_fa \
  --quiet \
  --un $s03_repeats_filtered_fa \
  2> $s03_repeats_bowtielog_file

# 9. Align mirnas
info "step 04: aligning filtered reads to mirna sequences..."

s04_mirna_filtered_fa="${sample_id}_s04_mirna_filtered.fa"
s04_mirna_aligned_sam="${sample_id}_s04_mirna_aligned.sam"

s04_mirna_bowtielog_file="${bowtielog_dir}/${sample_id}_s04_mirna_filtered_bowtielog.txt"

srun bowtie bra-mirna \
  --threads ${cpus} \
  --best \
  -v 0 \
  -f $s03_repeats_filtered_fa \
  --quiet \
  -S $s04_mirna_aligned_sam \
  --un $s04_mirna_filtered_fa \
  2> $s04_mirna_bowtielog_file

# 10. Align to CDS
info "step 05: aligning non-mirna reads to coding sequences..."

s05_cds_filtered_fa="${sample_id}_s05_cds_filtered.fa"
s05_cds_bowtielog_file="${bowtielog_dir}/${sample_id}_s05_cds_filtered_bowtielog.txt"

srun bowtie bra-cdna \
  --threads ${cpus} \
  -f $s04_mirna_filtered_fa \
  --quiet \
  --un $s05_cds_filtered_fa \
  2> $s05_cds_bowtielog_file

# 11. Align to genome again
info "step 06: align unidentified sequences to genome..."

s06_unidentified_aligned_fa="${sample_id}_s06_unidentified_aligned.fa"
s06_unidentified_aligned_sam="${sample_id}_s06_unidentified_aligned.sam"
s06_unidentified_bowtielog_file="${bowtielog_dir}/${sample_id}_s06_unidentified_aligned_bowtielog.txt"

srun bowtie ro18-genome \
  --threads ${cpus} \
  -f $s05_cds_filtered_fa \
  --quiet \
  --al $s06_unidentified_aligned_fa \
  -S $s06_unidentified_aligned_sam \
  2> $s06_unidentified_bowtielog_file

# make .bam files and copy to final_data
info "convert unidentified reads to .bam files..."

s06_unidentified_aligned_bam="${sample_id}_s06_unidentified_aligned.bam"
srun samtools view -bS $s06_unidentified_aligned_sam | samtools sort -o $s06_unidentified_aligned_bam -O bam -

srun cp $s06_unidentified_aligned_bam $mirna_bam_dir/

# 12. Get mirna readcounts
info "get mirna readcounts..."

s04_mirna_aligned_bam="${sample_id}_s04_mirna_aligned.bam"
s04_mirna_counts="${sample_id}_s04_mirna_counts.txt"

samtools view -bS $s04_mirna_aligned_sam | samtools sort -o $s04_mirna_aligned_bam -O bam - 
    
samtools index $s04_mirna_aligned_bam
    
samtools idxstats $s04_mirna_aligned_bam | cut -f 1,3 > $s04_mirna_counts

srun cp $s04_mirna_counts $mirna_counts_dir/
srun cp $s03_repeats_filtered_fa $filtered_fa/

conda deactivate

info "completed alignments and mirna read counts!"
