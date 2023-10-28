#!/bin/sh
#SBATCH --job-name n10x_count
#SBATCH --mem 50G
#SBATCH --tasks 50
#SBATCH --output logs/%x.%j.stdout
#SBATCH --error logs/%x.%j.stderr
#SBATCH --partition cpu-med
#SBATCH --time 08:00:00 ## HH:MM:SS


# change directory to your desired output folder
cd /research/labs/neurology/fryer/m214960/PLX5622_mice/counts

# source settings
source $HOME/.bash_profile

# get cellranger version
cellranger -V

# print sample
sample=$1
echo "sample: $sample"

# run cellranger
cellranger count \
	--id=$sample \
	--sample=$sample \
	--fastqs=/research/labs/neurology/fryer/projects/PLX5622_diet \
	--transcriptome=/research/labs/neurology/fryer/projects/references/mouse/refdata-gex-mm10-2020-A \
	--localcores=16 \
	--localmem=50 


