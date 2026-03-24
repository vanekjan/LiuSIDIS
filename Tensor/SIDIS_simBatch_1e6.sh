#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=SIDIS_sim
#SBATCH --time=0-05:00:00
#SBATCH --mem=4G
#SBATCH --qos=normal
#SBATCH --output=/w/hallc-scshelf2102/tensor/vanekjan/LiuSIDIS/Tensor/output/jobs/%x-%j-%Ni_%a.out
#SBATCH --error=/w/hallc-scshelf2102/tensor/vanekjan/LiuSIDIS/Tensor/output/jobs/%x-%j-%N_%a.err

./analysis_deuteron 1000000 0 
 
