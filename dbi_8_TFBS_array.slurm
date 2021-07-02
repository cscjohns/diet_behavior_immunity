#!/bin/bash

## Job Name
#SBATCH --job-name=TFBS

## Allocation Definition
#SBATCH --account=csde-ckpt
#SBATCH --partition=ckpt

## Nodes
#SBATCH --nodes=1
#SBATCH --spread-job

## Wall Time
#SBATCH --time=4:00:00

## Memory per Node
#SBATCH --mem=120G

## Stdout and Stderr
#SBATCH --output=./slurm_%j.out
#SBATCH --error=./slurm_%j.err

## Email Notification
#SBATCH --mail-user=cscjohns@uw.edu
#SBATCH --mail-type=ALL

## sbatch -p ckpt -A csde-ckpt --array=1-8 ../../scripts/TFBS_array.slurm

## load module for homer

module load contrib/homer/4.10 

## Define varibles for focal and background gene sets as well as results directory names

focal=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./focal_gene_sets.txt)
background=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./bg_gene_sets.txt)
dir_name=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./directories.txt)
upstream=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./upstream.txt)

## run findMotifs for each set of focal and background genes

findMotifs.pl $focal Mfas_5.0_100kb $dir_name -start $upstream -end 2000 -bg $background -mset vertebrates

exit;
