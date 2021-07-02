#!/bin/bash

## Job Name
#SBATCH --job-name=pval_distribution

## Allocation Definition
#SBATCH --account=csde-ckpt
#SBATCH --partition=ckpt

## Wall Time
#SBATCH --time=10:00:00

## Stdout and Stderr
##SBATCH --output=./pval_dist.out
##SBATCH --error=./pval_dist.err

## Email Notification
#SBATCH --mail-user=cscjohns@uw.edu
#SBATCH --mail-type=ALL

## sbatch --array=1-100 ../scripts/dbi_4_pval_DAB_distribution.slurm

## define variables
taskID=${SLURM_ARRAY_TASK_ID}

sbatch -p ckpt -A csde-ckpt --time=4:00:00 --mem=120G ../scripts/dbi_4_perm_DAB_pvals.R $taskID
