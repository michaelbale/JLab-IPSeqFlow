#! /bin/bash -l
 
#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=callFlow # change as desired
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=mib4004@med.cornell.edu # change
#SBATCH --time=120:00:00   # HH/MM/SS; 5 days should be sufficient. If you have a lot of data, consider breaking it down into subfolders
#SBATCH --mem=16G   # memory requested, units available: K,M,G,T
 
 
 
# Follow your own needs for replacing {OPTIONS} according to your data
# source bashrc may or may not be necessary depending on how you call the script
source ~/.bashrc

# if nextflow is installed via conda, the spack load is not necessary
# if conda needs to be activated (i.e. you don't see a (base) next to your ID) use condaon
#condaon
#spack load nextflow

nextflow run michaelbale/JLab-Flow {OPTIONS}
 
exit
