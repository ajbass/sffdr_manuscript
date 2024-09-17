#!/bin/bash
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J ukbb
#! Which project should be charged (NB Wilkes2 projects end in '-GPU'):
#SBATCH -A CWALLACE-INTREPID-SL2-CPU
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total?
#! Note probably this should not exceed the total number of GPUs in use.
#SBATCH --ntasks=4
#! How much wallclock time will be required?
#SBATCH --time=04:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=NONE
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#SBATCH --mem=60000mb
#SBATCH --array=1-22%2

#! Do not change:
#SBATCH -p cclake

#! sbatch directives end here (put any additional directives above this line)

#! Notes:
#! Charging is determined by GPU number*walltime.

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
#! module purge                               # Removes all modules still loaded
#! module load ab3105/default-ccl              # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:
module load plink/2.00-alpha

#! Full path to application executable:
application=""

#! Run options for the application:
options=""

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"  # The value of SLURM_SUBMIT_DIR sets workdir to the directory
                             # in which sbatch is run.

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 56:
export OMP_NUM_THREADS=1

#! Number of MPI tasks to be started by the application per node and in total (do not change):
#!np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
#!export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
#!export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.


#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
#!CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"


###############################################################
### You should not have to change anything below this line ####
###############################################################

#!cd $workdir
#!echo -e "Changed directory to `pwd`.\n"

#!JOBID=$SLURM_JOB_ID

#!echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
chr=${SLURM_ARRAY_TASK_ID}

for id in {1..4}; do \
plink2 --bfile /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2 --fam /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2_s488131.fam --bim /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr${chr}_v2.bim --keep /home/ab3105/rds/hpc-work/ukbiobank/obesity/subset/pgwas${id}.tab  --hwe 1e-10 --maf 0.001 --geno .05 --linear hide-covar --pheno /home/ab3105/rds/hpc-work/ukbiobank/obesity/obesity_phenotypes.tab --covar /home/ab3105/rds/hpc-work/ukbiobank/obesity/obesity_covariates.tab --out /home/ab3105/rds/hpc-work/ukbiobank/obesity/assoc/chr${chr}_psubset${id};
done

for id in {1..4}; do \
plink2 --bfile /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2  --fam /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2_s488131.fam --bim /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr${chr}_v2.bim --keep /home/ab3105/rds/hpc-work/ukbiobank/obesity/subset/cgwas${id}.tab  --hwe 1e-10 --maf 0.001 --geno .05 --linear hide-covar --pheno /home/ab3105/rds/hpc-work/ukbiobank/obesity/obesity_phenotypes.tab --covar /home/ab3105/rds/hpc-work/ukbiobank/obesity/obesity_covariates.tab --out /home/ab3105/rds/hpc-work/ukbiobank/obesity/assoc/chr${chr}_csubset${id};
done

plink2 --bfile /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2 --fam /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb22418_c${chr}_b0_v2_s488131.fam --bim /home/ab3105/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes/ukb_snp_chr${chr}_v2.bim --keep /home/ab3105/rds/hpc-work/ukbiobank/obesity/subset/cgwas.tab  --hwe 1e-10 --maf 0.001 --geno .05 --linear hide-covar --pheno /home/ab3105/rds/hpc-work/ukbiobank/obesity/obesity_phenotypes.tab --covar /home/ab3105/rds/hpc-work/ukbiobank/obesity/obesity_covariates.tab --out /home/ab3105/rds/hpc-work/ukbiobank/obesity/assoc/chr${chr}_csubset0;
