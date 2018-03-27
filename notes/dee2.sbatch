#!/bin/bash
# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'

# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=mmu_chn2

# To set a project account for credit charging, 
#SBATCH --account=bd17

# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12

# Memory usage (MB)
#SBATCH --mem-per-cpu=4000

# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=0-4:00:00

# SBATCH --partition=m3a
# SBATCH --reservation=epigenetics

# Set the file for output (stdout)
#SBATCH --output=/projects/bd17/mziemann/dee2/MyJob-%j.out

# Set the file for error log (stderr)
#SBATCH --error=/projects/bd17/mziemann/dee2/MyJob-%j.err

set -x

clear_shared_mem(){
for SHMEM_ID in $(ipcs | awk '$5>1000 && $3~/mziemann/ {print $2}') ; do
 ipcrm --shmem-id $SHMEM_ID
done
}
export -f clear_shared_mem

clear_shared_mem

module add singularity/2.4.2

cd /scratch/bd17/mziemann/singularity/working
touch a

for i in * ; do
  CONTAINER_AGE=$(expr $(date +%s) - $(stat -c %Y $i) )
  echo $i modified $CONTAINER_AGE seconds ago 
  if [ $CONTAINER_AGE -gt 14400 ] ; then
    rm -rf $i
    echo deleted $i
  fi
done


TIMESTAMP="$(date +'%s')_$(echo {a..z}{0..9} | tr ' ' '\n'  | shuf -n1)"
CWD=dee2_container_${TIMESTAMP}
mkdir $CWD
cd $CWD
cp -r /projects/bd17/mziemann/dee2/sing/singularity_mmusculus/* .
pwd

while true ; do
 clear_shared_mem
 singularity run -w -B tmp:/tmp mziemann_tallyup-2018-01-23-3a04c7fd01c6.img mmusculus
# udocker run --rm mziemann/tallyup mmusculus
 sleep 10
done &

sleep 110m
killall -9 udocker
clear_shared_mem
sbatch /projects/bd17/mziemann/dee2/dee2_scratch_mmu_sing.sbatch

