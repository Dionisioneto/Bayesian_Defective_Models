#PBS -N monte_carlo_K_Gompertz_cens15
#PBS -l ncpus=1
#PBS -l walltime=300:00:00
#PBS -m ae
#PBS -M dionisioneto@usp.br
cd /home/dionisio/EULER_simulacao_MC_K_Gompertz
module load gcc/4.9.2
module load R/3.3.3
R CMD BATCH monte_carlo_K_Gompertz.R