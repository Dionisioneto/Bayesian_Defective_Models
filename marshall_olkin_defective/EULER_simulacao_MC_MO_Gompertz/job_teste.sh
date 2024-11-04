#PBS -N monte_carlo_MO_Gompertz_cens45
#PBS -l ncpus=1
#PBS -l walltime=200:00:00
#PBS -m ae
#PBS -M dionisioneto@usp.br
cd /home/dionisio/EULER_simulacao_MC_MO_Gompertz
module load gcc/4.9.2
module load R/3.3.3
R CMD BATCH monte_carlo_MO_Gompertz.R