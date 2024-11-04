#PBS -N Teste_programaR
#PBS -l ncpus=1
#PBS -l walltime=96:00:00
#PBS -m ae
#PBS -M dionisioneto@usp.br
cd /home/dionisio/teste_simulacao_EULER
module load gcc/4.9.2
module load R/3.3.3
R CMD BATCH monte_carlo_MO_Gompertz.R