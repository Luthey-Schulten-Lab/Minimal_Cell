#!/bin/bash
mpirun -np 15 python3 MinCell_CMEODE_mpi_two_.py -t 1

for ((i = $1; i < 2; i++));do
	echo "restart level is:" $i
	#python3 MC_wRestart.py -procid ${SGE_TASK_ID} -t 1 -iter $i
	python3 MCrestartLoop.py -procid $2 -t 1 -iter $i
	rm ./sims/scan125-ribo/out-${i}rs${i}.lm
	rm ./sims/scan125-ribo/log-${i}rs${i}.log
	#wait
done
