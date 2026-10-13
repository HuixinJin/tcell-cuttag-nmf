#!/bin/bash

sbatch -c 5 --mem=40G --partition=compute_plus -t 7-0:00 \
	-o computematrix_test01.o.txt \
	-e computematrix_test01.e.txt \
		t_deeptools.sh
