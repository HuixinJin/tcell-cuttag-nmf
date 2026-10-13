#!/bin/bash

sbatch -c 10 --mem=400G --partition=compute_fat -t 7-0:00 \
	-o nmf_test03.o.txt \
	-e nmf_test03.e.txt \
		do_python.sh
