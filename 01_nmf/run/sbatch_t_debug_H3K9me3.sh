#!/bin/bash

sbatch -c 10 --mem=100G --partition=compute_fat -t 7-0:00 \
	-o nmf_test04.o.txt \
	-e nmf_test04.e.txt \
		debug_python.sh
