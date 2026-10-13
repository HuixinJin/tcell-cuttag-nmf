#!/bin/bash


#for h in "H3K27ac" "H3K4me1" "H3K4me3";
for h in "H3K9me3_rank3";
do
	directory="t/result/${h}/cluster_bed"
	log_dir="deeptools/${h}/log"
	mkdir -p $log_dir

	shopt -s nullglob
	for file in "$directory"/*.bed;do
		base_name="${file##*/}"  # Get the filename without the path  
		cluster_name="${base_name%%_peaks.bed}"  # Remove the suffix "_peaks.bed"  
		
		job_id=$(sbatch -c 8 --mem=20G --partition=compute_plus -t 7-0:00 \
			-o "$log_dir/${cluster_name}.o.txt" \
			-e "$log_dir/${cluster_name}.e.txt" \
			t_deeptools-H3K9me3-rank3.sh H3K9me3 "$cluster_name")
		if [ $? -ne 0 ]; then  
                	echo "Failed to submit job for $h; $file"  
	            else  
        	        echo "Submitted job: $job_id; $h; $file"  
                fi  
	done
done
