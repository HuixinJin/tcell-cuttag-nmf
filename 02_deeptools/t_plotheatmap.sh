#!/bin/bash

export PATH=/dshare/xielab/home/jinhx/miniconda3/envs/deeptools/bin:$PATH

plotHeatmap -m t/result/H3K27ac/matrix_H3K27ac_cluster_2_center.gz \
	--samplesLabel "Arm-Tcm" "Arm-Tem" "Arm-Trm" "Arm-Teff" "Arm-MP" "TdLN-Ttsm" "TdLN-Tpex" "Cl-Tran" "Cl-Tex" "Cl-Tpex" "Cl-Teff" "Cl-Tpex" "Naive" "TME-Tpex" "TME-Ttran" "TME-Tex" \
	--plotTitle H3K27ac_cluster_2 \
	-x "" \
	--heatmapHeight 15 \
	--regionsLabel peaks \
	-out t/result/H3K27ac/heatmap_H3K27ac_cluster_2_center.pdf
