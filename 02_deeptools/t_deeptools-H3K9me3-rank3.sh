#!/bin/bash

export PATH=/dshare/xielab/home/jinhx/miniconda3/envs/deeptools/bin:$PATH

h=$1
cluster=$2
h_dir=H3K9me3_rank3
bed_dir=/dshare/xielab/analysisdata/ThreeDGene/zhangjk/TCellProject/CUT_TAG/pool/
num_processors=$(nproc)

mkdir -p deeptools/${h_dir}

computeMatrix reference-point \
	--referencePoint center \
	-b 2000 -a 2000 \
	-R t/result/${h_dir}/cluster_bed/${cluster}_peaks.bed \
	-S ${bed_dir}/Naive_1-Naive_${h}/Naive_1-Naive_${h}.cpm.norm.bw ${bed_dir}/AVI-d8_3-Arm-MP_${h}/AVI-d8_3-Arm-MP_${h}.cpm.norm.bw ${bed_dir}/AVI-d70_4-Arm-Tcm_${h}/AVI-d70_4-Arm-Tcm_${h}.cpm.norm.bw ${bed_dir}/AVI-d70_6-Arm-Trm_${h}/AVI-d70_6-Arm-Trm_${h}.cpm.norm.bw ${bed_dir}/AVI-d70_5-Arm-Tem_${h}/AVI-d70_5-Arm-Tem_${h}.cpm.norm.bw ${bed_dir}/AVI-d8_2-Arm-Teff_${h}/AVI-d8_2-Arm-Teff_${h}.cpm.norm.bw ${bed_dir}/B16GP-TdLN_12-TdLN-Ttsm_${h}/B16GP-TdLN_12-TdLN-Ttsm_${h}.cpm.norm.bw ${bed_dir}/B16GP-TdLN_13-TdLN-Tpex_${h}/B16GP-TdLN_13-TdLN-Tpex_${h}.cpm.norm.bw ${bed_dir}/TME_14-TME-Tpex_${h}/TME_14-TME-Tpex_${h}.cpm.norm.bw ${bed_dir}/TME_15-TME-Ttran_${h}/TME_15-TME-Ttran_${h}.cpm.norm.bw ${bed_dir}/TME_16-TME-Tex_${h}/TME_16-TME-Tex_${h}.cpm.norm.bw ${bed_dir}/CVI-d21_9-Cl-13-Tpex_${h}/CVI-d21_9-Cl-13-Tpex_${h}.cpm.norm.bw ${bed_dir}/CVI-d7_8-Cl-13-Tpex_${h}/CVI-d7_8-Cl-13-Tpex_${h}.cpm.norm.bw ${bed_dir}/CVI-d21_10-Cl-13-Ttran_${h}/CVI-d21_10-Cl-13-Ttran_${h}.cpm.norm.bw ${bed_dir}/CVI-d21_11-Cl-13-Tex_${h}/CVI-d21_11-Cl-13-Tex_${h}.cpm.norm.bw ${bed_dir}/CVI-d7_7-Cl-13-Tef_${h}/CVI-d7_7-Cl-13-Tef_${h}.cpm.norm.bw \
	--skipZeros \
	-p $num_processors \
	-o deeptools/${h_dir}/matrix_${h_dir}_${cluster}_center.gz \
	--outFileSortedRegions deeptools/${h_dir}/regions_${h_dir}_${cluster}_center.bed

plotHeatmap -m deeptools/${h_dir}/matrix_${h_dir}_${cluster}_center.gz \
	--samplesLabel "Naive" "Arm-MP" "Arm-Tcm" "Arm-Trm" "Arm-Tem" "Arm-Teff" "TdLN-Ttsm" "TdLN-Tpex" "TME-Tpex" "TME-Ttran" "TME-Tex" "Cl-Tpex" "Cl-Tpex" "Cl-Tran" "Cl-Tex" "Cl-Teff" \
	--plotTitle ${h_dir}_${cluster} \
	-x "" \
	--heatmapHeight 15 \
	--regionsLabel peaks \
	--colorMap coolwarm \
	-out deeptools/${h_dir}/heatmap_${h_dir}_${cluster}_center.pdf
