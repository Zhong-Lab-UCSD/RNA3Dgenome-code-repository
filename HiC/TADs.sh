############ juicer
main_dir="/mnt/extraids/OceanStor-SysCmn-2/rcalandrelli/phase_separation/HiC/tads/"
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

for i in "${chromosomes[@]}"; do
	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar arrowhead \
	-c chr$i \
	-m 2000 \
	-r 50000 \
	--threads 40 \
	-k KR \
	/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
	$main_dir"H1_control_merged/TADs_chr"$i

	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar arrowhead \
	-c chr$i \
	-m 2000 \
	-r 10000 \
	--threads 40 \
	-k KR \
	/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
	$main_dir"H1_control_merged/TADs_chr"$i
done

for i in "${chromosomes[@]}"; do
	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar arrowhead \
	-c chr$i \
	-m 2000 \
	-r 50000 \
	--threads 40 \
	-k KR \
	/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
	$main_dir"H1_NH4OAc_merged/TADs_chr"$i

	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar arrowhead \
	-c chr$i \
	-m 2000 \
	-r 10000 \
	--threads 40 \
	-k KR \
	/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
	$main_dir"H1_NH4OAc_merged/TADs_chr"$i
done

for i in "${chromosomes[@]}"; do
	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar arrowhead \
	-c chr$i \
	-m 2000 \
	-r 50000 \
	--threads 40 \
	-k KR \
	/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
	$main_dir"H1_FL_merged/TADs_chr"$i

	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar arrowhead \
	-c chr$i \
	-m 2000 \
	-r 10000 \
	--threads 40 \
	-k KR \
	/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
	$main_dir"H1_FL_merged/TADs_chr"$i
done

for i in "${chromosomes[@]}"; do
	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar arrowhead \
	-c chr$i \
	-m 2000 \
	-r 50000 \
	--threads 40 \
	-k KR \
	/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
	$main_dir"H1_RNase_merged/TADs_chr"$i

	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar arrowhead \
	-c chr$i \
	-m 2000 \
	-r 10000 \
	--threads 40 \
	-k KR \
	/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
	$main_dir"H1_RNase_merged/TADs_chr"$i
done