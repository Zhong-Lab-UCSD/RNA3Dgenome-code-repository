############ juicer
main_dir="/dataOS/rcalandrelli/phase_separation/HiC/compartments/"
chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

mkdir $main_dir"H1_control_merged/2500000"
mkdir $main_dir"H1_control_merged/1000000"
mkdir $main_dir"H1_control_merged/500000"
for i in "${chromosomes[@]}"; do
	
    java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	KR \
	/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
	chr$i \
	BP \
	2500000 \
	$main_dir"H1_control_merged/2500000/eigen_chr"$i".txt"

    # java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	# KR \
	# /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
	# chr$i \
	# BP \
	# 1000000 \
	# $main_dir"H1_control_merged/1000000/eigen_chr"$i".txt"
    
    # java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	# KR \
	# /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
	# chr$i \
	# BP \
	# 500000 \
	# $main_dir"H1_control_merged/500000/eigen_chr"$i".txt"
done

mkdir $main_dir"H1_NH4OAc_merged/2500000"
mkdir $main_dir"H1_NH4OAc_merged/1000000"
mkdir $main_dir"H1_NH4OAc_merged/500000"
for i in "${chromosomes[@]}"; do

    java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	KR \
	/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
	chr$i \
	BP \
	2500000 \
	$main_dir"H1_NH4OAc_merged/2500000/eigen_chr"$i".txt"

    # java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	# KR \
	# /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
	# chr$i \
	# BP \
	# 1000000 \
	# $main_dir"H1_NH4OAc_merged/1000000/eigen_chr"$i".txt"

    # java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	# KR \
	# /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
	# chr$i \
	# BP \
	# 500000 \
	# $main_dir"H1_NH4OAc_merged/500000/eigen_chr"$i".txt"
done

mkdir $main_dir"H1_FL_merged/2500000"
mkdir $main_dir"H1_FL_merged/1000000"
mkdir $main_dir"H1_FL_merged/500000"
for i in "${chromosomes[@]}"; do
	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	KR \
	/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
	chr$i \
	BP \
	2500000 \
	$main_dir"H1_FL_merged/2500000/eigen_chr"$i".txt"

	# java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	# KR \
	# /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
	# chr$i \
	# BP \
	# 1000000 \
	# $main_dir"H1_FL_merged/1000000/eigen_chr"$i".txt"

    # java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	# KR \
	# /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
	# chr$i \
	# BP \
	# 500000 \
	# $main_dir"H1_FL_merged/500000/eigen_chr"$i".txt"
done

mkdir $main_dir"H1_RNase_merged/2500000"
mkdir $main_dir"H1_RNase_merged/1000000"
mkdir $main_dir"H1_RNase_merged/500000"
for i in "${chromosomes[@]}"; do
	
	java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	KR \
	/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
	chr$i \
	BP \
	2500000 \
	$main_dir"H1_RNase_merged/2500000/eigen_chr"$i".txt"

	# java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	# KR \
	# /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
	# chr$i \
	# BP \
	# 1000000 \
	# $main_dir"H1_RNase_merged/1000000/eigen_chr"$i".txt"

    # java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
	# KR \
	# /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
	# chr$i \
	# BP \
	# 500000 \
	# $main_dir"H1_RNase_merged/500000/eigen_chr"$i".txt"
done


mkdir $main_dir"H1_RNaseCtrl_merged/2500000"
mkdir $main_dir"H1_RNaseCtrl_merged/1000000"
mkdir $main_dir"H1_RNaseCtrl_merged/500000"
for i in "${chromosomes[@]}"; do
java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
KR \
/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNaseCtrl/HiC_H1_RNaseCtrl.hic \
chr$i \
BP \
2500000 \
$main_dir"H1_RNaseCtrl_merged/2500000/eigen_chr"$i".txt"

java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
KR \
/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNaseCtrl/HiC_H1_RNaseCtrl.hic \
chr$i \
BP \
1000000 \
$main_dir"H1_RNaseCtrl_merged/1000000/eigen_chr"$i".txt"

java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar eigenvector \
KR \
/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNaseCtrl/HiC_H1_RNaseCtrl.hic \
chr$i \
BP \
500000 \
$main_dir"H1_RNaseCtrl_merged/500000/eigen_chr"$i".txt"
done