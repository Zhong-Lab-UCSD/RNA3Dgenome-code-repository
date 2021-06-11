############ juicer
main_dir="/dataOS/rcalandrelli/phase_separation/HiC/loops/"

# running hiccups --cpu for loops
# java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar hiccups \
# --cpu \
# -m 512 \
# -r 5000,10000 \
# -k KR \
# -f .1,.1 \
# -p 4,2 \
# -i 7,5 \
# -t 0.02,1.5,1.75,2 \
# -d 20000,20000 \
# /mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
# $main_dir"hiccups_H1_control_merged" \
# --threads 40

# java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar hiccups \
# --cpu \
# -m 512 \
# -r 5000,10000 \
# -k KR \
# -f .1,.1 \
# -p 4,2 \
# -i 7,5 \
# -t 0.02,1.5,1.75,2 \
# -d 20000,20000 \
# /mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
# $main_dir"hiccups_H1_NH4OAc_merged" \
# --threads 40

# java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar hiccups \
# --cpu \
# -m 512 \
# -r 5000,10000 \
# -k KR \
# -f .1,.1 \
# -p 4,2 \
# -i 7,5 \
# -t 0.02,1.5,1.75,2 \
# -d 20000,20000 \
# /mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
# $main_dir"hiccups_H1_FL_merged" \
# --threads 40

# java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar hiccups \
# --cpu \
# -m 512 \
# -r 5000,10000 \
# -k KR \
# -f .1,.1 \
# -p 4,2 \
# -i 7,5 \
# -t 0.02,1.5,1.75,2 \
# -d 20000,20000 \
# /mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
# $main_dir"hiccups_H1_RNase_merged" \
# --threads 40

java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar hiccups \
--cpu \
-m 512 \
-r 5000,10000 \
-k KR \
-f .1,.1 \
-p 4,2 \
-i 7,5 \
-t 0.02,1.5,1.75,2 \
-d 20000,20000 \
/dataOS/wenxingzhao/database/4DN/Micro-C/4DNFI2TK7L2F.hic \
$main_dir"hiccups_H1_ES_MicroC" \
--threads 40

java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar hiccups \
--cpu \
-m 512 \
-r 5000,10000 \
-k KR \
-f .1,.1 \
-p 4,2 \
-i 7,5 \
-t 0.02,1.5,1.75,2 \
-d 20000,20000 \
/dataOS/wenxingzhao/database/4DN/Micro-C/H1-endoderm/4DNFIPIP38S7.hic \
$main_dir"hiccups_H1_endoderm_MicroC" \
--threads 40

java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar hiccups \
--cpu \
-m 512 \
-r 5000,10000 \
-k KR \
-f .1,.1 \
-p 4,2 \
-i 7,5 \
-t 0.02,1.5,1.75,2 \
-d 20000,20000 \
/dataOS/rcalandrelli/phase_separation/HiC/4DNFIMROE6N4_HFFc6.hic \
$main_dir"hiccups_HFF" \
--threads 40

java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar hiccups \
--cpu \
-m 512 \
-r 5000,10000 \
-k KR \
-f .1,.1 \
-p 4,2 \
-i 7,5 \
-t 0.02,1.5,1.75,2 \
-d 20000,20000 \
/dataOS/rcalandrelli/phase_separation/HiC/4DNFITUOMFUQ_K562.hic \
$main_dir"hiccups_K562" \
--threads 40


### Aggregative peak analysis (APA)

# java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
# /mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
# /dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_control_merged/merged_loops.bedpe \
# /dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_control_merged \
# -u

# java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
# /mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
# /dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_NH4OAc_merged/merged_loops.bedpe \
# /dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_NH4OAc_merged \
# -u

# java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
# /mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
# /dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_FL_merged/merged_loops.bedpe \
# /dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_FL_merged \
# -u

# java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
# /mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
# /dataOS/rcalandrelli/phase_separation/HiC/loops/hiccups_H1_RNase_merged/merged_loops.bedpe \
# /dataOS/rcalandrelli/phase_separation/HiC/loops/APA/H1_RNase_merged \
# -u

### Aggregate peak analysis on union of loops

java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe \
/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_control_merged \
-u

java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe \
/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_NH4OAc_merged \
-u

java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe \
/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_FL_merged \
-u

java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
/mnt/extraids/OceanStor-SysCmn-5/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe \
/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/H1_RNase_merged \
-u

### Aggregate peak analysis on union of loops with n=0 (no minimum distance from the diagonal)

java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
-u \
-n 0 \
-r 10000,5000 \
/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_control/HiC_H1_control.mapq30.hic \
/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe \
/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/ALL/H1_control_merged

java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
-u \
-n 0 \
-r 10000,5000 \
/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_NH4OAc/HiC_H1_NH4OAc.mapq30.hic \
/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe \
/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/ALL/H1_NH4OAc_merged

java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
-u \
-n 0 \
-r 10000,5000 \
/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_FL/HiC_H1_FL.mapq30.hic \
/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe \
/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/ALL/H1_FL_merged

java -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar apa \
-u \
-n 0 \
-r 10000,5000 \
/dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_H1_RNase/HiC_H1_RNase.mapq30.hic \
/dataOS/rcalandrelli/phase_separation/HiC/loops/union_of_loops.bedpe \
/dataOS/rcalandrelli/phase_separation/HiC/loops/APA_union_loops/ALL/H1_RNase_merged