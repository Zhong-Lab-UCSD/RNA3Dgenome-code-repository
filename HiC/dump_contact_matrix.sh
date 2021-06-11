############ juicer 
main_dir="/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/"

for sample_dir in H1_control H1_NH4OAc H1_FL H1_RNase; do
#mkdir $main_dir${sample_dir}_merged
for resolution in 2500000 5000000 1000000 500000 50000 10000; do
mkdir $main_dir${sample_dir}_merged/${resolution}
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar dump observed KR /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_$sample_dir/HiC_$sample_dir.mapq30.hic $chr $chr BP $resolution $main_dir${sample_dir}_merged/${resolution}"/chr"$chr"_chr"$chr"_"$resolution".txt"
done
done
done

sample_dir=H1_RNaseCtrl
mkdir $main_dir${sample_dir}_merged
for resolution in 2500000 5000000 1000000 500000 50000 10000; do
mkdir $main_dir${sample_dir}_merged/${resolution}
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar dump observed KR /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_$sample_dir/HiC_$sample_dir.hic $chr $chr BP $resolution $main_dir${sample_dir}_merged/${resolution}"/chr"$chr"_chr"$chr"_"$resolution".txt"
done
done


####### Interchromosomal nuclear maps
sample_dir=H1_control
resolution=500000
for chr1 in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
for chr2 in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar dump observed KR /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_$sample_dir/HiC_$sample_dir.mapq30.hic $chr1 $chr2 BP $resolution $main_dir$sample_dir"_merged/"$resolution"/chr"$chr1"_chr"$chr2"_"$resolution".txt"
done
done


######### RNase interchromosomal mitochondrial
sample_dir=H1_RNase
resolution=1000
chr1=M
chr2=17

java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar dump observed KR /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_$sample_dir/HiC_$sample_dir.mapq30.hic $chr1 $chr2 BP $resolution $main_dir$sample_dir"_merged/"$resolution"/chr"$chr1"_chr"$chr2"_"$resolution".txt"


####### Oe KR maps

main_dir="/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/oe_KR/"

for sample_dir in H1_control H1_NH4OAc H1_FL H1_RNase
do
#mkdir $main_dir${sample_dir}_merged
for resolution in 2500000 5000000 1000000 500000; do
mkdir $main_dir${sample_dir}_merged/${resolution}
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar dump oe KR /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_$sample_dir/HiC_$sample_dir.mapq30.hic $chr $chr BP $resolution $main_dir${sample_dir}_merged/${resolution}"/chr"$chr"_chr"$chr"_"$resolution".txt"
done
done
done

sample_dir=H1_RNaseCtrl
mkdir $main_dir${sample_dir}_merged
for resolution in 2500000 5000000 1000000 500000; do
mkdir $main_dir${sample_dir}_merged/${resolution}
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar dump oe KR /dataOS/frankyan/phase_separation/std_results/merged_HiC/HiC_$sample_dir/HiC_$sample_dir.hic $chr $chr BP $resolution $main_dir${sample_dir}_merged/${resolution}"/chr"$chr"_chr"$chr"_"$resolution".txt"
done
done


########################### Dump matrices at the biological replicate level
main_dir="/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/KR/"

for sample_dir in H1_control H1_NH4OAc H1_FL H1_RNase; do
	for biological_rep in 1 2; do
	mkdir $main_dir${sample_dir}_$biological_rep

	for resolution in 5000000 1000000 500000 50000 10000; do
		mkdir $main_dir${sample_dir}_$biological_rep/${resolution}

		for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do

			java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar dump observed KR /dataOS/frankyan/phase_separation/std_results/HiC/HiC_${sample_dir}_$biological_rep/strict_MAPQ30/HiC_${sample_dir}_$biological_rep.mapq30.hic $chr $chr BP $resolution $main_dir${sample_dir}_$biological_rep/${resolution}"/chr"$chr"_chr"$chr"_"$resolution".txt"

		done
	done
	done
done


main_dir="/dataOS/rcalandrelli/phase_separation/HiC/HiC_contact_matrices/oe_KR/"

for sample_dir in H1_control H1_NH4OAc H1_FL H1_RNase; do
	for biological_rep in 1 2; do
	mkdir $main_dir${sample_dir}_$biological_rep

	for resolution in 5000000 1000000 500000 50000 10000; do
		mkdir $main_dir${sample_dir}_$biological_rep/${resolution}

		for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do

			java -Xms512m -Xmx2048m -jar /dataOS/rcalandrelli/Software/juicer_tools_1.21.01.jar dump oe KR /dataOS/frankyan/phase_separation/std_results/HiC/HiC_${sample_dir}_$biological_rep/strict_MAPQ30/HiC_${sample_dir}_$biological_rep.mapq30.hic $chr $chr BP $resolution $main_dir${sample_dir}_$biological_rep/${resolution}"/chr"$chr"_chr"$chr"_"$resolution".txt"

		done
	done
	done
done