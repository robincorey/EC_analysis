#!/bin/bash

CD=/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses
dir=$CD/PyLipID_poses2/Sites/
mkdir -p $dir

for pdb in 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
	data=$CD/PyLipID_poses2/$pdb/lipid_interactions/Interaction_CARD/Binding_Sites_CARD/
	cp $data/BindingSites_Info_CARD.txt $dir/${pdb}_info.txt
	for frame in `ls $data/Binding_Poses/*0.gro`
	do
		num=`echo $frame | awk -F '_No0' '{print $1}' | awk -F 'BSid' '{print $2}'`
		gmx editconf -f $frame -o $dir/${pdb}_site$num.pdb
	done
	rm -f $dir/*#*
done
