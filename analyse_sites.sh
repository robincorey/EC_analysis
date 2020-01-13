#!/bin/bash

get_site_info () {
:
}

strip_residues () {
:
}

get_pdb () {
:
}

cluster () {

}

rm -f sites_all.txt
count=0
for pdb in 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
	grep "BS Residence Time" $pdb/lipid_interactions/Interaction_CARD/interaction_network_CARD/BindingSites_Info_CARD.txt | while read -r line
	do
		res=`echo $line | awk '{print $4}'`
		R=`echo $line | awk '{print $9}'`
		echo $pdb $count $res $R >> sites_all.txt
		(( count++ ))
	done
	count_all=`grep -c $pdb sites_all.txt`
	echo $pdb $count_all
done
