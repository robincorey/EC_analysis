#!/bin/bash

CD="/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/Sites"
mkdir -p $CD/pdb_no_chain

for pdb in 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
	in=$CD/$pdb/lipid_interactions/Interaction_CARD/coords_Duration.pdb
	gmx editconf -f $in -o $CD/pdb_no_chain/$pdb.pdb -label ' '
	rm -f $CD/pdb_no_chain/*#*
done
