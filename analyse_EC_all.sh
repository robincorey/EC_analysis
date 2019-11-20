#!/bin/bash

CD=`/sansom/s156a/bioc1535/EC_MEMPROT/Data/5us_analyses`

cattrj () {
mkdir -p $1
if [[ ! -f $1/$1.$2.xtc ]]
then
	echo -e Protein '\n' System | gmx trjcat -f ../$1/md_$1_$2*.xtc -o $1/$1.$2.xtc -cat
fi
}

initialise () {
conda env remove --name PyLipID
rm -r PyLipID
git clone https://github.com/wlsong/PyLipID
cd PyLipID
conda env create -f env.yml
}


run_lipid_analysis () {
mkdir -p Sites/$1
echo Protein | gmx editconf -f ../$1/md_$1_1.tpr.gro -o Sites/$1/prot.pdb -n ../$1/leaflets.ndx
python3 PyLipID/pylipid.py -f $1/$1.1.xtc $1/$1.2.xtc $1/$1.3.xtc $1/$1.4.xtc $1/$1.5.xtc -c ../$1/md_$1_1.tpr.gro ../$1/md_$1_1.tpr.gro ../$1/md_$1_1.tpr.gro ../$1/md_$1_1.tpr.gro ../$1/md_$1_1.tpr.gro -lipids CARD -lipid_atoms GL0 PO1 PO2 -save_dir Sites/$1/lipid_interactions -cutoffs 0.5 1 -nprot 1 -resi_offset 1 -pdb Sites/$1/prot.pdb
}

leaflet_analysis () {
for leaflet in lower upper
do
        mkdir -p Leaflets/$leaflet/$1
	rm -f Leaflets/$leaflet/$1/$leaflet*.txt
	for lipid in POPE POPG CARD
	do
		echo 'group "'$lipid'_'$leaflet'" and same residue as within 0.6 of group Protein' > Leaflets/$leaflet/select_$lipid.dat
		gmx "select" -f $1/$1.$2.xtc -s ../$1/md_$1_$2.tpr -os Leaflets/$leaflet/$1/$lipid.$2.xvg -sf Leaflets/$leaflet/select_$lipid.dat -n ../$1/leaflets.ndx -tu ns -b 500 -seltype res_cog -pbc -dt 10
	done
done
}

plot_all () {
python EC_analysis/PLOT_leaflets.py
}

cd $CD

#initialise

#conda init bash
#source ~/.bashrc
#conda activate PyLipID

for pdb in 5OQT #4JR9 2HI7 3O7P 1NEK 3ZE3 5MRW 5OC0 1ZCD 1PV6 3OB6 4IU8 1KF6 4KX6 3QE7 1KPK 5SV0 1U77 1FFT 5AJI 4ZP0 5AZC 1Q16 1FX8 6AL2 2QFI 3DHW 1IWG 3K07 4GD3 5JWY 3B5D 4Q65 1PW4 2R6G 4DJI 2WSX 1NEK 1L7V 1RC2 5ZUG 2IC8
do
	for num in 1 2 3 4 5 
	do
		cattrj $pdb $num
		leaflet_analysis $pdb $num
	done
        run_lipid_analysis $pdb
done

plot_all
