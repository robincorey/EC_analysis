#!/bin/bash

CD="/sansom/s156a/bioc1535/EC_MEMPROT/Data/5us_analyses"

cattrj () {
# combine all traj
mkdir -p $1
if [[ ! -f $1/$1.$2.xtc ]]
then
	echo -e Protein '\n' System | gmx trjcat -f ../$1/md_$1_$2*.xtc -o $1/$1.$2.xtc -cat -dt 1000
fi
}

initialise () {
# for PyLipID
conda env remove --name PyLipID
rm -r PyLipID
git clone https://github.com/wlsong/PyLipID
cd PyLipID
conda env create -f env.yml
}


run_lipid_analysis () {
# do PyLipID for site
#initialise
#conda init bash
#source ~/.bashrc
#conda activate PyLipID
mkdir -p Sites/$1
gmx editconf -f ../$1/md_$1_1.tpr -o ../$1/md_$1_1.tpr.gro
echo Protein | gmx editconf -f ../$1/md_$1_1.tpr.gro -o Sites/$1/prot.pdb -n ../$1/leaflets.ndx
python3 PyLipID/pylipid.py -f $1/$1.1.xtc $1/$1.2.xtc $1/$1.3.xtc $1/$1.4.xtc $1/$1.5.xtc -c ../$1/md_$1_1.tpr.gro ../$1/md_$1_1.tpr.gro ../$1/md_$1_1.tpr.gro ../$1/md_$1_1.tpr.gro ../$1/md_$1_1.tpr.gro -lipids CARD -lipid_atoms GL0 PO1 PO2 -save_dir Sites/$1/lipid_interactions -cutoffs 0.5 1 -nprot 1 -resi_offset 1 -pdb Sites/$1/prot.pdb
}

leaflet_analysis () {
# get data for leaflet boxplot
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

phi_analysis () {
# get data for phi BS
mkdir -p Phi_analysis
mkdir -p Phi_analysis/$1
for lipid in CARD POPE POPG
do
	python lipid-contact_phi.py ../$1/md_$1_$2.tpr $1/$1.$2.xtc Phi_analysis/$1/$1.$2.$lipid $lipid
done

}


kinetics () {
# convergence analysis
:
}

site_predict () {
# get coarse site
mkdir -p Predicting_site/raw
mkdir -p Predicting_site/res
for component in 'GL0' 'PO2 PO1' 'GL2 GL1 GL3 GL4' 'C* D*'
do
	name=`echo ${component} | tr -d " "*`
	python ../lipid-contact_construct.py ../$1/md_$1_1.tpr.gro $1/$1.$2.xtc Predicting_site/raw/$1_$2_$name "$component"
#paste -d ',' ${pdb}_*_${lipid}.csv | awk -F, '{for(i=1;i<=NF;i++) t+=$i; print $1","t/5; t=0}' > ${pdb}_$lipid.txt
done
}

residue_distribution () {
# residue distribution across membrane
mkdir -p Residue_distribution/$1
echo -e aPO4 '\n' q | gmx make_ndx -f ../${1}/md_${1}_1.tpr -o Residue_distribution/$1/$1_po4.ndx
echo -e PO4 '\n' System | gmx trjconv -f ${1}/${1}.${2}.xtc -s ../${1}/md_${1}_${2}.tpr -boxcenter zero -b 0.5 -tu us -skip 1000 -o Residue_distribution/$1/md_${1}_${2}.po4.xtc -n Residue_distribution/$1/${1}_po4.ndx -center
echo -e System | gmx trjconv -f Residue_distribution/$1/md_${1}_${2}.po4.xtc -s ../${1}/md_${1}_${2}.tpr -tu us -dump 5 -b 4.8 -o Residue_distribution/$1/md_${1}_${2}.po4.pdb
python ../res_contact_z.py Residue_distribution/$1/md_${1}_${2}.po4.pdb Residue_distribution/$1/md_${1}_${2}.po4.xtc Residue_distribution/$1/md_${1}_${2}
awk -F, '{print $1","($2+$3+$4)}' Residue_distribution/$1/md_${1}_${2}.csv > Residue_distribution/$1/md_${1}_${2}_comb.csv
perl /sansom/s137/bioc1535/Desktop/CG_KIT/make_b_factor.pl Residue_distribution/$1/md_${1}_${2}.po4.pdb Residue_distribution/$1/md_${1}_${2}_comb.csv 1 Residue_distribution/$1/z_${1}_${2}.pdb
grep -e BB -e SC* Residue_distribution/$1/z_${1}_${2}.pdb | awk '{print $4","$8","$10}' > Residue_distribution/z_${1}_${2}.pdb
rm -f Residue_distribution/$1/z_${1}_${2}.pdb
}

phi_combine () {
cd Phi_analysis
for i in 1 2 3 4 5
do 
	for lipid in CARD POPG POPE
	do
		paste -d ',' */*.$i.$lipid.csv | awk -F, '{for(i=1;i<=NF;i++) t+=$i; print $1","t/5; t=0}' > $lipid.$i.txt
	done
	paste -d ',' CARD.$i.txt POPE.$i.txt POPG.$i.txt > LIP.$i.txt 
done


}

plot_all () {
python EC_analysis/PLOT_leaflets.py
python EC_analysis/PLOT_predicted_site.py
python EC_analysis/PLOT_z_analysis.py ARG LYS HIS
python EC_analysis/PLOT_z_analysis.py ASP GLU
python EC_analysis/PLOT_z_analysis.py PHE ILE LEU VAL ALA
python EC_analysis/PLOT_z_analysis.py ASN SER GLN THR
python EC_analysis/PLOT_z_analysis.py CYS GLY PRO MET
python EC_analysis/PLOT_z_analysis.py TRP TYR
python EC_analysis/PLOT_phi.py
}

cd $CD

#initialise

for pdb in 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
#for pdb in 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
	for num in 1 2 3 4 5 
	do
#		cattrj $pdb $num
		#leaflet_analysis $pdb $num
		#site_predict $pdb $num
		#residue_distribution $pdb $num
		: #		phi_analysis $pdb $num		
	done
	cd $CD
 	echo $pdb
	#phi_combine $pdb
#	run_lipid_analysis $pdb
done

phi_combine

cd $CD
#plot_all
