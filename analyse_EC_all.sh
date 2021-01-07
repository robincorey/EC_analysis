#!/bin/bash

# dir with data present
CD="/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses"

# define functions here:

# simple function to combine data into single trajectory
cattrj () {
mkdir -p $1
if [[ ! -f $1/$1.$2.xtc ]]
then
	echo -e Protein '\n' System | gmx trjcat -f ../$1/md_$1_$2*.xtc -o $1/$1.$2.xtc -dt 1000 #-cat
fi
echo -e Protein '\n' System | gmx trjconv -f $1/$1.$2.xtc -s ../$1/md_$1_$2.tpr -o $1/$1.$2.pbc.xtc -pbc mol -center
}

# function for analysis using PyLipID
# see https://github.com/wlsong/PyLipID

run_lipid_analysis () {
mkdir -p PyLipID_poses2/$1
for i in {1..5}; do gmx editconf -f ../$1/md_$1_$i.tpr -o ../$1/md_$1_$i.tpr.gro ; done
python3 $GIT/PyLipID/pylipid.py -f $1/$1.{1..5}.pbc.xtc -c ../$1/md_$1_{1..5}.tpr.gro -lipids CARD -lipid_atoms GL0 PO1 PO2 -save_dir PyLipID_poses2/$1/lipid_interactions -cutoffs 0.55 1 -nprot 1 -pymol_gui False -gen_binding_poses 10 -score_weights GL0:10 PO1:10 PO2:10 -stride 10
}

# analysis of lipid contacts across the leaflets
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

# not used in final study
phi_analysis () {
mkdir -p Phi_analysis
mkdir -p Phi_analysis/$1
for lipid in CARD POPE POPG
do
	python lipid-contact_phi.py ../$1/md_$1_$2.tpr $1/$1.$2.xtc Phi_analysis/$1/$1.$2.$lipid $lipid
done
}

# not used in final study
phi_combine () {
mkdir -p Phi_analysis/res
cd Phi_analysis
for lipid in CARD #POPG POPE
do
        rm -f res/$lipid.$1.txt
        for res in ARG LYS ALA ILE VAL LEU TRP PHE TYR HIS GLY PRO CYS MET ASP GLU ASN GLN SER THR
        do
                tot=0
                for i in 1 2 3 4 5
                do
                        val=`grep $res $1/$1.$i.$lipid.csv | awk -F, '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}'`
                        tot=`echo "$tot + $val" | bc`
                done
                ave=`echo "scale=4; $tot / 5" | bc`
                echo $res,$ave >> res/$lipid.$1.txt
        done
        #paste -d ',' */*.$i.$lipid.csv | awk -F, '{for(i=1;i<=NF;i++) t+=$i; print $1","t/5; t=0}' > $lipid.$i.txt
done
#paste -d ',' CARD.$i.txt POPE.$i.txt POPG.$i.txt > LIP.$i.txt 
}

# generation of rough site
site_predict_CARD () {
mkdir -p Predicting_site_CARD/raw
for component in 'GL0' 'PO2 PO1' 'GL2 GL1 GL3 GL4' 'C* D*'
do
	name=`echo ${component} | tr -d " "*`
	python ../lipid-contact_construct_CARD.py ../$1/md_$1_1.tpr.gro $1/$1.$2.xtc Predicting_site_CARD/raw/$1_$2_$name "$component"
done
}

site_predict_POPE () {
mkdir -p Predicting_site_POPE/raw
for component in 'NH3' 'PO4' 'GL2 GL1 GL3' 'C* D*'
do
        name=`echo ${component} | tr -d " "*`
        python ../lipid-contact_construct_POPE.py ../$1/md_$1_1.tpr.gro $1/$1.$2.xtc Predicting_site_POPE/raw/$1_$2_$name "$component"
done
}

site_predict_POPG () {
mkdir -p Predicting_site_POPG/raw
for component in 'GL0' 'PO4' 'GL2 GL1 GL3' 'C* D*' 
do
        name=`echo ${component} | tr -d " "*`
        python ../lipid-contact_construct_POPG.py ../$1/md_$1_1.tpr.gro $1/$1.$2.xtc Predicting_site_POPG/raw/$1_$2_$name "$component"
	#paste -d ',' ${pdb}_*_${lipid}.csv | awk -F, '{for(i=1;i<=NF;i++) t+=$i; print $1","t/5; t=0}' > ${pdb}_$lipid.txt
done
}

# residue distribution across membrane
residue_distribution () {
mkdir -p Residue_distribution/$1
echo -e aPO4 '\n' q | gmx make_ndx -f ../${1}/md_${1}_1.tpr -o Residue_distribution/$1/$1_po4.ndx
echo -e PO4 '\n' System | gmx trjconv -f ${1}/${1}.${2}.xtc -s ../${1}/md_${1}_${2}.tpr -boxcenter zero -tu us -skip 1000 -o Residue_distribution/$1/md_${1}_${2}.po4.xtc -n Residue_distribution/$1/${1}_po4.ndx -center
echo -e System | gmx trjconv -f Residue_distribution/$1/md_${1}_${2}.po4.xtc -s ../${1}/md_${1}_${2}.tpr -tu us -dump 1 -b 0 -o Residue_distribution/$1/md_${1}_${2}.po4.pdb
python ../res_contact_z.py Residue_distribution/$1/md_${1}_${2}.po4.pdb Residue_distribution/$1/md_${1}_${2}.po4.xtc Residue_distribution/$1/md_${1}_${2}
awk -F, '{print $1","($2+$3+$4)}' Residue_distribution/$1/md_${1}_${2}.csv > Residue_distribution/$1/md_${1}_${2}_comb.csv
perl /sansom/s137/bioc1535/Desktop/CG_KIT/make_b_factor.pl Residue_distribution/$1/md_${1}_${2}.po4.pdb Residue_distribution/$1/md_${1}_${2}_comb.csv 1 Residue_distribution/$1/z_${1}_${2}.pdb
grep -e BB -e SC* Residue_distribution/$1/z_${1}_${2}.pdb |  sed 's/[0-9]-/[0-9] /g' | awk '{print $4","$8","$10}' > Residue_distribution/z_${1}_${2}.pdb
grep -e BB -e SC* Residue_distribution/$1/md_${1}_${2}.po4.pdb |  sed 's/[0-9]-/[0-9] /g' | awk '{print $4","$8","$10}' > Residue_distribution/all_${1}_${2}.pdb
rm -f Residue_distribution/$1/z_${1}_${2}.pdb
}

residue_distribution_CDL () {
# CDL GL0+PO*:Arg/Lys and GL0+PO*:Ser/Gly/Thr contacts along Z axis
mkdir -p Residue_distribution_CDL/$1
# input unchaged, output to new file
###python /sansom/s156a/bioc1535/EC_MEMPROT/res_contact_z.cdl.py Residue_distribution/$1/md_${1}_${2}.po4.pdb Residue_distribution/$1/md_${1}_${2}.po4.xtc Residue_distribution_CDL/$1/md_${1}_${2}
awk -F, '{print $1","($2+$3+$4)}' Residue_distribution_CDL/$1/md_${1}_${2}.csv > Residue_distribution_CDL/$1/md_${1}_${2}_comb.csv
perl /sansom/s137/bioc1535/Desktop/CG_KIT/make_b_factor.pl Residue_distribution/$1/md_${1}_${2}.po4.pdb Residue_distribution_CDL/$1/md_${1}_${2}_comb.csv 1 Residue_distribution_CDL/$1/z_${1}_${2}.pdb
grep -e BB -e SC* Residue_distribution_CDL/$1/z_${1}_${2}.pdb |  sed 's/[0-9]-/[0-9] /g' | awk '{print $4","$8","$10}' > Residue_distribution_CDL/z_${1}_${2}.pdb
#grep -e BB -e SC* Residue_distribution_CDL/$1/md_${1}_${2}.po4.pdb |  sed 's/[0-9]-/[0-9] /g' | awk '{print $4","$8","$10}' > Residue_distribution_CDL/all_${1}_${2}.pdb
rm -f Residue_distribution_CDL/$1/z_${1}_${2}.pdb
}

# PLOT final data
plot_all () {
#python EC_analysis/PLOT_leaflets.py
#python EC_analysis/PLOT_predicted_site.py
#python EC_analysis/PLOT_predicted_site_allres_CARD.py
#python EC_analysis/PLOT_predicted_site_allres_POPE.py
#python EC_analysis/PLOT_predicted_site_allres_POPG.py
#python EC_analysis/PLOT_z_analysis.py ARG LYS
#python EC_analysis/PLOT_z_analysis_allres.py ARG LYS
#python EC_analysis/PLOT_z_analysis.py ASP GLU ASN GLN CYS
#python EC_analysis/PLOT_z_analysis.py ILE LEU VAL ALA MET
#python EC_analysis/PLOT_z_analysis.py SER THR GLY PRO 
#python EC_analysis/PLOT_z_analysis.py TRP TYR PHE HIS
#python EC_analysis/PLOT_phi.py
python EC_analysis/PLOT_z_analysis.cdl.py ARG LYS GLY SER THR HIS
}

cd $CD
# loop through all PDBs analyses
for pdb in 1NEK 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
	for num in 1 2 3 4 5 
	do
		#cattrj $pdb $num
		#leaflet_analysis $pdb $num
#		site_predict_CARD $pdb $num
#		site_predict_POPE $pdb $num
#		site_predict_POPG $pdb $num
#		residue_distribution $pdb $num
		residue_distribution_CDL $pdb $num	
		#phi_analysis $pdb $num		
	done
	cd $CD
#	rm -f Phi_analysis/res/$pdb*txt
#	phi_combine $pdb
#	run_lipid_analysis $pdb
done

#phi_combine
cd $CD
plot_all
