#!/bin/bash

CD=`pwd`

residue_distribution () {
mkdir -p Residue_distribution/$1
# make ndx for PO4 only
echo -e aPO4 '\n' q | gmx make_ndx -f ../${1}/md_${1}_1.tpr -o Residue_distribution/$1/$1_po4.ndx
# align traj by PO4 beads
echo -e PO4 '\n' System | gmx trjconv -f ${1}/${1}.${2}.xtc -s ../${1}/md_${1}_${2}.tpr -boxcenter zero -tu us -skip 1000 -o Residue_distribution/$1/md_${1}_${2}.po4.xtc -n Residue_distribution/$1/${1}_po4.ndx -center
# dump 1 µs frame
echo -e System | gmx trjconv -f Residue_distribution/$1/md_${1}_${2}.po4.xtc -s ../${1}/md_${1}_${2}.tpr -tu us -dump 1 -b 0 -o Residue_distribution/$1/md_${1}_${2}.po4.pdb
# run contact analysis script
python ../res_contact_z.py Residue_distribution/$1/md_${1}_${2}.po4.pdb Residue_distribution/$1/md_${1}_${2}.po4.xtc Residue_distribution/$1/md_${1}_${2}
# combine lipids into single array (residual from diff analysis)
awk -F, '{print $1","($2+$3+$4)}' Residue_distribution/$1/md_${1}_${2}.csv > Residue_distribution/$1/md_${1}_${2}_comb.csv
# convert contacts to PDB - only residues which contact lipids
perl /sansom/s137/bioc1535/Desktop/CG_KIT/make_b_factor.pl Residue_distribution/$1/md_${1}_${2}.po4.pdb Residue_distribution/$1/md_${1}_${2}_comb.csv 1 Residue_distribution/$1/z_${1}_${2}.pdb
# extract z coords of residues
grep -e BB -e SC* Residue_distribution/$1/z_${1}_${2}.pdb |  sed 's/[0-9]-/[0-9] /g' | awk '{print $4","$8","$10}' > Residue_distribution/z_${1}_${2}.pdb
rm -f Residue_distribution/$1/z_${1}_${2}.pdb
}

cd $CD
# loop through all PDBs 
for pdb in .... 
do
        for num in 1 2 3 4 5
        do
               residue_distribution $pdb $num
	done
done
