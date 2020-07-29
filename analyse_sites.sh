#!/bin/bash

# extracting site info and running analyses

# taking on residues with >=10% occ of site occ
refine_sites () {
rm -f site_refine.txt
mkdir -p refine_sites
rm -f refine_sites/site_refine.txt
for pdb in 1FFT 1FX8 1KF6 1KPK 1NEK 5OQT 4JR9 2HI7 3O7P 3ZE3 1ZCD 5OC0 1PV6 3OB6 5MRW 5AZC 1Q16 2QFI 2IC8 1RC2 1IWG 2WSX 5JWY 3B5D 3DHW 1PW4 4Q65 4DJI 2R6G 4GD3 5ZUG 6AL2 1L7V 4IU8 4KX6 3QE7 5SV0 1U77 5AJI 4ZP0 3K07 1KQF
do
        for i in {0..39}
        do
                if grep -q "Binding site $i" $pdb/lipid_interactions/Interaction_CARD/interaction_network_CARD/BindingSites_Info_CARD.txt
                then
			site_occ=`sed -n "/Binding site $i$/,/^$/p" $pdb/lipid_interactions/Interaction_CARD/interaction_network_CARD/BindingSites_Info_CARD.txt | grep "BS Lipid Occupancy" | awk '{print $4}'`
                        koff=`sed -n "/Binding site $i$/,/^$/p" $pdb/lipid_interactions/Interaction_CARD/interaction_network_CARD/BindingSites_Info_CARD.txt | grep "BS Residence Time" | awk '{print $4}'`
                        R=`sed -n "/Binding site $i$/,/^$/p" $pdb/lipid_interactions/Interaction_CARD/interaction_network_CARD/BindingSites_Info_CARD.txt | grep "BS Residence Time" | awk '{print $9}'`
			cutoff=`echo "scale=4; $site_occ / 10" | bc`
			list_prune=`sed -n "/Binding site $i$/,/^$/p" $pdb/lipid_interactions/Interaction_CARD/interaction_network_CARD/BindingSites_Info_CARD.txt | tail -n+7 | awk -v var=$cutoff '$6>var' | awk '{print $1}' | tr -d [A-Z] | sed ':a;N;$!ba;s/\n/,/g'`
			list_prune_res=`sed -n "/Binding site $i$/,/^$/p" $pdb/lipid_interactions/Interaction_CARD/interaction_network_CARD/BindingSites_Info_CARD.txt | tail -n+7 | awk -v var=$cutoff '$6>var' | awk '{print $1}' | tr -d [0-9] | sed ':a;N;$!ba;s/\n/,/g'`
                        echo $pdb $i dur=$koff R=$R occ=$site_occ resid=\"$list_prune\" res=\"$list_prune_res\"  >> refine_sites/site_refine.txt
                fi
        done
done
#sort -r -k3,3 site_resid.txt | grep -v -e dur=0.00 -e dur=nan > sites_above_10ns.txt
}

# remove sites less than 10 ns (within diffusion)
refine_stats () {
sort -r -k3,3 refine_sites/site_refine.txt | grep -v -e dur=0.00 -e dur=nan > refine_sites/sites_above_10ns.txt
# analyse pruned sites
grep -e ARG -e LYS refine_sites/sites_above_10ns.txt | awk -F',' '{print $0" resnum="NF-1}' > refine_sites/above_10ns/site_res_basic.txt
grep -e 'ARG.*LYS' -e 'LYS.*ARG' -e 'ARG.*ARG' -e 'LYS.*LYS'  refine_sites/site_refine.txt | awk -F',' '{print $0" resnum="NF-1}' > refine_sites/above_10ns/site_res_double_basic.txt
grep -v -e ARG -e LYS refine_sites/sites_above_10ns.txt | awk -F',' '{print $0" resnum="NF-1}' > refine_sites/above_10ns/site_res_no_basic.txt
awk '{print $5}' refine_sites/above_10ns/site_res_no_basic.txt | tr -d 'occ=' > refine_sites/above_10ns/site_res_no_basic_occ.txt
awk '{print $5}' refine_sites/above_10ns/site_res_basic.txt | tr -d 'occ=' > refine_sites/above_10ns/site_res_basic_occ.txt
awk '{print $5}' refine_sites/above_10ns/site_res_double_basic.txt | tr -d 'occ=' > refine_sites/above_10ns/site_res_double_basic_occ.txt
# analyse raw sites
grep -e ARG -e LYS refine_sites/site_refine.txt | awk -F',' '{print $0" resnum="NF-1}' > refine_sites/all_ns/site_res_basic.txt
grep -e 'ARG.*LYS' -e 'LYS.*ARG' -e 'ARG.*ARG' -e 'LYS.*LYS'  refine_sites/site_refine.txt | awk -F',' '{print $0" resnum="NF-1}' > refine_sites/all_ns/site_res_double_basic.txt
grep -v -e ARG -e LYS refine_sites/site_refine.txt | awk -F',' '{print $0" resnum="NF-1}' > refine_sites/all_ns/site_res_no_basic.txt
awk '{print $5}' refine_sites/all_ns/site_res_no_basic.txt | tr -d 'occ=' > refine_sites/all_ns/site_res_no_basic_occ.txt
awk '{print $5}' refine_sites/all_ns/site_res_basic.txt | tr -d 'occ=' > refine_sites/all_ns/site_res_basic_occ.txt
awk '{print $5}' refine_sites/all_ns/site_res_double_basic.txt | tr -d 'occ=' > refine_sites/all_ns/site_res_double_basic_occ.txt
}

get_site_size (){
rm -f refine_sites/above_10ns/size.*basic.txt
awk '{print $6}' refine_sites/above_10ns/site_res_basic.txt | awk -F , '{print NF}' >> refine_sites/above_10ns/size.basic.txt
awk '{print $6}' refine_sites/above_10ns/site_res_no_basic.txt | awk -F , '{print NF}' >> refine_sites/above_10ns/size.no_basic.txt
}

get_pdb_refine () {
cd refine_sites
mkdir -p indv_sites
while read -r line
do
        #echo $line
        pdb=`echo $line | awk '{print $1}'`
        site=`echo $line | awk '{print $2}'`
        dur=`echo $line | awk '{print $3}'`
        R=`echo $line | awk '{print $4}'`
        occ=`echo $line | awk '{print $5}'`
        res=`echo $line | awk '{print $6}' | tr -d =res\"`
        rm -f indv_sites/$pdb.$site.txt
        for i in $(echo $res | sed 's/,/ /g')
        do
                coords=`grep " $i " ../$pdb/lipid_interactions/Interaction_CARD/coords_Duration.pdb | grep BB | tail -n 1 | awk '{print $4","$(NF-2)}'`
                echo $coords >> indv_sites/$pdb.$site.txt
        done
        if grep -q - indv_sites/$pdb.$site.txt
        then
                sort -t ',' -n -k2,2 indv_sites/$pdb.$site.txt > indv_sites/$pdb.$site.sort.txt
        else
                sort -r -t ',' -n -k2,2 indv_sites/$pdb.$site.txt > indv_sites/$pdb.$site.sort.txt
        fi
        rm -f indv_sites/$pdb.$site.txt
done < sites_above_10ns.txt 
}

get_site_2y () {
cd refine_sites
mkdir -p site_2y
while read -r line
do
	pdb=`echo $line | awk '{print $1}'`
        site=`echo $line | awk '{print $2}'`
        dur=`echo $line | awk '{print $3}'`
        R=`echo $line | awk '{print $4}'`
        occ=`echo $line | awk '{print $5}'`
        res=`echo $line | awk '{print $6}' | tr -d =resid\"`
       	echo $occ > site_2y/$pdb.$site.txt
	for i in $(echo $res | sed 's/,/ /g')
        do
		# pylipid and itp misalign by 1
               	struc=`awk -v var=$((i-1)) '{if ($3 == var) print $0;}' ../../../../../Ecoli_patch/full_complement/chosen/$pdb/Protein.itp | grep BB | head -n 1 | awk '{print $4","$(NF)}'`
		echo $struc >> site_2y/$pdb.$site.txt
        done
done < above_10ns/site_res_basic.txt
}

analyse_2y () {
cd refine_sites
rm -f site_2y/with_helix.txt
rm -f site_2y/just_with_helix.txt
rm -f site_2y/without_helix.txt
rm -f site_2y/just_without_helix.txt
rm -f site_2y/both.txt
for file in `ls site_2y/*txt`
do
	# get all with structure
	if grep -q -e "ARG,H" -e "LYS,H" -e "ARG,1" -e "LYS,1" -e "ARG,2" -e "LYS,2" -e "ARG,3" -e "LYS,3" $file
	then
		head -n 1 $file | tr -d =occ >> site_2y/with_helix.txt
		if ! grep -q -e "ARG,E" -e "LYS,E" -e "ARG,T" -e "LYS,T" -e "ARG,S" -e "LYS,S" -e "ARG,C" -e "LYS,C" $file 
		then
			head -n 1 $file | tr -d =occ >> site_2y/just_with_helix.txt
		fi
	fi
	if grep -q -e "ARG,E" -e "LYS,E" -e "ARG,T" -e "LYS,T" -e "ARG,S" -e "LYS,S" -e "ARG,C" -e "LYS,C" $file
	then
		head -n 1 $file | tr -d =occ >> site_2y/without_helix.txt
		if ! grep -q -e "ARG,H" -e "LYS,H" -e "ARG,1" -e "LYS,1" -e "ARG,2" -e "LYS,2" -e "ARG,3" -e "LYS,3" $file 
		then
			head -n 1 $file | tr -d =occ >> site_2y/just_without_helix.txt
		fi
	fi
	if grep -q -e "ARG,H" -e "LYS,H" -e "ARG,1" -e "LYS,1" -e "ARG,2" -e "LYS,2" -e "ARG,3" -e "LYS,3" -e "ARG,E" -e "LYS,E" -e "ARG,T" -e "LYS,T" -e "ARG,S" -e "LYS,S" -e "ARG,C" -e "LYS,C" $file
	then
		head -n 1 $file | tr -d =occ >> site_2y/both.txt
	fi	
done
}


CD="/sansom/s156a/bioc1535/EC_MEMPROT/5us_analyses/Sites"
cd $CD

#refine_sites
#refine_stats
#get_site_size
#get_pdb_refine
#get_site_2y
analyse_2y
