#!/bin/bash

convert () {
mkdir -p Convert/$1
cd Convert/$1
python3 $GIT/cg2at/cg2at.py -c ../../ -a cryo_em2.pdb -vvv -clean -w tip3p -ff charmm36-mar2019 -nt -ct -fg martini_2-2_charmm36
}

for i in `ls *gro`
do
	convert $i
done
