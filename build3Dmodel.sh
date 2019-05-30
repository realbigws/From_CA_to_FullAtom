#!/bin/bash

if [ $# -lt 5 ]
then
	echo "Usage: ./build3Dmodel.sh <fasta> <qnam> <template_root> <BdMod_root> <FIXorNOT> "
	exit
fi

#----------- run DeepAlign and build3Dmodel to select best template -----------------#
#function DeepThreader_EPAD()
#{
	#-- input --#
	fasta=${1}           #-> 1st input is the alignment fasta
	qnam=${2}            #-> 2nd input is the qnam
	template_root=${3}   #-> 3rd input is the template root (e.g., pdb_BC100)
	BdMod_root=${4}      #-> 4th input is the build3Dmodel root
	FIXorNOT=${5}        #-> 5th input is to use FIX or RELAX mode

	#--- export mod9v8 root ----#
	if true
	then
		export MODINSTALL9v8=$BdMod_root/modeller9v8
		root=`head -n1 $MODINSTALL9v8/modlib/modeller/config.py | awk -F "'" '{print $2}'`
		if [ "$root" != "$MODINSTALL9v8" ]
		then
			echo "create mod9v8 install_dir in $MODINSTALL9v8/modlib/modeller/config.py"
			echo "install_dir = r'$MODINSTALL9v8'" > $MODINSTALL9v8/modlib/modeller/config.py
			echo "license = 'MODELIRANJE'" >> $MODINSTALL9v8/modlib/modeller/config.py
		fi
	fi
	
	#--- build 3D models ----#
	if [ $FIXorNOT -eq 1 ]
	then
		$BdMod_root/buildModel/build3Dmodel_long_distance -i $fasta -q $qnam -d $template_root -m $BdMod_root/modeller9v8/bin/mod9v8
	else
		$BdMod_root/buildModel/build3Dmodel -i $fasta -q $qnam -d $template_root -m $BdMod_root/modeller9v8/bin/mod9v8
	fi

	#------ exit 0 ------#
	exit 0
#}

