#!/bin/bash


# ----- usage ------ #
usage()
{
	echo "CA_to_FullAtom v0.10 [May-30-2019] "
	echo "    Given CA trace in PDB format, restore full atoms. "
	echo ""
	echo "USAGE:  ./CA_to_FullAtom.sh <-i CA_trance> [-o out_file] [-F fix_or_not] "
	echo "                           [-K remove_tmp] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i CA_trance    : Input CA trace in PDB format. "
	echo ""
	echo "***** optional arguments *****"
	echo "-o out_file     : Default output file would be '\${input_name}_FullAtom.pdb'.] "
	echo ""
	echo "-F fix_or_not   : Use FIX or RELAX mode to construct 3D model [default = 1 for FIX] "
	echo ""
	echo "-K remove_tmp   : Remove temporary folder or not. [default = 1 to remove] "
	echo ""
	echo "***** home directory *****"
	echo "-H home         : home directory of From_CA_to_FullAtom. "
	echo "                  [default = `dirname $0`] "
	echo ""
	exit 1
}


#------------------------------------------------------------#
##### ===== get pwd and check BlastSearchHome ====== #########
#------------------------------------------------------------#

#------ current directory ------#
curdir="$(pwd)"

#-------- check usage -------#
if [ $# -lt 1 ];
then
        usage
fi


#---------------------------------------------------------#
##### ===== All arguments are defined here ====== #########
#---------------------------------------------------------#


# ----- get arguments ----- #
#-> required arguments
input=""
#-> others
out_file=""
fix_or_not=1        #-> default: 1 for FIX mode (0 for RELAX mode)
kill_tmp=1          #-> default: kill temporary root
home=`dirname $0`   #-> default: dirname of the 0-th argument

#-> parse arguments
while getopts ":i:o:F:K:H:" opt;
do
	case $opt in
	#-> required arguments
	i)
		input=$OPTARG
		;;
	#-> optional arguments
	o)
		out_file=$OPTARG
		;;
	F)
		fix_or_not=$OPTARG
		;;
	K)
		kill_tmp=$OPTARG
		;;
	H)
		home=$OPTARG
		;;
	#-> others
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done



#---------------------------------------------------------#
##### ===== Part 0: initial argument check ====== #########
#---------------------------------------------------------#

# ------ check home directory ---------- #
if [ ! -d "$home" ]
then
	echo "home directory $home not exist " >&2
	exit 1
fi
home=`readlink -f $home`

# ------ check input fasta ------#
if [ ! -s "$input" ]
then
	echo "input $input not found !!" >&2
	exit 1
fi
input=`readlink -f $input`
fulnam=`basename $input`
relnam=${fulnam%.*}

# ------ check output file ------#
if [ "$out_file" == "" ]
then
	out_file=${relnam}_FullAtom.pdb
fi


#---------------- initialization ----------------#

# --- default directories --#
bin=$home/bin
util=$home/util

# --- create temporary folder --#
DATE=`date '+%Y_%m_%d_%H_%M_%S'`
tmp_root="TMP_FullAtom_${relnam}_${RANDOM}_${DATE}"
mkdir -p $tmp_root



#------------------------------------------------------------#
##### ===== Part 1: Restore Heavy Atom and CB ====== #########
#------------------------------------------------------------#

#-> restore heavy atomand CB
$bin/PDB_Tool -i $input -R 1 -N -1 -o $tmp_root/${relnam}_HeavyAtom.pdb
#-> generate a pseudo model
$util/PDB_Add_Chain $tmp_root/${relnam}_HeavyAtom.pdb A $tmp_root/1pdbA.pdb


#----------------------------------------------------------------#
##### ===== Part 2: Run Modeller to restore FullAtom ====== ######
#----------------------------------------------------------------#


#-> extract SEQRES from original input
$util/PDB_To_SEQ_miss $input $tmp_root/${relnam}.fasta
#-> extract ATOM from HeavyAtom model
$util/PDB_To_SEQ_miss $tmp_root/1pdbA.pdb $tmp_root/1pdbA.fasta
#-> perform DynaProg
$util/Protein_DynaProg $tmp_root/${relnam}.fasta $tmp_root/1pdbA.fasta \
	$tmp_root/${relnam}-1pdbA.fasta
#-> run modeller
$home/build3Dmodel.sh $tmp_root/${relnam}-1pdbA.fasta ${relnam} $tmp_root $home $fix_or_not



#============= post process ============#
mv 1pdbA_${relnam}.B99990001.pdb $out_file
$bin/DeepScore $out_file $input -P 0

#============= kill tmp ================#
if [ $kill_tmp -eq 1 ]
then
	rm -rf $tmp_root
else
	mv $tmp_root TMP_FullAtom_${relnam}
fi

#=============== exit ===============#
exit 0


