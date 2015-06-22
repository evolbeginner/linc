#! /bin/bash

########################################################################
function usage(){
	echo -e "USAGE of `basename $0`"
	cat <<EOF
--in_blast	e.g., for lincRNA blast results in the tabular blast format
--in_gff	e.g., for gff file of the genome (ATH) in the gff3 format
--outdir
--clean|force	remove outdir if outdir exists already?
--feature
--bias
-h|--help	usage
EOF
	echo
	exit
}

########################################################################
[ $# -eq 0 ] && usage

ARGS=`getopt -o h --long in_blast:,in_gff:,ori_linc_gff:,outdir:,clean,force,feature:,bias:,help, -- "$@"`
[ $? -ne 0 ] && usage
eval set -- ${ARGS}

while true; do
	case $1 in
		--in_blast)
			in_blast="$2"
			shift
			;;
		--in_gff)
			in_gff="$2"
			shift
			;;
		--ori_linc_gff)
			ori_linc_gff="$2"
			shift
			;;
		--outdir)
			outdir="$2"
			shift
			;;
		--force|--clean)
			force=true
			;;
		--bias)
			bias="$2"
			shift
			;;
		--feature)
			feature="$2"
			shift
			;;
		-h|--help)
			usage
			break
			;;
		--)
			shift
			break
			;;
		*)
			shift
			echo -e "\nparameter error!\n"
			usage
			exit
	esac	
	shift
done

[ -z $feature ] && feature="gene"
[ -z $outdir ] && outdir="."
[ -z $ori_linc_gff ] && ori_linc_gff=~/project/linc/data/gff/liu_linc.gff


if [ -e $outdir ]; then
	if [ ! -z $force ]; then
		rm -rf $outdir
	else
		echo "outdir $outdir already exists! Exiting ......"
		exit
	fi
fi

if [ $outdir != "." ]; then
	mkdir -p $outdir
fi


dirname=`dirname $0`
generate_pseudo_blastoutput_gff="$dirname/generate_pseudo_blastoutput_gff.sh"
get_gene_from_gff="$dirname/get_gene_from_gff.sh"

speciesB_gff="$outdir/speciesB.biased_chr.gff"
new_linc_gff="$outdir/blast_result.pseudo.gff"
final_linc_gff="$outdir/final_linc.gff"

########################################################################
bash $generate_pseudo_blastoutput_gff -i $in_blast --bias $bias --outdir $outdir

bash $get_gene_from_gff --feature $feature --bias $bias -i $in_gff > $speciesB_gff

cat $ori_linc_gff $new_linc_gff > $final_linc_gff


