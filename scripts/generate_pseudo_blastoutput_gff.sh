#! /bin/bash

source ~/program/shell/colors_list.sh 

function usage(){
	basename=`basename $0`
	echo
	echo "The usage of $basename"
	echo -e "${Red}bash $basename <-i|--in|--input>${NC}"
	echo -e "${Blue}Options:${NC}
--outdir:	default ./
--new_gff:	(gff0 format)
--new_blast_output
--bias:		default 0
"
	echo
	exit 
}


###################################################################################
ARGS=`getopt -o i:h --long in:,input:,help,outdir:,new_gff:,new_blast:,new_blast_output:,bias:, -- "$@"`
[ $? -ne 0 ] && usage
eval set -- ${ARGS}

[ $# -eq 0 ] && usage

###################################################################################
bias=0
outdir="."

while true; do
	case $1 in
		-i|--in|--input)
			original_blast_output="$2"
			shift
			;;
		--new_gff)
			new_gff="$2"
			shift
			;;
		--new_blast|--new_blast_output)
			new_blast_output="$2"
			shift
			;;
		--bias)
			bias="$2"
			shift
			;;
		--outdir)
			outdir="$2"
			shift
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


[ -z $original_blast_output ] && echo "original_blast_output has to be given. Exiting ......" && exit

original_basename=`basename $original_blast_output`

[ -z $new_gff ] && new_gff=$outdir/$original_basename.pseudo.gff
[ -z $new_blast_output ] && new_blast_output=$outdir/$original_basename.pseudo.blast


###################################################################################
#awk -F"\t" 'BEGIN{OFS="\t"}{new2=$2":"$9"-"$10":"$1; print $2, new2, $9, $10 }' $original_blast_output > $new_gff
#awk -F"\t" 'BEGIN{OFS="\t"}{$2=$2":"$9"-"$10":"$1; print $0}' $original_blast_output > $new_blast_output

awk -F"\t" -v "bias=$bias" 'BEGIN{OFS="\t"}{$2=substr($2,match($2,/[0-9]+/))+bias; new2=$2":"$9"-"$10":"$1; print $2,new2,$9,$10}' $original_blast_output > $new_gff

awk -F"\t" -v "bias=$bias" 'BEGIN{OFS="\t"}{$2=substr($2,match($2,/[0-9]+/))+bias;   $2=$2":"$9"-"$10":"$1; print $0}' $original_blast_output > $new_blast_output


