#! /bin/bash

ARGS=`getopt -o i:h --long in:,input:,help,new_gff:,feature:,bias:, -- "$@"`
[ $? -ne 0 ] && usage
eval set -- ${ARGS}

while true; do
	case $1 in
		-i|--in|--input)
			input="$2"
			shift
			;;
		--bias)
			bias="$2"
			shift
			;;
		--feature)
			feature="$2"
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

######################################################
[ -z $feature ] && feature="gene"

awk -v "bias=$bias" -v "feature=$feature" -F "\t" 'BEGIN{OFS="\t"}{if($3==feature){$1=substr($1,match($1,/[0-9]/))+bias; gene2="gene"; print $0}}' $input |sort -n -k1



