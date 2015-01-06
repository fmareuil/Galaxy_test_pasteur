#! /bin/sh
#$ -S /bin/sh
#$ -cwd -sync y -now n
#$ -o $JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e $JOB_NAME.$JOB_ID.$TASK_ID.err

err=0
dir=/local/gensoft/libexec/taxoptimizer
tmpnotaxo="tmpnotaxo.$$.${SGE_TASK_ID}"
outputnotaxo=""
# arguments: <outfile> <inlist> <taxOopts...>

output=$1
shift
input=$(cat $1 | head -n ${SGE_TASK_ID} | tail -n 1)
shift

# Run the blast command on the expected sequence chunk
if test "$1" = '-f'
   then 
	shift
	outputnotaxo=$1 
	shift
	#echo "1- ${dir}/seqpos input | taxoptimizer -i stdin -f ${tmpnotaxo} ${1+"$@"} || err=1"
	${dir}/seqpos ${input} | taxoptimizer -i stdin -f ${tmpnotaxo} ${1+"$@"} || err=2
   else
	#echo "2- ${dir}/seqpos ${input} | taxoptimizer -i stdin ${1+"$@"} || err=1"
        ${dir}/seqpos ${input} | taxoptimizer -i stdin ${1+"$@"} || err=3
fi

if test "${outputnotaxo}" != "" ; then
  ${dir}/resmerge ${outputnotaxo} ${tmpnotaxo} || err=6
fi

# Merge results if exists ...
if test -s ${SGE_STDOUT_PATH}; then
  ${dir}/resmerge ${output} ${SGE_STDOUT_PATH} || err=4
fi
if test -s ${SGE_STDERR_PATH}; then
  ${dir}/resmerge ${output}.err ${SGE_STDERR_PATH} || err=5
fi


# Cleanup
rm -f ${SGE_STDOUT_PATH} ${SGE_STDERR_PATH}
rm -f ${tmpnotaxo}

exit $err
