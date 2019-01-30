#!/bin/bash

# Usage info
show_help() {
cat <<EOF

NAME
	Network construction following Connor et al. (2017) for SLURM jod sheduler
	
SYNOPSIS
	Usage: ${0##*/} [-h] [-a account_name] [-d expected_depth][-o minimum_occurrence_percent] [-r minimum_read_count] INPUT_OTU_MATRIX

DESCRIPTION
	-h	display this help and exit
	
	-a account_name
		Account name for SLURM sbatch -A option.
	
	-b number_of_bootstrap
		Number of bootstraped random noise addition. Default: 1000
	
	-d expected_depth
		Expected sequencing depth to normalize read counts. Default: 0.5 * median read count. Values between 0 and 1 will be use as median read count ratio. Values above 1 will be used as integer read counts.
	
	-e environmental_parameters
		[EXPERIMENTAL OPTION!] Environmental parameter table in Tab separeted format, with parameters as column and samples as rows. The first column have to contain identical sample names as in the OTU table, the first row contains the parameter names (and thus contain one field less). This will include the environmental parameters in the observed matrix Spearman's correlation calculations and include them in the final networks.
	
	-n null_model
		Select the randomization algorithm between cells shuffling over the whole matrix (0), constrained among samples (1), constrained among OTUs (2), or individuals shuffling over the whole matrix (none), with fixed sample read counts (rows), with fixed OTU read counts (columns) or with fixed sample and OTU read counts (both) using the R vegan permatfull function with default parameters. For all models modifying the sample read counts, the read count is re-normalized using the -d parameter. Default: 0
	
	-o minimum_occurrence_percent
		Minimum occurrence percentage threshold to keep an OTU. Default: 0.1 * number of samples. Values between 0 and 1 will be use as the minimum sample number ratio. Values above 1 will be used as the minimum number of samples.
	
	-r minimum_read_count
		Minimum read count threshold to keep a sample. Default: 0.1 * median read count per default. Values between 0 and 1 will be use as median read count ratio. Values above 1 will be used as integer read counts.

AUTHOR
	Guillaume Lentendu

REFERENCE
	Connor, N., Barberán, A., & Clauset, A. (2017). Using null models to infer microbial co-occurrence networks. PLOS ONE, 12(5), e0176751. doi:10.1371/journal.pone.0176751

REPORTING BUGS
	Submit suggestions and bug-reports at <https://github.com/lentendu/NetworkNullHPC/issues>, or send an email to <lentendu@rhrk.uni-kl.de>.
	
COPYRIGHT
	MIT License
	Copyright (c) 2018-2019 Guillaume Lentendu

EOF
}

#set default options
BOOTSTRAP=1000
DEPTH=0.5
ENVMAT="NA"
MINOCC=0.1
MINCOUNT=0.1
NULLM=0

# get options
while getopts ":a:b:d:e:hn:o:r:" opt
do
	case $opt in
		h)	show_help | fmt -s -w $(tput cols)
			exit 1;;
		a)	SLURMACCOUNT=$(echo "#SBATCH -A $OPTARG");;
		b)	BOOTSTRAP=$OPTARG;;
		d)	DEPTH=$OPTARG;;
		e)	ENVMAT=$(readlink -f $OPTARG);;
		n)	NULLM=$OPTARG;;
		o)	MINOCC=$OPTARG;;
		r)	MINCOUNT=$OPTARG;;
		c)	CLEAN=no;;
		\?)	echo "# Error" >&2
			echo "# Invalid option: -$OPTARG" >&2
			show_help | fmt -s -w $(tput cols) >&2
			exit 1;;
		:)	echo "# Error" >&2
			echo "# Option -$OPTARG requires an argument." >&2
			show_help | fmt -s -w $(tput cols) >&2
			exit 1;;
	esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

if [ -z "$1" ]
then 
	echo "# Error: Input OTU matrix is missing." >&2
	show_help | fmt -s -w $(tput cols) >&2
	exit 1
else
	FULLINPUT=$(readlink -f $1) ; shift
	INPUT=${FULLINPUT##*/}
fi

# check me
MYPATH=$(which ${0##*/})
if [ -z "$MYPATH" ]
then
	echo "${0##*/} is not available in any directory listed in the PATH variable."
	echo "Please correct and re-launch."
	echo "Aborting." && exit 1
fi
MYSD=$(dirname $MYPATH})

# check R
RMODULE=$(module -l list 2>&1 | grep -o "^R/[^ ]*")
if [ -z "$(command -v Rscript)" ] || [ -z "$RMODULE" ]
then
	echo "R could not be found. Load a R module first."
	echo "Aborting." && exit 1
fi

# Prepare directories and configuration file
OPTIONS=("$FULLINPUT $ENVMAT $BOOTSTRAP $DEPTH $MINOCC $MINCOUNT $NULLM")
MYCK=$(echo ${OPTIONS[@]} | cat - $FULLINPUT | cksum | awk '{print $1}')
if [ -d "NetworkNullHPC.$MYCK" ]
then
	echo "${0##*/} was already executed on the same input matrix $INPUT with the same options."
	echo "Check outputs with the label $MYCK."
	echo "Aborting"
	exit 1
fi
mkdir NetworkNullHPC.$MYCK && cd NetworkNullHPC.$MYCK
mkdir spearman_noise_r spearman_noise_p spearman_rand_r
cat <(echo "cksum mat env nboot depth minocc mincount nullm") <(echo "$MYCK ${OPTIONS[@]}") | tr " " "\t" > config

# Normalize OTU matrix and get its size
Rscript --vanilla $MYSD/rscripts/clean_mat.R > log.clean_mat.out 2> log.clean_mat.err

# Environmental parameter matrix check
if [ $? -eq 1 ]
then
	echo "Something went wrong in the initial check of the OTU matrix."
	echo "Check that your input matrix have the correct format, as described in the README."
	echo "Also check the content of the initial R script logs below:"
	echo " ### NetworkNullHPC.$MYCK/log.clean_mat.out ###"
	cat log.clean_mat.out
	echo " ######"
	echo " ### NetworkNullHPC.$MYCK/log.clean_mat.err ###"
	cat log.clean_mat.err
	echo " ######"
	echo "Once you identify and correct the issue, delete the NetworkNullHPC.$MYCK directory before re-executing NetworkNullHPC."
	echo " ######"
	echo "Aborting."
	cd ..
	exit 1
elif [ $? -eq 2 ]
then
	echo "The number of samples in the OTU matrix and the environmental parameter tables do not match. Please correct"
	echo "Aborting."
	cd ..
	exit 1
elif [ $? -eq 3 ]
then
	echo "Not all environmental parameters are numeric. Please remove character and/or factoriel parameter(s)."
	echo "Aborting."
	cd ..
	exit 1
elif [ $? -eq 4 ]
then
	echo "The sample names in the OTU matrix and the environmental parameter tables do not match. Please correct"
	echo "Aborting."
	cd ..
	exit 1
fi

# Calculate number of parallel jobs and the amount of memory and time to request
if [ $ENVMAT != "NA" ]
then
	matsize=$(( $(cat nbotu) + $(cat nbenv) ))	
else
	matsize=$(cat nbotu)
fi
pairsize=$((matsize*(matsize-1)/2))
memsize=$(awk -v M=$pairsize 'BEGIN{mem=M/5000000; if(mem!=int(mem)){mem=mem+1};print int(mem)+1}')
blocks=$(( (pairsize/10000+9)/10 ))
if [ $blocks -eq 0 ]; then blocks=1 ; fi
reqtime=$(awk -v M=$pairsize 'BEGIN{T=M*0.0000005+1; if(T!=int(T)){T=T+1};print int(T)}')
if [ $((reqtime*2)) -ge 60 ]
then
	array=1-$BOOTSTRAP
	reqtime2=$reqtime
	PREVSPEAR="Rscript --vanilla $MYSD/rscripts/spearman.R \$SLURM_ARRAY_TASK_ID"
	NULLSPEAR="Rscript --vanilla $MYSD/rscripts/rand_network.R \$SLURM_ARRAY_TASK_ID"
else
	step=$((60/reqtime))
	seqlen=$(seq 1 $step $BOOTSTRAP | sed -n '$=')
	if [ $seqlen -lt 16 ]
	then
		step=$(awk -v B=$BOOTSTRAP 'BEGIN{s=B/16;if(s!=int(s)){print int(s)+1} else print s}')
	fi
	array=1-${BOOTSTRAP}:${step}
	reqtime2=$((reqtime*step))
	MAXSEQ="maxseq=\$(if [ \$((SLURM_ARRAY_TASK_ID + $step)) -gt $BOOTSTRAP ] ; then echo $BOOTSTRAP ; else echo \$((SLURM_ARRAY_TASK_ID+$step-1)) ; fi )"
	PREVSPEAR="for i in \$(seq \$SLURM_ARRAY_TASK_ID \$maxseq); do Rscript --vanilla $MYSD/rscripts/spearman.R \$i ; done"
	NULLSPEAR="for i in \$(seq \$SLURM_ARRAY_TASK_ID \$maxseq); do Rscript --vanilla $MYSD/rscripts/rand_network.R \$i ; done"
fi

# check previous computation(s) and symlink spearman's rho of observed matrix if identical
md5sum $FULLINPUT > md5input
for i in ../NetworkNullHPC.*
do
	if [ "$(basename $i)" != "$(basename $PWD)" ]
	then
		TEST=$(cut -d " " -f 1 $i/md5input | paste - <(echo $FULLINPUT) | md5sum -c --quiet 2> /dev/null)
		if [ -z "$TEST" ]
		then
			PREV=$(readlink -f $i)
			PREVBOOT=$(cut -f 3 $PREV/config | sed -n '2p')
			PREVDEPTH=$(cut -f 4 $PREV/config | sed -n '2p')
			if [ "$PREVBOOT" == "$(ls $PREV/spearman_noise_r/[0-9]*.h5 | wc -l)" ] && [ "$PREVDEPTH" == "$DEPTH" ]
			then
				ln -s $PREV/spearman_noise_r/[0-9]*.h5 $PWD/spearman_noise_r
				ln -s $PREV/spearman_noise_p/[0-9]*.h5 $PWD/spearman_noise_p
				INFOPREV="Spearman's correlations of the observed matrix $INPUT were already calculated in a previous execution (${PREV##*/}), these calculations will be skipped here."
				unset PREVSPEAR
				break
			fi
		fi
	fi
done


# Print infos
cat > info <<EOF

The initial OTU matrix contains $(cat nbsamp_ori) samples and $(cat nbotu_ori) OTUs.
The normalied matrix used for network calculation now contains $(cat nbsamp) samples with a minimum read count of $(cat mincount) and $matsize OTUs with a minimum occurrence of $(cat minocc).
$INFOPREV
EOF
if [ $ENVMAT != "NA" ]
then
	cat >> info <<EOF

The environmental parameter table contains $(cat nbenv) variables.
EOF
fi
cat info

# Estimate the correlation thresholds from one random matrix
cat > sub_range <<EOF
#!/bin/bash

#SBATCH -J range_$MYCK
#SBATCH -o log.%x.out
#SBATCH -e log.%x.err
#SBATCH -t ${reqtime2}
#SBATCH --mem=${memsize}G
$SLURMACCOUNT

module load $RMODULE
Rscript --vanilla $MYSD/rscripts/threshold_range.R

EOF

# Observed matrix Spearman's rho calculation, then find largest connected component in random network at different thresholds
cat > sub_spearman <<EOF
#!/bin/bash

#SBATCH -J spearman_$MYCK
#SBATCH -a $array
#SBATCH -o log.%x.out
#SBATCH -e log.%x.err
#SBATCH --open-mode=append
#SBATCH -t ${reqtime2}
#SBATCH -n 1
#SBATCH --mem=${memsize}G
$SLURMACCOUNT

module load $RMODULE
$MAXSEQ
$PREVSPEAR
$NULLSPEAR

EOF

# find threshold
cat > sub_threshold <<EOF
#!/bin/bash

#SBATCH -J threshold_$MYCK
#SBATCH -o log.%x.out
#SBATCH -e log.%x.err
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
$SLURMACCOUNT

cd spearman_rand_r
parallel -j \$SLURM_CPUS_PER_TASK 'a={} ; eval paste cc_{\$a..\$((a+19))} > paste_cc_\$a' ::: {1..$BOOTSTRAP..20}
parallel -j \$SLURM_CPUS_PER_TASK 'a={} ; eval paste cc_ex_{\$a..\$((a+19))} > paste_cc_ex_\$a' ::: {1..$BOOTSTRAP..20}
wait
paste paste_cc_{1..$BOOTSTRAP..20} > ../rand_cc
paste paste_cc_ex_{1..$BOOTSTRAP..20} > ../rand_cc_ex
rm paste_cc_[0-9]*
rm paste_cc_ex_[0-9]*
cd ..
# at which threshold the largest connected component contain less than 1% of total OTU number in at least 90 % of the random matrices?
POS_TRESH=\$(awk '{print \$1*100}' pos_tresh_range)
paste <(seq \$((POS_TRESH-5)) \$((POS_TRESH+15)) | awk '{print \$1/100}') rand_cc | awk '{count=0;for(i=2;i<=NF;i++){if(\$i<=0.01*$matsize){count+=1}};if(count>=$BOOTSTRAP*0.9){print \$1;exit}}' > threshold
NEG_TRESH=\$(awk '{print \$1*100}' neg_tresh_range)
paste <(seq \$((NEG_TRESH+5)) -1 \$((NEG_TRESH-15)) | awk '{print \$1/100}') rand_cc_ex | awk '{count=0;for(i=2;i<=NF;i++){if(\$i<=0.01*$matsize){count+=1}};if(count>=$BOOTSTRAP*0.9){print \$1;exit}}' > ex_threshold

EOF

# retrieve significant correlations above/below threshold
cat > sub_edges <<EOF
#!/bin/bash

#SBATCH -J edges_$MYCK
#SBATCH -a 1-${blocks}
#SBATCH -o log.%x.out
#SBATCH -e log.%x.err
#SBATCH --open-mode=append
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH --mem=12G
$SLURMACCOUNT

module load $RMODULE
Rscript --vanilla $MYSD/rscripts/edges_r.R \$SLURM_ARRAY_TASK_ID
Rscript --vanilla $MYSD/rscripts/edges_p.R \$SLURM_ARRAY_TASK_ID 
Rscript --vanilla $MYSD/rscripts/edges_p_ex.R \$SLURM_ARRAY_TASK_ID
EOF

# create the final network in form of an edge list
cat > sub_network <<EOF
#!/bin/bash

#SBATCH -J network_$MYCK
#SBATCH -o log.%x.out
#SBATCH -e log.%x.err
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
$SLURMACCOUNT

module load $RMODULE
Rscript --vanilla $MYSD/rscripts/network.R \$SLURM_CPUS_PER_TASK ${blocks}
REXIT=\$?

if [ \$REXIT == 1 ]
then
	echo "Error during network step, exiting"
	exit 1
else
	COOCT=\$(cat threshold)
	COEXT=\$(cat ex_threshold)
fi

if [ \$REXIT == 0 ] || [ \$REXIT == 3 ]
then
	COOCE=\$(sed -n '$=' ../cooccurrence.$MYCK.${INPUT%.*}.txt)
	COOCN=\$(cut -d " " -f 1-2 ../cooccurrence.$MYCK.${INPUT%.*}.txt | tr " " "\n" | sort -u | wc -l)
	COOCSTAT=\$(echo "The co-occurrence network contains \$COOCE edges involving \$COOCN OTUs.")
else
	COOCSTAT=\$(echo "The co-occurrence network is empty.")
fi

if [ \$REXIT == 0 ] || [ \$REXIT == 2 ]
then
	COEXE=\$(sed -n '$=' ../coexclusion.$MYCK.${INPUT%.*}.txt)
	COEXN=\$(cut -d " " -f 1-2 ../coexclusion.$MYCK.${INPUT%.*}.txt | tr " " "\n" | sort -u | wc -l)
	COEXSTAT=\$(echo "The co-exclusion network contains \$COEXE edges involving \$COEXN OTUs.")
else
	COEXSTAT=\$(echo "The co-exclusion network is empty.")
fi

cat > info_start <<EOF2
${0##*/}

## INPUT ##

The input matrix is ${INPUT}.

The input options are:
EOF2

cat > info_end <<EOF3

## OUTPUT ##

The Spearman's rank correlation threshold was set to \$COOCT for the co-occurrence network.
\$COOCSTAT

The Spearman's rank correlation threshold was set to \$COEXT for the co-exclusion network.
\$COEXSTAT

EOF3

cat info_start <(column -t -s \$'\t' config) info info_end > ../info.NetworkNullHPC.$MYCK.${INPUT%.*}.txt

EOF


# Submit to queue
jobid_range=$(sbatch --parsable sub_range)
jobid_spearman=$(sbatch --parsable -d afterok:${jobid_range} sub_spearman)
jobid_threshold=$(sbatch --parsable -d afterok:${jobid_spearman} sub_threshold)
jobid_edges=$(sbatch --parsable -d afterok:${jobid_threshold} sub_edges)
jobid_network=$(sbatch --parsable -d afterok:${jobid_edges} sub_network)

echo ""
echo "The following jobs were submitted to the queue: $jobid_range $jobid_spearman $jobid_threshold $jobid_edges $jobid_network."
echo "Use <squeue> to view the jobs."
echo "The output networks will be available in the files cooccurrence.$MYCK.${INPUT%.*}.txt and coexclusion.$MYCK.${INPUT%.*}.txt once the last job will be completed."
echo ""

cd ..
