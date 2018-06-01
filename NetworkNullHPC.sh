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
	
	-d expected_depth
		Expected sequencing depth to normalize read counts. Default: 0.5 * median read count. Values between 0 and 1 will be use as median read count ratio. Values above 1 will be used as integer read counts.
	
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
	Copyright (c) 2018 Guillaume Lentendu

EOF
}

# get options
while getopts ":a:d:ho:r:" opt
do
	case $opt in
		h)	show_help | fmt -s -w $(tput cols)
			exit 1;;
		a)	SLURMACCOUNT=$(echo "#SBATCH -A $OPTARG");;
		d)	DEPTH=$OPTARG;;
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
i=$(cat <(date +%s) $FULLINPUT | cksum | awk '{print $1}')
mkdir NetworkNull.$i && cd NetworkNull.$i
mkdir spearman_noise_r spearman_noise_p spearman_rand_r
TAB=$'\t'
cat > config <<EOF
cksum${TAB}mat${TAB}depth${TAB}minocc${TAB}mincount
$i${TAB}$FULLINPUT${TAB}${DEPTH:-0.5}${TAB}${MINOCC:-0.1}${TAB}${MINCOUNT:-0.1}
EOF

# Normalize OTU matrix and get its size
Rscript --vanilla $MYSD/rscripts/clean_mat.R > log.clean_mat.out 2> log.clean_mat.err

# Calculate number of parallel jobs and the amount of memory to request
matsize=`cat nbotu`
pairsize=$((matsize*(matsize-1)/2))
memsize=$((pairsize/5000000+1))
blocks=$(( (pairsize/10000+9)/10 ))
if [ $(echo ${MINOCC:-0.1}"<=1" | bc) -eq 1 ]; then MINOCCINFO=" * number of samples" ; fi
if [ $(echo ${MINOCC:-0.1}"<=1" | bc) -eq 1 ]; then MINCOUNTINFO=" * the median read count" ; fi
cat > info <<EOF

The initial OTU matrix contains $(cat nbsamp_ori) samples and $(cat nbotu_ori) OTUs.
The normalied matrix used for network calculation now contains $(cat nbsamp) samples with a minimum read count of ${MINCOUNT:-0.1}$MINCOUNTINFO and $matsize OTUs with a minimum occurrence of ${MINOCC:-0.1}$MINOCCINFO.
EOF
cat info

# Estimate the correlation thresholds from one random matrix
cat > sub_range <<EOF
#!/bin/bash

#SBATCH -J estimate_range_$i
#SBATCH -o log.estimate_range.%j.out
#SBATCH -e log.estimate_range.%j.err
#SBATCH -t 01:00:00
#SBATCH --mem=${memsize}G
$SLURMACCOUNT

module load $RMODULE
Rscript --vanilla $MYSD/rscripts/threshold_range.R

EOF

# Observed matrix Spearman's rho calculation, then find largest connected component in random network at different thresholds
cat > sub_spearman <<EOF
#!/bin/bash

#SBATCH -J spearman_$i
#SBATCH -a 1-1000%400
#SBATCH -o log.%x.%A.out
#SBATCH -e log.%x.%A.err
#SBATCH --open-mode=append
#SBATCH -t 03:00:00
#SBATCH -N 1
#SBATCH --mem=${memsize}G
$SLURMACCOUNT

module load $RMODULE
Rscript --vanilla $MYSD/rscripts/spearman.R \$SLURM_ARRAY_TASK_ID
Rscript --vanilla $MYSD/rscripts/rand_network.R \$SLURM_ARRAY_TASK_ID

EOF

# find threshold
cat > sub_threshold <<EOF
#!/bin/bash

#SBATCH -J threshold_$i
#SBATCH -o log.%x.%j.out
#SBATCH -e log.%x.%j.err
#SBATCH --open-mode=append
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
$SLURMACCOUNT

cd spearman_rand_r
parallel -j \$SLURM_CPUS_PER_TASK 'a={} ; eval paste cc_{\$a..\$((a+19))} > paste_cc_\$a' ::: {1..1000..20}
parallel -j \$SLURM_CPUS_PER_TASK 'a={} ; eval paste cc_ex_{\$a..\$((a+19))} > paste_cc_ex_\$a' ::: {1..1000..20}
wait
paste paste_cc_{1..1000..20} > ../rand_cc
paste paste_cc_ex_{1..1000..20} > ../rand_cc_ex
rm paste_cc_[0-9]*
rm paste_cc_ex_[0-9]*
cd ..
# at which threshold the largest connected component contain less than 1% of total OTU number in 95 % of the random matrices?
POS_TRESH=\$(awk '{print \$1*100}' pos_tresh_range)
paste <(seq \$((POS_TRESH-10)) \$((POS_TRESH+10)) | awk '{print \$1/100}') rand_cc | awk -v S=\$(cat $i_size) '{count=0;for(i=2;i<=NF;i++){if(\$i<=0.01*S){count+=1}};if(count>=1000*0.95){print \$1;exit}}' > threshold
NEG_TRESH=\$(awk '{print \$1*100}' neg_tresh_range)
paste <(seq \$((NEG_TRESH-10)) \$((NEG_TRESH+10)) | awk '{print -\$1/100}') rand_cc_ex | awk -v S=\$(cat $i_size) '{count=0;for(i=2;i<=NF;i++){if(\$i<=0.01*S){count+=1}};if(count>=1000*0.95){print \$1;exit}}' > ex_threshold

EOF

# retrieve significant correlations above/below threshold
cat > sub_edges <<EOF
#!/bin/bash

#SBATCH -J edges_$i
#SBATCH -a 1-${blocks}
#SBATCH -o log.%x.%A.out
#SBATCH -e log.%x.%A.err
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

#SBATCH -J network_$i
#SBATCH -o log.%x.%j.out
#SBATCH -e log.%x.%j.err
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
$SLURMACCOUNT

module load $RMODULE
Rscript --vanilla $MYSD/rscripts/network.R $i \$SLURM_CPUS_PER_TASK ${blocks}
	cat > info_start <<EOF2
		${0##*/}
		
		## INPUT ##
		
		The input matrix is ${INPUT}.
		
		The input options are:
		column -t -s \$'\t' config
EOF2
	cat > info_end <<EOF3
		
		## OUTPUT ##
		
		The Spearman's rank correlation threshold was set to \$(cat threshold) for the co-occurrence network.
		The co-occurrence network contains \${wc -l ../cooccurrence.$i.${INPUT%.*}.txt} edges involving \${cut -d " " -f 1-2 ../cooccurrence.$i.${INPUT%.*}.txt | tr " " "\n" | sort -u | wc -l} OTUs.
		
		The Spearman's rank correlation threshold was set to \$(cat ex_threshold) for the co-exclusion network.
		The co-exclusion network contains \${wc -l ../coexclusion.$i.${INPUT%.*}.txt} edges involving \${cut -d " " -f 1-2 ../coexclusion.$i.${INPUT%.*}.txt | tr " " "\n" | sort -u | wc -l} OTUs.
EOF3
cat info_start info info_end > ../info.NetworkNullHPC.$i.${INPUT%.*}.txt
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
echo "The output networks will be available in the files cooccurrence.$i.${INPUT%.*}.txt and coexclusion.$i.${INPUT%.*}.txt once the last job will be completed."
echo ""

cd ..
