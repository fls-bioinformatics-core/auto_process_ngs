#!/bin/bash
#
# Script to process NextSeq data using bcl2fastq2
#
# Use this script to prep an analysis directory, run the
# fastq generation and then make the analysis directories
# in a way that is compatible with auto_process.py.
#
# For example:
#
# process_nextseq.sh setup 151124_NBXXXXX_002_ABCF SampleSheet.csv
# cd 151124_NBXXXXX_002_ABCF_analysis
# process_nextseq.sh make_fastqs
# vi projects.info
# process_nextseq.sh setup_analysis_dirs
#
# Then do:
#
# auto_process.py run_qc etc
#
function usage() {
    # 1: CMD
    case "$1" in
	setup)
	    echo $(basename $0) setup DATADIR SAMPLESHEET
	    ;;
	make_fastqs)
	    echo $(basename $0) make_fastqs
	    ;;
	update_stats)
	    echo $(basename $0) update_stats
	    ;;
	setup_analysis_dirs)
	    echo $(basename $0) setup_analysis_dirs
	    ;;
    esac
    exit
}
function full_path() {
    # Convert relative path to full path by prepending PWD
    if [ -z "$1" ] ; then
	echo $1
    fi
    local is_abs=$(echo $1 | grep "^/")
    if [ -z "$is_abs" ] ; then
	if [ "$1" == "." ] ; then
	    local path=$(pwd)
	else
	    local path=$(pwd)/${1#./}
	fi
    else
	local path=$1
    fi
    echo ${path%/}
}
# Setup a new NextSeq analysis dir
function setup() {
    local datadir=$(full_path $1)
    if [ ! -d "$datadir" ] ; then
	echo No source directory $datadir
	exit 1
    fi
    local samplesheet=$2
    if [ ! -f "$samplesheet" ] ; then
	echo No sample sheet $samplesheet
	exit 1
    fi
    # Create analysis dir
    local analysis_dir=$(basename $datadir)_analysis
    if [ -d "$analysis_dir" ] ; then
	echo Analysis dir $analysis alreay exists
	exit 1
    fi
    echo Creating $analysis_dir
    mkdir $analysis_dir
    cd $analysis_dir
    mkdir logs
    mkdir ScriptCode
    # Set up sample sheet
    echo Fetching $samplesheet
    cp $samplesheet SampleSheet.orig.csv
    samplesheet=$(pwd)/custom_SampleSheet.csv
    cp SampleSheet.orig.csv $samplesheet
    dos2unix $samplesheet ; mac2unix $samplesheet
    # Make auto_process.info
    cat > auto_process.info <<EOF
analysis_dir	$analysis_dir
bases_mask	.
data_dir	$datadir
primary_data_dir	.
project_metadata	projects.info
sample_sheet	$samplesheet
stats_file	statistics.info
unaligned_dir	bcl2fastq
EOF
    # Make metadata.info
    cat > metadata.info <<EOF
assay	.
bcl2fastq_software	.
platform	nextseq
run_number	.
source	.
EOF
}
# Generate Fastqs
function make_fastqs() {
    # 1: bcl2fastq dir
    local bcl2fastq=$(which bcl2fastq)
    if [ -z "$bcl2fastq" ] ; then
	echo No bcl2fastq
	exit 1
    fi
    local bcl2fastq_version=$(bcl2fastq --version 2>&1 | grep ^bcl2fastq | cut -d" " -f2 | cut -c2-)
    mkdir $1
    # Update metadata file
    grep -v ^bcl2fastq_software metadata.info > metadata.info.tmp
    cat >> metadata.info.tmp <<EOF
bcl2fastq_software	('$bcl2fastq', '$(basename $bcl2fastq)', '$bcl2fastq_version')
EOF
    mv metadata.info.tmp metadata.info
    # Get the data and sample sheet locations
    local datadir=$(grep ^data_dir auto_process.info | cut -f 2)
    local samplesheet=$(grep ^sample_sheet auto_process.info | cut -f 2)
    local cmd="bcl2fastq -R $datadir -o $1 --sample-sheet $samplesheet"
    echo Running $cmd
    $cmd 1>logs/bcl2fastq.log 2>&1
    if [ $? -ne 0 ] ; then
	echo Non-zero exit status from bcl2fastq
	exit 1
    fi
}
# Get projects from the bcl2fastq directory
function get_projects() {
    # 1: bcl2fastq dir
    local projects=
    for d in $(ls $1) ; do
	if [ -d $1/$d ] ; then
	    if [ $d != "Reports" ] && [ $d != "Stats" ] ; then
		projects="$projects $d"
	    fi
	fi
    done
    echo $projects
}
# Get undetermined fastqs from the bcl2fastq directory
function get_undetermined() {
    # 1: bcl2fastq dir
    local undetermined=
    for f in $(ls bcl2fastq) ; do
	if [ -f bcl2fastq/$f ] && [ ! -z "$(echo $f | grep ^Undetermined_S)" ] ; then
            undetermined="$undetermined $f"
	fi
    done
    echo $undetermined
}
# Derive sample name from fastq file
function get_sample_name() {
    # 1: fastq file
    local fq=$1
    local sample_name=
    for field in $(echo $(basename ${fq%%.*}) | tr _ " ") ; do
	 if [ -z "$(echo $field | grep ^S[0-9])" ] ; then
            if [ -z "$sample_name" ] ; then
		sample_name=$field
            else
		sample_name="${sample_name}_${field}"
            fi
         else
            break
         fi
    done
    echo $sample_name
}
# Count reads in fastq file
function count_reads() {
    # 1: fastq file
    if [ -z "$(echo $1 | grep '.gz$')" ] ; then
	echo $(($(wc -l $1 | cut -d" " -f1)/4))
    else
	echo $(($(zcat $1 | wc -l)/4))
    fi
}
# Determine if run is paired end
function is_paired_end() {
    # 1: bcl2fastq dir
    if [ -z "$(echo $(get_undetermined $1) | grep _R2_)" ] ; then
	echo N
    else
	echo Y
    fi
}
# Generate the statistics.info file
function generate_stats() {
    # 1: bcl2fastq dir
    echo "Generating statistics.info file"
    paired_end=$(is_paired_end $1)
    cat > statistics.info <<EOF
#Project	Sample	Fastq	Size	Nreads	Paired_end
EOF
    for project in $(get_projects $1) ; do
	for fq in $(ls $1/$project) ; do
	    sample_name=$(get_sample_name $fq)
	    size=$(du -sh $1/$project/$fq | cut -f1)
	    nreads=$(count_reads $1/$project/$fq)
	    cat >> statistics.info <<EOF
${project}	${sample_name}	$(basename $fq)	$size	$nreads	$paired_end
EOF
	done
    done
    for fq in $(get_undetermined $1) ; do
	sample_name=$(get_sample_name $fq)
	size=$(du -sh $1/$fq | cut -f1)
	nreads=$(count_reads $1/$fq)
	cat >> statistics.info <<EOF
Undetermined_indices	${sample_name}	$(basename $fq)	$size	$nreads	$paired_end
EOF
    done
}
# Generate the projects.info file
function make_projects_info() {
    # 1: bcl2fastq dir
    echo "Generating projects.info file"
    cat > projects.info <<EOF
#Project	Samples	User	Library	Organism	PI	Comments
EOF
    for project in $(get_projects $1) ; do
	echo $project
	local samples=
	for fq in $(ls $1/$project) ; do
	    sample_name=$(get_sample_name $fq)
	    if [ -z "$samples" ] ; then
		samples="${sample_name}"
	    else
		samples="${samples},${sample_name}"
	    fi
	done
	cat >> projects.info <<EOF
${project}	${samples}	.	.	.	.	.
EOF
    done
}
# Set up the analysis dirs
function setup_analysis_dirs() {
    # Populate projects
    while read line ; do
	if [ -z "$(echo $line | grep ^#)" ] ; then
	    # Gather metadata
	    project=$(echo $line | cut -d" " -f1)
	    samples=$(echo $line | cut -d" " -f2)
	    user=$(echo $line | cut -d" " -f3)
	    library=$(echo $line | cut -d" " -f4)
	    organism=$(echo $line | cut -d" " -f5)
	    pi=$(echo $line | cut -d" " -f6)
	    comments=$(echo $line | cut -d" " -f7-)
	    nsamples=$(echo $samples | tr -cd , | wc -m)
	    echo $project
	    # Create the project
	    mkdir $project
	    mkdir $project/fastqs
	    mkdir $project/ScriptCode
	    ln $1/$project/*.fastq.gz $project/fastqs/
	    cat >> $project/README.info <<EOF
Run	.
Platform	nextseq
User	$user
PI	$pi
Organism	$organism
Library type	$library
Paired_end	$(is_paired_end $1)
Samples	$nsamples sample (${samples})
Comments	$comments
EOF
	fi
    done < projects.info
    # Undetermined
    mkdir undetermined
    mkdir undetermined/fastqs
    for fq in $(get_undetermined $1) ; do
	ln $1/$fq undetermined/fastqs/
	    cat >> $project/README.info <<EOF
Run	.
Platform	nextseq
User	.
PI	.
Organism	.
Library type	.
Paired_end	$(is_paired_end $1)
Comments	.
EOF
    done
}
# Main script
# Sort out command line
if [ $# -lt 1 ] ; then
   echo Usage: $(basename $0) CMD ARGS
   echo 
   echo CMD can be setup, make_fastqs or setup_analysis_dirs
   exit 1
fi
cmd=$1
#echo $cmd
#echo $#
case "$cmd" in
    setup)
	if [ $# -ne 3 ] ; then
	    usage $cmd
	fi
	setup "$2" "$3"
	;;
    make_fastqs)
	if [ $# -ne 1 ] ; then
	    usage $cmd
	fi
	make_fastqs bcl2fastq
	generate_stats bcl2fastq
	make_projects_info bcl2fastq
	;;
    update_stats)
	if [ $# -ne 1 ] ; then
	    usage $cmd
	fi
	generate_stats bcl2fastq
	;;
    setup_analysis_dirs)
	if [ $# -ne 1 ] ; then
	    usage $cmd
	fi
	setup_analysis_dirs bcl2fastq
	;;
    *)
	echo Unrecognised command: $cmd
	exit 1
	;;
esac
#generate_stats bcl2fastq
#make_projects_info bcl2fastq
##setup_analysis_dirs bcl2fastq
