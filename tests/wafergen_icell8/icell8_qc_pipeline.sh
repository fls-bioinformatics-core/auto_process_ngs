#!/bin/bash -e
#
# ICell8 QC pipeline implemented as a shell script for testing

# Command line options
FASTQ_NAMES="icell8"
FASTQ_SCREEN_MAMMALIAN_CONF=
FASTQ_SCREEN_CONTAMINANTS_CONF=
OUTDIR="test.out"

function usage() {
    echo "Usage: $(basename $0) options"
    echo "Options"
    echo "  -h"
    echo "  -n FASTQ_NAMES"
    echo "  -m FASTQ_SCREEN_MAMMALIAN_CONF"
    echo "  -c FASTQ_SCREEN_CONTAMINANTS_CONF"
    echo "  -o OUTDIR"
}

while getopts ":n:m:c:o:h" opt ; do
    case "$opt" in
	n)
	    FASTQ_NAMES="${OPTARG}"
	    ;;
	m)
	    FASTQ_SCREEN_MAMMALIAN_CONF="${OPTARG}"
	    ;;
	c)
	    FASTQ_SCREEN_CONTAMINANTS_CONF="${OPTARG}"
	    ;;
	o)
	    OUTDIR="{OPTARG}"
	    ;;
	h)
	    usage
	    exit 0
	    ;;
	*)
	    usage
	    echo "Unrecognised option supplied" >&2
	    exit 1
	    ;;
    esac
done

# Check inputs
for conf in "$FASTQ_SCREEN_MAMMALIAN_CONF" "$FASTQ_SCREEN_CONTAMINANTS_CONF" ; do
    if [ -z "$conf" ] ; then
	echo Missing Fastq Screen conf file >&2
	exit 1
    elif [ ! -f $conf ] ; then
	echo $conf: not found >&2
	exit 1
    fi
done

# Check programs
for rq in cutadapt fastq_screen ; do
    if [ -z "$(which $rq)" ] ; then
	echo Missing $rq
	exit 1
    fi
done

# Set up Fastq Screen "nohits" strings
# Produces strings of the form "FQST:0000"
# The number of zeroes for each should match the
# number of databases declared in the conf file
function get_no_hits() {
    local nohits=
    while read line ; do
	if [ ! -z "$(echo $line | grep -w ^DATABASE)" ]  ; then
	    nohits="${nohits}0"
	fi
    done <$1
    echo "FQST:${nohits}"
}
MAMMALIAN_NOHITS=$(get_no_hits $FASTQ_SCREEN_MAMMALIAN_CONF)
CONTAMINANTS_NOHITS=$(get_no_hits $FASTQ_SCREEN_CONTAMINANTS_CONF)
echo Mammalian nohits: $MAMMALIAN_NOHITS
echo Contaminant nohits: $CONTAMINANTS_NOHITS

# Count barcodes and distinct UMIs
function fastq_stats() {
    local barcodes=$(cat $@ | sed -n '2~4p' | cut -c1-11 | sort -u)
    local nreads=$(($(cat $@ | wc -l)/4))
    #echo $@
    #echo "#reads $nreads"
    #echo "Barcodes:"
    for barcode in $barcodes ; do
	local nreads=$(cat $@ | sed -n '2~4p' | grep "^$barcode" | wc -l)
	local umis=$(cat $@ | sed -n '2~4p' | grep "^$barcode" | cut -c12-21 | sort -u | wc -l)
	echo -e "${barcode}\t${nreads}\t${umis}"
    done
}

# Return list of read numbers with barcode quality < 10
function find_bad_barcodes() {
    sed -n '4~4p' $1 | cut -c1-11 | grep -n '[!"#$%&*''()+]' | cut -d":" -f1
}

# Return list of read numbers with UMI quality < 30
function find_bad_umis() {
    sed -n '4~4p' $R1 | cut -c12-21 | grep -n '[^@ABCDEFGHIJ]' | cut -d":" -f1
}

# Remove reads
# Must be supplied with FASTQ file and a list of read
# numbers to be removed
function remove_reads() {
    local i
    local sed_expr
    local start
    local end
    if [ ! -z "$2" ] ; then
	for i in $2 ; do
	    end=$(($i*4))
	    start=$(($end-3))
	    sed_expr="${sed_expr}${start},${end}d;"
	done
	sed -e "${sed_expr}" $1
    else
	cat $1
    fi
}

# Remove read pairs with low barcode and/or UMI quality scores
function quality_filter_fastq_pair() {
    echo Quality filter: $1
    local bad_barcodes=$(find_bad_barcodes $1)
    remove_reads "$1" "$bad_barcodes" > ${3}.tmp
    remove_reads "$2" "$bad_barcodes" > ${4}.tmp
    local bad_umis=$(find_bad_umis $3)
    remove_reads "${3}.tmp" "$bad_umis" > $3
    remove_reads "${4}.tmp" "$bad_umis" > $4
    rm -f ${3}.tmp ${4}.tmp
}

# Trim reads using cutadapt
function trim_fastq_pair() {
    echo Trim reads: $1
    cutadapt \
	-a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
        -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
        -a AAAAAAAA \
        -a TTTTTTTT \
        -m 20 \
        --trim-n \
        --max-n 0.7 \
        -q 25 \
	-o $4 -p $3 \
	$2 $1
}

# Filter reads failing contamination checks
function filter_contaminants() {
    echo Filter contaminants: $1
    local fastq_base=${2%.*}
    fastq_screen \
	--subset 0 \
	--conf $FASTQ_SCREEN_CONTAMINANTS_CONF \
	--tag \
	--force \
	$2
    mv ${fastq_base}.tagged.fastq ${fastq_base}.contaminants.tagged.fastq
    
    local contaminants_hits=$(sed -n '1~4p' ${fastq_base}.contaminants.tagged.fastq | grep -n -v ${CONTAMINANTS_NOHITS}'$' | cut -d":" -f1)
    ##echo "CONTAMINANT HITS: $contaminant_hits"
    fastq_screen \
	--subset 0 \
	--conf $FASTQ_SCREEN_MAMMALIAN_CONF \
	--tag \
	--force \
	$2
    mv ${fastq_base}.tagged.fastq ${fastq_base}.mammalian.tagged.fastq
    local mammalian_hits=$(sed -n '1~4p' ${fastq_base}.mammalian.tagged.fastq | grep -n -v ${MAMMALIAN_NOHITS}'$' | cut -d":" -f1)
    ##echo "MAMMALIAN HITS: $mammalian_hits"
    local contaminated_reads
    for ii in $contaminants_hits ; do
	local remove_read=yes
	for jj in $mammalian_hits ; do
	    if [ $ii -eq $jj ] ; then
		remove_read=
		break
	    fi
	done
	if [ ! -z "$remove_read" ] ; then
	    contaminated_reads="$contaminated_reads $ii"
	fi
    done
    remove_reads "$1" "$contaminated_reads" > $3
    remove_reads "$2" "$contaminated_reads" > $4
    for ext in screen.html screen.png screen.txt ; do
	rm -f ${fastq_base}_${ext}
    done
}
# Main script 

# Temporary working dir
wd=$(mktemp -d --suffix=".verify" --tmpdir=$(pwd))
for fq_base in $FASTQ_NAMES ; do
    R1=${fq_base}.r1.fastq
    R2=${fq_base}.r2.fastq
    zcat $(dirname $0)/test-data/$R1.gz >$wd/$R1
    zcat $(dirname $0)/test-data/$R2.gz >$wd/$R2
done
cd $wd

# Process supplied FASTQs
FASTQS_IN=
FILTERED_FASTQS=
TRIMMED_FASTQS=
CONTAMINANT_FILTERED_FASTQS=
for fq_base in $FASTQ_NAMES ; do
    # Initial fastqs
    R1=${fq_base}.r1.fastq
    R2=${fq_base}.r2.fastq
    FASTQS_IN="$FASTQS_IN $R1"
    # Quality filter
    FILTERED_R1=${fq_base}.filtered.r1.fastq
    FILTERED_R2=${fq_base}.filtered.r2.fastq
    FILTERED_FASTQS="$FILTERED_FASTQS $FILTERED_R1"
    quality_filter_fastq_pair $R1 $R2 $FILTERED_R1 $FILTERED_R2
    # Read trimming
    TRIMMED_R1=${fq_base}.trimmed.r1.fastq
    TRIMMED_R2=${fq_base}.trimmed.r2.fastq
    TRIMMED_FASTQS="$TRIMMED_FASTQS $TRIMMED_R1"
    trim_fastq_pair $FILTERED_R1 $FILTERED_R2 $TRIMMED_R1 $TRIMMED_R2
    # Contamination filter
    CONTAMINANT_FILTERED_R1=${fq_base}.contaminant_filtered.r1.fastq
    CONTAMINANT_FILTERED_R2=${fq_base}.contaminant_filtered.r2.fastq
    CONTAMINANT_FILTERED_FASTQS="$CONTAMINANT_FILTERED_FASTQS $CONTAMINANT_FILTERED_R1"
    filter_contaminants $TRIMMED_R1 $TRIMMED_R2 \
	$CONTAMINANT_FILTERED_R1 $CONTAMINANT_FILTERED_R2
done

# Collect counts for each stage
echo "*** Initial statistics ***"
fastq_stats $FASTQS_IN >stats.initial
echo "*** Post quality filter statistics ***"
fastq_stats $FILTERED_FASTQS >stats.filtered
echo "*** Post trimming statistics ***"
fastq_stats $TRIMMED_FASTQS >stats.trimmed
echo "*** Post contaminant filter statistics ***"
fastq_stats $CONTAMINANT_FILTERED_FASTQS >stats.contaminants

# Combine stats
echo -e "#Barcodes\tNreads_initial\tUMIs_initial\tNreads_filtered\tUMIs_filtered\tNreads_trimmed\tUMIs_trimmed\tNreads_final\tUMIs_final" >stats
paste \
    stats.initial \
    stats.filtered \
    stats.trimmed \
    stats.contaminants | \
    cut -f1-3,5-6,8-9,11-12 >>stats
cd ..

# Copy to final location
if [ -d "$OUTDIR" ] ; then
    echo Removing existing directory $OUTDIR
    rm -f $OUTDIR/*
    rmdir $OUTDIR
fi
echo Copy outputs to $OUTDIR
mkdir $OUTDIR
mv $wd/* $OUTDIR
rmdir $wd
##
#
