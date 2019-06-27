#!/bin/bash
#
# Count barcodes and distinct UMIs in arbitrary ICell8 R1 FASTQs

function usage() {
    echo "Usage: $(basename $0) WELL_LIST FASTQ_R1 [FASTQ_R1...]"
}

# Command line
well_list=$1
fastqs=${@:2}
if [ -z "$well_list" ] || [ -z "$fastqs" ] ; then
    usage
    exit 1
fi

# Get expected barcodes
expected_barcodes=$(cut -f6 $well_list | grep -v -w "Barcode")

# Version of 'cat' to use
for fastq in $fastqs ; do
    ext=${fastq##*.}
    if [ "$ext" == "gz" ] ; then
	CAT=zcat
    else
	CAT=cat
    fi
    break
done

# Extract barcode and UMI data
barcodes_and_umis_file=$(mktemp --tmpdir=$(pwd) --suffix=.barcodes_and_umi)
barcodes_file=$(mktemp --tmpdir=$(pwd) --suffix=.barcodes)
$CAT $fastqs | sed -n '2~4p' | cut -c1-21 | sort > $barcodes_and_umis_file
cat $barcodes_and_umis_file | cut -c1-11 | sort > $barcodes_file

# Report counts
barcodes=$(cut -f6 $well_list | grep -v -w "Barcode")
for barcode in $barcodes ; do
    nreads=$(grep -c "^$barcode" $barcodes_file)
    umis=$(grep "^$barcode" $barcodes_and_umis_file | sort -u | wc -l)
    echo -e "${barcode}\t${nreads}\t${umis}"
done

# Remove temporary files
rm -f $barcodes_file $barcodes_and_umis_file
##
#
