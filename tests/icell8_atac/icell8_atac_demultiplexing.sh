#!/bin/bash -e
#
# Run ICELL8 ATAC demultiplexing on test data

TEST_DATA_DIR=$(dirname $0)/test-data
WELL_LIST=${TEST_DATA_DIR}/well_list.txt

# Clean up output directories
OUTDIRS=\
"demultiplex.icell8_atac.samples
demultiplex.icell8_atac.barcodes
no_demultiplexing.icell8_atac"
        
for d in $OUTDIRS ; do
    rm -rf $d
done

# Demultiplex into samples
demultiplex_icell8_atac.py \
    --mode=samples \
    --batch_size=2 \
    --swap-i1-i2 \
    --reverse-complement=i1 \
    --update-read-headers \
    --output-dir demultiplex.icell8_atac.samples \
    ${WELL_LIST} ${TEST_DATA_DIR}/*.fastq.gz

# Demultiplex into barcodes
demultiplex_icell8_atac.py \
    --mode=barcodes \
    --batch_size=2 \
    --swap-i1-i2 \
    --reverse-complement=i1 \
    --update-read-headers \
    --output-dir demultiplex.icell8_atac.barcodes \
    ${WELL_LIST} ${TEST_DATA_DIR}/*.fastq.gz

# No demultiplexing
demultiplex_icell8_atac.py \
    --batch_size=2 \
    --swap-i1-i2 \
    --reverse-complement=i1 \
    --output-dir no_demultiplexing.icell8_atac \
    --no-demultiplexing \
    ${WELL_LIST} ${TEST_DATA_DIR}/*.fastq.gz
##
#
