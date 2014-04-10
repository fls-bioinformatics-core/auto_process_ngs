#!/bin/sh
#
# Manage fastqs:
# - report individal file sizes
# - create checksums
# - make zip archive
# - copy to another location
#
if [ $# -lt 1 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ] ; then
    echo Usage $(basename $0) DIR PROJECT CMD
    echo Optional CMD can be:
    echo zip: make zip archive of fastq.gzs
    echo md5sums: generate MD5sums
    echo scp LOCATION: scp fastq.gzs to LOCATION
    exit
fi
# Analysis dir
DIR=$1
# Name of project
PROJECT=$2
# Command
CMD=$3
case "$CMD" in
    zip)
	do_md5sums=yes
	do_zip=yes
	;;
    scp)
	do_md5sums=yes
	do_scp=yes
	LOCATION=$4
	;;
    md5sums)
	do_md5sums=yes
	;;
    *)
	;;
esac
# 
# Do checks
echo -n Checking for analysis dir $DIR...
if [ ! -d "$DIR" ] ; then
    not found
    echo FAILED: directory $DIR not found >&2
    exit 1
fi
echo ok
echo -n Checking for data dir...
BCL2FASTQ_DIR=
for d in bcl2fastq unaligned ; do
    if [ -z "$BCL2FASTQ_DIR" ] ; then
	if [ -d $DIR/$d ] ; then
	    BCL2FASTQ_DIR=$d
	    break
	fi
    fi
done
if [ -z "$BCL2FASTQ_DIR" ] ; then
    echo not found
    echo FAILED: no data dir found >&2
    exit 1
fi
echo $BCL2FASTQ_DIR
# List all project names
if [ -z "$PROJECT" ] ; then
    echo Projects:
    find $DIR/$BCL2FASTQ_DIR -maxdepth 1 -name "Project_*" -type d -exec basename {} \; | sed 's/^Project_//'
    exit
fi
echo -n Checking for project $PROJECT...
if [ ! -d $DIR/$BCL2FASTQ_DIR/Project_$PROJECT ] ; then
    echo not found
    echo FAILED: cannot find project $PROJECT >&2
    exit 1
fi
echo ok
echo -n Checking for stats...
STATS=$DIR/statistics.info
if [ ! -f $STATS ] ; then
    echo not found
    echo FAILED: no statistics.info file found >&2
    exit 1
fi
echo ok
# File sizes
grep ^$PROJECT $STATS | cut -f2-4
mb_size=$(find $DIR/$BCL2FASTQ_DIR/Project_$PROJECT -name "*.fastq.gz" -exec du -B1048576 {} \; | cut -f1 | paste -s -d"+" | bc)
if [ $mb_size -gt 1000 ] ; then
    size=$(echo "scale=1; $mb_size.0/1000.0" | bc)
    units=G
else
    size=${mb_size}
    units=M
fi
echo Total$'\t'$size$units
# Other info
n_files=$(grep ^$PROJECT $STATS | wc -l)
n_R2=$(grep ^$PROJECT $STATS | grep _R2_ | wc -l)
echo $n_files fastqs
if [ $n_R2 -gt 0 ] ; then
    echo Paired end data
fi
# Generate MD5 checksums
if [ ! -z "$do_md5sums" ] ; then
    MD5SUMS=$PROJECT.chksums
    if [ ! -f $MD5SUMS ] ; then
	echo Making checksums file $MD5SUMS
	if [ -f $MD5SUMS.part ] ; then
	    /bin/rm -f $MD5SUMS.part
	fi
	find $DIR/$BCL2FASTQ_DIR/Project_$PROJECT -name "*.fastq.gz" -exec md5sum {} \; > $MD5SUMS.part
	sed -i 's/\/.*\/Sample_.*\///g' $MD5SUMS.part
	/bin/mv $MD5SUMS.part $MD5SUMS
	echo Done
    else
	echo Checksums file $MD5SUMS already exists
    fi
fi
# Make zip file
if [ ! -z "$do_zip" ] ; then
    ZIP_FILE=$PROJECT.zip
    if [ ! -f $ZIP_FILE ] ; then
	echo Making zip file $ZIP_FILE
	if [ -f $ZIP_FILE.part ] ; then
	    /bin/rm -f $ZIP_FILE.part
	fi
	find $DIR/$BCL2FASTQ_DIR/Project_$PROJECT -name "*.fastq.gz" -exec zip -gj $ZIP_FILE.part {} \;
	zip -gj $ZIP_FILE.part $MD5SUMS
	/bin/mv $ZIP_FILE.part $ZIP_FILE
    else
	echo Zip file $ZIP_FILE already exists
    fi
fi
# Copy files with scp
if [ ! -z "$do_scp" ] ; then
    if [ -z "$LOCATION" ] ; then
	echo "Can't scp: no location supplied" >&2
	exit 1
    fi
    echo Scping fastq.gz files to $LOCATION
    find $DIR/$BCL2FASTQ_DIR/Project_$PROJECT -name "*.fastq.gz" -exec scp {} $LOCATION \;
    scp $MD5SUMS $LOCATION
    echo Done
fi
##
#
	
