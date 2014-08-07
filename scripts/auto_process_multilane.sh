#!/bin/sh -e
#
# First attempt at a wrapper to process each lane in a sequencing
# run separately
#
# Get lanes
lanes=$(grep -v ^FCID, custom_SampleSheet.csv | cut -d, -f2 | sort -u)
#
# Get projects
projects=$(grep -v ^FCID, custom_SampleSheet.csv | cut -d, -f10 | sort -u)
#
# Fetch primary data manually
data_dir=$(grep ^data_dir auto_process.info | cut -f2)
echo Source $data_dir
if [ ! -z "$(echo $data_dir | grep '@')" ] || \
    [ !  -z "$(echo $data_dir | grep ':')" ] ; then
    remote=yes
fi
if [ ! -d primary_data ] ; then
    mkdir primary_data
fi
if [ ! -z "$remote" ] ; then
    echo -n Fetching primary data
    auto_process.py make_fastqs --only-fetch-primary-data
else
    echo -n Making link to primary data...
    ln -s $data_dir primary_data/$(basename $data_dir)
    echo done
fi
echo -n Adding primary_data dir to auto_process.info...
sed -i 's/^primary_data_dir\t.*/primary_data_dir\tprimary_data/' auto_process.info
echo done
#
# Set up and run processing for each lane
if [ ! -d processing ] ; then
    mkdir processing
    mkdir processing/logs
fi
for lane in $lanes ; do
    echo -n Start processing lane $lane in background...
    prep_sample_sheet.py \
	--include-lanes=$lane \
	-o processing/custom_SampleSheet.lane$lane.csv \
	custom_SampleSheet.csv \
	>processing/logs/prep_sample_sheet.lane$lane.log
    auto_process.py make_fastqs \
	--no-save \
	--skip-rsync \
	--sample-sheet=processing/custom_SampleSheet.lane$lane.csv \
	--output-dir=processing/bcl2fastq.lane$lane \
	--stats-file=processing/statistics.lane$lane.info \
        --report-barcodes --barcodes-file=processing/index_sequences.lane$lane.txt \
	>processing/logs/$(basename $0).lane$lane 2>&1 &
    echo pid $!
    sleep 2
done
echo -n Waiting for processing to finish...
wait
echo done
#
# Processing completed - rebuild stats files
echo -n Combining stats for individual lanes...
echo "#Project"$'\t'"Sample"$'\t'"Fastq"$'\t'"Size"$'\t'"Nreads"$'\t'"Paired_end" >statistics.info.tmp
cat processing/statistics.lane*.info | grep -v ^#Project >>statistics.info.tmp
LC_COLLATE=C sort statistics.info.tmp >statistics.info
rm -f statistics.info.tmp
echo done
echo -n Regenerate per-lane statistics...
echo "#Lane"$'\t'"Total reads"$'\t'"Assigned reads"$'\t'"Unassigned reads" >per_lane_stats.info
for lane in $lanes ; do
    stats_file=processing/statistics.lane$lane.info
    total=$(grep -v ^# $stats_file | cut -f5 | paste -s -d+ | bc)
    assigned=$(grep -v ^# $stats_file | grep -v ^Undetermined_indices | cut -f5 | paste -s -d+ | bc)
    unassigned=$(grep ^Undetermined_indices $stats_file | cut -f5)
    echo "Lane $lane"$'\t'$total$'\t'$assigned$'\t'$unassigned >>per_lane_stats.info
done
echo done
#
# Put all 'unaligned' outputs into a single directory
echo Collecting outputs into master bcl2fastq directory
mkdir bcl2fastq
for lane in $lanes ; do
    echo Examining lane $lane
    for project_dir in $(ls -d processing/bcl2fastq.lane$lane/Project_* processing/bcl2fastq.lane$lane/Undetermined_indices) ; do
	project=$(basename $project_dir)
	echo Project $project
	mkdir -p bcl2fastq/$project
	for sample_dir in $(ls -d $project_dir/Sample_*) ; do
	    sample=$(basename $sample_dir)
	    echo Sample $sample
	    cp -r $sample_dir bcl2fastq/$project/$sample
	done
    done
done
echo Outputs collected
#
# Remake the projects.info file
echo -n Remaking projects.info...
echo "#Project"$'\t'"Samples"$'\t'"User"$'\t'"Library"$'\t'"Organism"$'\t'"PI"$'\t'"Comments" >projects.info
for project in $projects ; do
    samples=$(grep ^$project$'\t' statistics.info | cut -f2 | paste -s -d,)
    echo "$project"$'\t'"$samples"$'\t'"."$'\t'"."$'\t'"."$'\t'"."$'\t'"." >>projects.info
done
echo done
##
#
