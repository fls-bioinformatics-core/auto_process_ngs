Structure of Analysis Directories
=================================

Top-level analysis directory
****************************

The top-level analysis directory is created by the ``setup`` command:

::

 140729_SN1234_XYZ_analysis/            Top-level analysis directory
          |
          +--- auto_process.info        Parameters used for processing
          |
          +--- statistics.info          Per-file statistics (size, no. of reads etc)
          |
          +--- per_lane_stats.info      Basic per-lane statistics
          |
          +--- (README.txt)             General README file with info on unusual processing
          |
          +--- SampleSheet.orig.csv     Original sample sheet file from sequencer
          |
          +--- custom_SampleSheet.csv   Converted/edited sample sheet used for FASTQ generation
          |
          +--- projects.info            Info on each project, used to generate project subdirs
          |
          +--- ScriptCode/
          |
	  +--- primary_data/            Temporary
	  |
	  +--- logs/
	  |
	  +--- bcl2fastq/               Output from CASAVA/bclToFastq conversion
	  |           |
          |           +--- Project_PeterBriggs/
          |           |            |
          |           |            +--- Sample_PB1/
          |           |            |           |
          |           |            |           +--- [Gzipped FASTQs]
          |           |            |
          |           |            +--- Sample_PB2/
          |           |            |
          |           |           ...
	  |           |
          |           +--- Project_AnotherUser/
          |           |            |
          |           |           ...
          |           |
          |          ...
          |
          +--- PeterBriggs/             Analysis project directory for 'PeterBriggs' project
          |           |
          |          ...
          |
          +--- AnotherUser/             Analysis project directory for 'AnotherUser' project
          |           |
          |          ...
          |
          +--- _undetermined/           Analysis project directory for undetermined reads
          |           |
          |          ...


Analysis Project Subdirectories
*******************************

The analysis project subdirectories are created within the top-level analysis
directory by the ``setup_analysis_dirs`` command:

::

 140729_SN1234_XYZ_analysis/            Top-level analysis directory
          |
          .
          .
          .
          |
          +--- PeterBriggs/             Analysis project directory for 'PeterBriggs' project
          |           |
          |           +--- README.info  Metadata for project (user, PI, organism etc)
          |           |
          |           +--- fastqs/
          |           |        |
          |           |        +--- [Gzipped FASTQS - hard links to 'bcl2fastq' copies]
          |           |
          |           +--- qc/
          |           |
          |           +--- ScriptCode/
          |           |
          |           +--- logs/
          |
	  |
          +--- AnotherUser/             Analysis project directory for 'AnotherUser' project
          |           |
          |          ...

