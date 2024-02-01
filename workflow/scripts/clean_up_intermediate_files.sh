#!/usr/bin/bash

# delete the intermediate files from smf pipeline to save space
# in particular, sometimes it is useful to inspect the sam files from the all_alignments
# but these take up mad space
# run this from the homedir (i.e. one above results)

rm results/*/*/tmp/*all_alignments.sam
rm results/*/*/tmp/*problematic_reads.sam
