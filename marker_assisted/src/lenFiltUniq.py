
######################################################################
#  PUBLIC DOMAIN NOTICE
#
#  This software is "United States Government Work" under the terms of the United
#  States Copyright Act. It was written as part of the authors' official duties
#  for the United States Government and thus cannot be copyrighted. This software
#  is freely available to the public for use without a copyright
#  notice. Restrictions cannot be placed on its present or future use.
#
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and associated data, the National Human Genome
#  Research Institute (NHGRI), National Institutes of Health (NIH) and the
#  U.S. Government do not and cannot warrant the performance or results that may
#  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
#  disclaim all warranties as to performance, merchantability or fitness for any
#  particular purpose.
#
#  Please cite the authors in any work or product based on this material.
######################################################################

# import statements
import sys

# make sure we have enough input (and quit if we don't)
if len(sys.argv) != 3:
	print(f"ERROR: {sys.argv[0]} expected 2 arguments: "
			"1. the length cutoff measured in Kb (e.g., 25). "
			"2. the percent identity cutoff in percent (i.e., 75). "
			"Integers only for both arguments. "
			"Input (the output from samToErrorRate) will be read from STDIN. "
			"Output will be written to STDOUT. The output will contain a single "
			"column of read (i.e., query) identifiers (the values in column "
			"number 1 (1-based) from the input). Output identifiers will not be "
			"duplicated, and each will be output only if one or more alignment "
			"records associated with the respective identifier passes the "
			"cutoffs provided as positional arguments to this program.", file=sys.stderr)
	sys.exit(1)

# grab the cutoffs from command line args
length_cutoff = int(sys.argv[1]) * -1   # e.g., 25 (kb)
pid_cutoff = int(sys.argv[2])           # e.g., 75 (%)

# initialize the set
rids = set() # "rids" means "read identifiers", i.e., the query identifiers from the alignments

# read through the file (i.e., the output from samToErrorRate)
for line in sys.stdin:
	# split line into fields and extract relevant ones
	fields = line.rstrip('\n').split('\t')
	rid = fields[0] # rid: "read id", i.e., query id
	score = int(fields[2])
	percent_identity = float(fields[3])

	# add rid to Set if passing cutoffs
	if score <= length_cutoff and percent_identity >= pid_cutoff:
		rids.add(rid)


# write output list
print('\n'.join(list(rids)))
#print('\n'.join(sorted(list(rids)))) # alternatively, sort the rids first (sorting not needed for SubFile program)

#############
### NOTES ###
#############
#
# Notes last updated on 22 March 2022
#
# Previously, filter_by_marker_nosplit.sh (line 130) did this:
#
#	$SCRIPT/src/samToErrorRate $target.filtered.sam $asm \
#		| awk -v l=$len_filt 'BEGIN{rid=""}{if($3 <= -1*l && $4 >= 75 && rid != $1){print $1; rid=$1}}' \
#		> $target.filteredList
#
# But $target.filtered.sam was not sorted by read ids. To avoid sorting, this
# script was implemented, and it is called like this:
#
#	$SCRIPT/src/samToErrorRate $target.filtered.sam $asm \
#		| python3 $SCRIPT/src/lenFiltUniq.py $len_filt 75 \
#		> $target.filteredList
#
# Internally, the cutoff checks performed in awk are now performed with this
# code (see lines 56-57):
#
#	if score <= length_cutoff and percent_identity >= pid_cutoff:
#		rids.add(rid)
#
# The variable "score" is $3 from awk. "length_cutoff" is -1*l from awk.
# "percent_identity" is $4 from awk. "pid_cutoff" was hardcoded as 75 in the awk
# script, but it is provided to this program via a command-line argument.
#
