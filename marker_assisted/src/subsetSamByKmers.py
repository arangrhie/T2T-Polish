import sys
import re
from sets import Set

f = open('%s'% (sys.argv[1]))
rids = {}
weights={}
start=0
end=0

for line in f:
#A00428:24:H5327DSXX:2:1101:10981:2566_0_0_58121944 9
   l = line.strip().split()
   weight=int(l[1])
   l = l[0].split("_")
   pos=int(l[-2])
   if start == 0 or start > pos:
      start=pos
   if end == 0 or end < pos:
      end=pos

   rid="%s_%s"%("_".join(l[0:-4]), l[-3])
   if rid not in rids:
      rids[rid]=pos
      weights[rid]=weight
      print >> sys.stderr, "Read %s at pos %d got input from score 0 because of %d, diff of %d"%(rid, pos, weight, weight)
   else:
      currWeight = weights[rid]
      if currWeight < weight:
         rids[rid] = pos
         weights[rid]=weight
         print >> sys.stderr, "Read %s got bumped to pos %d from score %d because of %d, diff of %d"%(rid, pos, currWeight, weight, weight-currWeight)
       
print >> sys.stderr, "Done readind %d records covering %d-%d"%(len(rids), start, end)
f.close()

f = open('%s'% (sys.argv[2]))
#A00428:24:H5327DSXX:2:2328:6216:10003	147	chrX_bothkpatchedin	57679805	60	39M1D83M28S	=	57679722	-206	GTAGAATCTACAAGTGGACATTTGGGGCTCTTTGAGGCCTATGTTGACAAAGGAACTATCTTGTCATGAAAACTACACAGAATCATTCTCAGAAACTTCTTTGTGATGTGTGCGTTCAACTCTGTGAGTTGAATGCACACATCACAAAGA	F:FFFFFFF:FF,FFFFFFFFFFFFFFF:F:FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:F:FF::FFFFFFFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	RX:Z:AAGGCTTGGTTCCTAA	QX:Z:FFFFFFFFFFFFFFFF	BC:Z:CTGTAACT	QT:Z:FFFFFFFF	XS:i:-79	AS:i:-27	XM:Z:0	AM:Z:0	XT:i:0	SA:Z:tig00000178|arrow|arrow,407355,+,38M112H,0,0;	RG:Z:CHM13_prep5:LibraryNotSpecified:1:unknown_fc:0	OM:i:60

for line in f:
   # output header lines
   if line.startswith("@"):
      print line.strip()
   else:
      l = line.strip().split()
      isReverse=int(l[1]) & 16 > 0
      # if this maping is outside of our defined region, output all of them
      #if (int(l[3]) < start or int(l[3]) > end):
         #print line.strip()
         # do nothing
      #else:
      if isReverse:
         rid="%s_%d"%(l[0], 1)
      else:
         rid="%s_%d"%(l[0], 0)
      if rid in rids:
         pos=rids[rid]
         if pos == int(l[3]):
            l[4]="60"
            print "\t".join(l)
f.close()
