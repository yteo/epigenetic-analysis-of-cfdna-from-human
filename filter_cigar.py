import pysam
import itertools
import sys
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
passed = pysam.AlignmentFile(sys.argv[2], "wb", template=samfile)
for read in samfile.fetch():
	a=[(sum(i[1] for i in group), key) for key, group in itertools.groupby(sorted(read.cigar, key = lambda i: i[0]), lambda i: i[0])]
	out_tup = [i for i in a if i[1] == 4]
	#print out_tup
	if len(out_tup)==0:
		passed.write(read)
	elif [i for i in out_tup if i[0] < 75]:
		passed.write(read)
passed.close()
samfile.close()