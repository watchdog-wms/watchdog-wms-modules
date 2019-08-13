#!/usr/bin/python3
import argparse
import subprocess
import os.path
import sys
from collections import defaultdict

### parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bedgraphFile", metavar="INPUT_BEDGRAPH_FILE", required=True, help="path to a sorted, not overlapping BEDGRAPH file")
parser.add_argument("-o", "--outputFile", metavar="OUTPUT_BEDGRAPH_FILE", required=True, help="path to the output file")
parser.add_argument("-e", "--expand", default = False, action='store_true', help="expand the ranges instead of shrinking them")
parser.add_argument("-s", "--genomeSize", metavar="INPUT_GENOMESIZE_FILE", required=False, help="path to file containing the size of the contigs")
parser.add_argument("-z", "--addZeroRanges",  default = False, action='store_true', help="adds ranges that are missing with a zero value")
parser.add_argument("-d", "--omitZeroRanges", default = False, action='store_true', help="suppress the output of ranges with a zero value")
args = parser.parse_args()

### var definition
processedChrs = set()

### global vars
lastChr = None
chrSize = dict()
###

### range definition
class Range:
	set = 0
	end = 0
	value = 0

	def __init__(self, start, end, value):
		self.start = start
		self.end = end
		self.value = value

	def length(self):
		return self.end - self.start + 1

	def getOutputFormat(self, chr):
		return '\t'.join([chr, str(self.start), str(self.end), str(self.value)])

	def split(self):
		t = []
		for x in range(self.start, self.end):
			t.append(Range(x, x+1, self.value))
		return t
###

### function definitions
# print error to stderr
def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)
###

# try to fill the buffer from file with maximal n entries from the chr that is currently processed
def fillBuffer(fh, buf, addZeroRanges = False, chrLength = None, n = 500000):
	global lastChr, chrSize

	lastPos = fh.tell()
	line = fh.readline()
	artificial = 0
	end = -1

	while line:
		tmp = line.strip().split("\t")
		chr = tmp[0]

		# check, if the same chr is processed
		if lastChr == chr or lastChr is None:
			start = int(tmp[1])
			end = int(tmp[2])
			value = int(float(tmp[3]))

			# add zero start region
			if addZeroRanges and lastChr is None:
				buf.append(Range(0, start, 0))
				artificial += 1

			# add region that should be added
			buf.append(Range(start, end, value))

			lastChr = chr
			lastPos = fh.tell()

			# test if buffer is full
			if len(buf) >= n:
				break
			else:
				# read next line
				line = fh.readline()
		else:
			fh.seek(lastPos)
			break

	# EOF was reached
	if len(line) == 0:
		print("Reached EOF")
		return False

	# test if some data was read
	if (len(buf)-artificial) > 0:
		print("loaded %i ranges from sample on %s" % ((len(buf)-artificial), lastChr))
	# EOF was not reached
	return True
###

# tests, if two ranges are consecutive
# if not a spacing element is returned
def getConsecutiveRanges(f, l):

	# test if the first two are connected
	if f.end != l.start:
		return Range(f.end, l.start, 0)

	return None
###


# writes a range to a file
def writeRangeToFile(org, chrom, fh, expand):
	# if in expand mode 
	if expand:
		newRanges = org.split()
		# ouput all of them
		for r in newRanges:
			outString = r.getOutputFormat(chrom)
			outFH.write(outString)
			outFH.write("\n")

	else:
			outString = org.getOutputFormat(chrom)
			outFH.write(outString)
			outFH.write("\n")
###

#######################################################
################# START OF MAIN PROGRAMM ##############
#######################################################

if args.addZeroRanges and args.omitZeroRanges:
	eprint("Parameters --addZeroRanges and --omitZeroRanges are exclusive.")
	exit(1)

print("mode: %s" %("expand" if args.expand else "shrink"))
print("add zero ranges: %s" %("yes" if args.addZeroRanges else "no"))
print("omit zero ranges: %s" %("yes" if args.omitZeroRanges else "no"))

### read in the genome size file
if not args.genomeSize is None:
	if os.path.isfile(args.genomeSize) and os.access(args.genomeSize, os.R_OK):
		with open(args.genomeSize, "r") as fh:
			print("genome size file: '%s'" % (args.genomeSize))

			for line in fh:
				tmp = line.strip().split("\t")
				chr = tmp[0]
				size = int(tmp[1])
				chrSize[chr] = size

		#########################################################
	else:
		eprint("genome size file '%s' was not found or is not readable." % (args.genomeSize))
		exit(1)


### open the bedgraph file
if os.path.isfile(args.bedgraphFile) and os.access(args.bedgraphFile, os.R_OK):
	fh = open(args.bedgraphFile, "r+")
	print("bedgraph file: '%s" % (args.bedgraphFile))
else:
	eprint("Bedgraph file '%s' was not found or is not readable." % (args.bedgraphFile))
	exit(1)


### check if base output folder exists
dirbase = os.path.dirname(args.outputFile)
if len(dirbase) > 0 and not os.path.exists(dirbase):
	os.makedirs(dirbase)

### open the output file and start with processing
with open(args.outputFile, "w") as outFH:
	print("output file: %s" % (args.outputFile))
	buffer = []
	startShrinkRange = None
	curRange = None
	nextRange = None
	lastChrIteration = True
	maxEnd = 0
	EOF = False
	while True:
		if len(buffer) > 0:
			# just for the case that len(buffer) is only 1
			if len(buffer) == 1:
				nextRange = buffer[0]
			else:
				# loop over ranges
				for i in range(0, (len(buffer)-1)):
					curRange = buffer[i]
					nextRange = buffer[i+1]

					# ensure that sorting is given
					if curRange.start < maxEnd:
						eprint("Bedgraph file is either not sorted or has overlapping regions: %s-%s:%s" % (lastChr, curRange.start, curRange.end))
						exit(1)

					# ommit zero ranges
					if (not args.omitZeroRanges) or curRange.value != 0:

						# test if a spacer must be added					
						spacer = getConsecutiveRanges(curRange, nextRange)

						# if in expand mode 
						if args.expand:
							writeRangeToFile(curRange, lastChr, outFH, True)
						# shrink mode					
						else:
							# test if both have the same value and no spacer is between them!
							if spacer is None and curRange.value == nextRange.value:
								if startShrinkRange is None:
									startShrinkRange = curRange
							else:
								# test if a merged range mustbe created
								if startShrinkRange is None:
									writeRangeToFile(curRange, lastChr, outFH, False)
								else:
									shrinkedRange = Range(startShrinkRange.start, curRange.end, curRange.value)
									writeRangeToFile(shrinkedRange, lastChr, outFH, False)
									startShrinkRange = None

						# output the spacer
						if args.addZeroRanges and not spacer is None:
							writeRangeToFile(spacer, lastChr, outFH, args.expand)

					# store max end of range
					if maxEnd < curRange.end:
						maxEnd = curRange.end
				
				### end of for loop
			# clear buffer
			buffer = []
		else:
			# EOF was reached
			if not EOF and not fillBuffer(fh, buffer, args.addZeroRanges):
				EOF = True
			
			# if buffer is empty, processing of chr has ended				
			if len(buffer) == 0:
				# add nextRange and dummy element
				if lastChrIteration:
					buffer = [nextRange]

					# test if last region should be added
					if args.addZeroRanges and lastChr in chrSize:
						le = chrSize[lastChr]
						if le > nextRange.end:
							buffer.append(Range(nextRange.end, le, 0))
						else:
							print("[WARN] Chromosome length of chr %s is longer than specified in the genome size file" % lastChr)

					# add dummy region that causes last region to be outputted
					if len(buffer) > 0:
						tmpRange = buffer[len(buffer)-1]
						buffer.append(Range(tmpRange.end, tmpRange.end+1, None))

				# start processing of new chr
				else:
					processedChrs.add(lastChr)
					lastChr = None
					startShrinkRange = None
					nextRange = None
					maxEnd = 0

					# end processing as chr is processed completely and file is empty
					if EOF:
						break

				# change value of last chr iteration
				lastChrIteration = False if lastChrIteration else True
			else:
				# add nextRange as it was not processed on current chr yet 
				if not nextRange is None:
					buffer.insert(0, nextRange)

		
### close the open file handles
fh.close()

### all was ok :)
print("Finished processing of %i chromosomes!" % (len(processedChrs)))
exit(0)
