## By Alex Dobin (modified by Jakub Orzechowski Westholm & Sophie ;) )
## Filters sam output for circular RNAs.
##
## Example usage: $ cat Chimeric.out.junction | awk -f filterCirc.awk | sort | uniq -c | sort -k1,1rn > circCollapsed.txt
## original version -> circRNA_finder

function cigarGenomicDist(cig)
{
    n=split(cig,L,/[A-Z]/)-1;
    split(cig,C,/[0-9]*/);
    g=0;
    for (ii=1;ii<=n;ii++) {//scan through CIGAR operations
	    if (C[ii+1]!="S" && C[ii+1]!="I") {
		g+=L[ii];
	    };
    };
    return g;
};
BEGIN {
    endTol=5;
};
{
    ### junction is outside mates AND on same chr AND on same strand AND read pair < 100Kbp apart in the right orientation
    if ( $7>=0 && $1==$4 && $3==$6 && (($3=="-" && $5>$2 && $5-$2<100000) || ($3=="+" && $2>$5 && $2-$5<100000)) )
    {
	## non-overlapping reads in pair.
	if ( ($3=="+" && $11+endTol>$5 && $13+cigarGenomicDist($14)-endTol<=$2) \
	     || ($3=="-" && $13+endTol>$2 && $11+cigarGenomicDist($12)-endTol<=$5) ) {
	     ## output: chrom, start = smaller coordinate, end = bigger coordinate -1 (0-based exclusive circ start and end), strand (STAR output = strand of the first read of a pair),
	     ## read_id, splice signal 0=no signal 1=on + strand 2=on - strand
	    print $1,($3=="+"?$5:$2),($3=="+"?$2-1:$5-1),$3,$10,$7;
	};
    };
};
