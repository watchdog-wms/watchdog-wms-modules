#! /usr/bin/perl -w

$dir=$ARGV[0];
$mergeDist=$ARGV[1];
$minLength=$ARGV[2];

$minScore=0;
if(@ARGV>3)
{
   $minScore=$ARGV[3];
}

$out=$dir."_combined_".$mergeDist."_".$minLength."_".$minScore.".bed";
$out=~s/\/_combined/_combined/;


opendir (DIR, $dir) or die $!;

open(OUT, ">".$out);


while (my $file = readdir(DIR)) 
{
    $prevmax=-1;
    $prevstart=-1;
    $prevend=-1;
    $previd="";
    
    open(IN, $dir."/".$file);
    @file=<IN>;
    close(IN);
    
    foreach (@file)
    {
       chomp();
       @line=split(/\t/, $_);
       
       if($prevmax==-1 || $line[1]-$prevend> $mergeDist)
       {
         if($prevmax!=-1 && $prevend-$prevstart+1>=$minLength&& $prevmax>=$minScore)
         {
	    print OUT $line[0], "\t", $prevstart, "\t", $prevend, "\t",  $previd, "\t", $prevmax, "\n", 
         }
        
	 $prevstart=$line[1];
	 $prevend=$line[2];
	 $previd=$line[3];
	 $prevmax=$line[4];
	 
       }
       else{
	 
	  $prevend=$line[2];
	  $previd.=",".$line[3];
	  if($line[4]>$prevmax)
	  {
	    $prevmax=$line[4];
	  }
       }
       
    }
    if($prevmax!=-1 && $prevend-$prevstart+1>=$minLength && $prevmax>=$minScore)
    {
      print OUT $line[0], "\t", $prevstart, "\t", $prevend, "\t",  $previd, "\t", $prevmax, "\n", 
    }
    

}
close(OUT);