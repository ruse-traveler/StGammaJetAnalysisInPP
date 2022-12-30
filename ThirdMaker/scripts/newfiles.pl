#! /usr/local/bin/perl

use Env;
use File::Basename;

@runs = ();
$nruns = 0;
$runlist = "tmpn.txt";
#"/phenix/workarea/saskia/run4/fast_track/run/prod14.lst";
open(RUNS,"< $runlist");
while(<RUNS>) {
    chomp;
    push(@runs, $_);
    $nruns++;
}

$filelist = "tmpnn.txt";
open(FILES,"< $filelist");
$numfiles = 0;
$noutfile = -1;

$outfile = "files.lst";
open(OUT,"> $outfile");
while(<FILES>) 
{

    $ifile = $_;
    chomp($ifile);
    
    $inlist = 0;
    foreach $runnumber(@runs) 
    {
	if ($runnumber eq $ifile) 
	{
	    $inlist = 1;
	    $numfiles++;
	    last;
	}
    }

    if ($inlist != 1) 
    {
	print OUT "$ifile\n";    
    }
    else {
	print "$ifile\n";
    }
    
}


