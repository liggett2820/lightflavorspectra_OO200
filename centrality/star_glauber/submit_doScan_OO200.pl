#! /usr/bin/perl
#----------------------------------------------------------------------------------------------------
# Step 2 for O+O 200 GeV: submit the NBD (npp, k, x, d) scan.
# One condor job per (npp, k, d); x is scanned inside each job (xbin points).
# Coarse pass below; after it, narrow the ranges around the minimum and run again.
# Usage: ./submit_doScan_OO200.pl        (prints the commands)
#        ./submit_doScan_OO200.pl -run   (submits)
#----------------------------------------------------------------------------------------------------
my $run = (defined $ARGV[0] && $ARGV[0] eq "-run") ? 1 : 0 ;

my $nppbin = 9;  my $nppmin = 1.6;  my $nppmax = 3.2;   # step 0.2
my $kbin   = 5;  my $kmin   = 1.0;  my $kmax   = 9.0;   # step 2.0
my $xbin   = 6;  my $xmin   = 0.08; my $xmax   = 0.28;  # step 0.04 (inside each job)
my $effbin = 4;  my $effmin = 0.00; my $effmax = 0.12;  # step 0.04; d=0 is a constant 98% efficiency

my $npp_step = ($nppbin==1) ? 0 : ($nppmax-$nppmin)/($nppbin-1);
my $k_step   = ($kbin==1)   ? 0 : ($kmax-$kmin)/($kbin-1);
my $eff_step = ($effbin==1) ? 0 : ($effmax-$effmin)/($effbin-1);

my $nJobs = $nppbin * $kbin * $effbin;
print "number of jobs: $nJobs  (each scans $xbin x values)\n";

for (my $i=0; $i<$nppbin; $i++) {
  for (my $j=0; $j<$kbin; $j++) {
    for (my $l=0; $l<$effbin; $l++) {
      my $npp = sprintf("%.3f", $nppmin + $npp_step*$i);
      my $k   = sprintf("%.3f", $kmin   + $k_step*$j);
      my $eff = sprintf("%.3f", $effmin + $eff_step*$l);
      my $submit = "./submit_condor.pl -s doScanX_OO200.csh $npp $k $eff $xbin $xmin $xmax";
      print "$submit\n";
      system("$submit") if $run;
    }
  }
}
