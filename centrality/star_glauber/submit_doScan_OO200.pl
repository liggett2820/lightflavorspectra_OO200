#! /usr/bin/perl
#----------------------------------------------------------------------------------------------------
# Step 2 for O+O 200 GeV: submit the NBD (npp, k, x, d) scan.
# One condor job per (npp, k, d); x is scanned inside each job (xbin points).
# Coarse pass below; after it, narrow the ranges around the minimum and run again.
# Usage: ./submit_doScan_OO200.pl        (prints the commands)
#        ./submit_doScan_OO200.pl -run   (submits)
#        add -fine for the finer grid around the coarse minimum (move LOG_Scan aside first)
#        add -local -multcut 10 (or 20) for the fit-range check around the fine minimum (logs in LOG_Scan_mc10/20)
#----------------------------------------------------------------------------------------------------
my $run  = (grep { $_ eq "-run" }  @ARGV) ? 1 : 0 ;
my $fine = (grep { $_ eq "-fine" } @ARGV) ? 1 : 0 ;   # finer grid around the coarse-scan minimum
my $local = (grep { $_ eq "-local" } @ARGV) ? 1 : 0 ; # small grid around the fine-scan minimum (for the multCut checks)
my $multcut = 15;                                      # fit range refMult >= multcut; -multcut N to change
for (my $a=0; $a<=$#ARGV; $a++) { $multcut = $ARGV[$a+1] if ($ARGV[$a] eq "-multcut" && defined $ARGV[$a+1]); }
my $logdir = ($multcut == 15) ? "LOG_Scan" : "LOG_Scan_mc$multcut";   # keep each fit range's logs apart

my $nppbin = 9;  my $nppmin = 1.6;  my $nppmax = 3.2;   # step 0.2
my $kbin   = 5;  my $kmin   = 1.0;  my $kmax   = 9.0;   # step 2.0
my $xbin   = 6;  my $xmin   = 0.08; my $xmax   = 0.28;  # step 0.04 (inside each job)
my $effbin = 4;  my $effmin = 0.00; my $effmax = 0.12;  # step 0.04; d=0 is a constant 98% efficiency

# -fine: around the coarse-scan best point (npp, k, x, d) = (2.8, 3, 0.16, 0.08), chi2/ndf 2.76 (Oct 6 2026)
if ($fine) {
  ($nppbin, $nppmin, $nppmax) = (9, 2.60, 3.00);   # step 0.05
  ($kbin,   $kmin,   $kmax)   = (5, 2.0,  4.0);    # step 0.5
  ($xbin,   $xmin,   $xmax)   = (9, 0.12, 0.20);   # step 0.01 (inside each job)
  ($effbin, $effmin, $effmax) = (5, 0.04, 0.12);   # step 0.02
}
# -local: around the fine-scan best point (2.75, 2.5, 0.17, 0.06), chi2/ndf 0.89 (Oct 6 2026); 45 jobs
if ($local) {
  ($nppbin, $nppmin, $nppmax) = (5, 2.65, 2.85);   # step 0.05
  ($kbin,   $kmin,   $kmax)   = (3, 2.0,  3.0);    # step 0.5
  ($xbin,   $xmin,   $xmax)   = (7, 0.14, 0.20);   # step 0.01 (inside each job)
  ($effbin, $effmin, $effmax) = (3, 0.04, 0.08);   # step 0.02
}
if ($run && !-d $logdir) { mkdir($logdir) or die "cannot create $logdir\n"; }

my $npp_step = ($nppbin==1) ? 0 : ($nppmax-$nppmin)/($nppbin-1);
my $k_step   = ($kbin==1)   ? 0 : ($kmax-$kmin)/($kbin-1);
my $eff_step = ($effbin==1) ? 0 : ($effmax-$effmin)/($effbin-1);

my $nJobs = $nppbin * $kbin * $effbin;
print "number of jobs: $nJobs  (each scans $xbin x values), fit refMult >= $multcut, logs in $logdir\n";

for (my $i=0; $i<$nppbin; $i++) {
  for (my $j=0; $j<$kbin; $j++) {
    for (my $l=0; $l<$effbin; $l++) {
      my $npp = sprintf("%.3f", $nppmin + $npp_step*$i);
      my $k   = sprintf("%.3f", $kmin   + $k_step*$j);
      my $eff = sprintf("%.3f", $effmin + $eff_step*$l);
      my $extra  = ($multcut == 15) ? "" : " $multcut";   # default jobs keep their old arguments and log names
      my $submit = "./submit_condor.pl -s -l $logdir doScanX_OO200.csh $npp $k $eff $xbin $xmin $xmax$extra";
      print "$submit\n";
      system("$submit") if $run;
    }
  }
}
