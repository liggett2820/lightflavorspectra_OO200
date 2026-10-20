#! /usr/bin/perl
#----------------------------------------------------------------------------------------------------
# Step 2 for O+O 200 GeV: submit the NBD (npp, k, x, d) scan.
# One condor job per (npp, k, d); x is scanned inside each job (xbin points).
# Coarse pass below; after it, narrow the ranges around the minimum and run again.
# Usage: ./submit_doScan_OO200.pl        (prints the commands)
#        ./submit_doScan_OO200.pl -run   (submits)
#        add -fine for the finer grid around the coarse minimum (move LOG_Scan aside first)
#        add -local -multcut 10 (or 20) for the fit-range check around the fine minimum (logs in LOG_Scan_mc10/20)
#        add -wide for the extended grid used on the corrected refMult (|Vz| <= 30 cm) histogram
#        add -totnmip -fine for the finer totnMIP grid around the coarse minimum (logs in LOG_Scan_totnmip_mcN_fine)
#        add -totnmip for the EPD totnMIP scan (data hTotnMIPCorr_OO200.root, fit totnMIP >= 200 by default,
#             logs in LOG_Scan_totnmip_mcN); combine with -multcut 150 / 250 for the fit-range check
#----------------------------------------------------------------------------------------------------
my $run  = (grep { $_ eq "-run" }  @ARGV) ? 1 : 0 ;
my $fine = (grep { $_ eq "-fine" } @ARGV) ? 1 : 0 ;   # finer grid around the coarse-scan minimum
my $local = (grep { $_ eq "-local" } @ARGV) ? 1 : 0 ; # small grid around the fine-scan minimum (for the multCut checks)
my $wide  = (grep { $_ eq "-wide" }  @ARGV) ? 1 : 0 ; # extended grid for the corrected refMult, |Vz| <= 30 cm
my $totnmip = (grep { $_ eq "-totnmip" } @ARGV) ? 1 : 0 ; # EPD totnMIP instead of refMult
my $multcut = $totnmip ? 200 : 15;                     # fit range mult >= multcut; -multcut N to change
for (my $a=0; $a<=$#ARGV; $a++) { $multcut = $ARGV[$a+1] if ($ARGV[$a] eq "-multcut" && defined $ARGV[$a+1]); }
my $logdir = ($multcut == 15) ? "LOG_Scan" : "LOG_Scan_mc$multcut";   # keep each fit range's logs apart
$logdir = "LOG_Scan_totnmip_mc$multcut" if $totnmip;
my $datafile = "hTotnMIPCorr_OO200.root";              # used only with -totnmip
if ($run && $totnmip && !-e $datafile) { die "-totnmip: $datafile not found here; copy it from the FitEpdCent_OO200.C output\n"; }

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
# -wide: corrected refMult, |Vz| <= 30 cm (Oct 13 2026). The -local scan of it put the best point (2.65, 2.5, 0.20, 0.04),
# chi2/ndf 1.00, on the npp, x and d edges, along the npp-x valley; extend all three. 135 jobs, 13 x values each.
if ($wide) {
  ($nppbin, $nppmin, $nppmax) = (9,  2.45, 2.85);  # step 0.05
  ($kbin,   $kmin,   $kmax)   = (3,  2.0,  3.0);   # step 0.5
  ($xbin,   $xmin,   $xmax)   = (13, 0.14, 0.26);  # step 0.01 (inside each job)
  ($effbin, $effmin, $effmax) = (5,  0.00, 0.08);  # step 0.02
}
# -totnmip: EPD totnMIP, coarse pass (Oct 20 2026). Starting guess: with O+O Glauber <Npart> ~ 9.5, <Ncoll> ~ 12
# (rough), <totnMIP> ~ 150 over the full cross section needs npp ~ 32 at x = 0 and ~ 24 at x = 0.2 (the same
# estimate gives <refMult> ~ 16 for the refMult best point). x starts at 0 because the EPD (2.1 < |eta| < 5.1)
# sits in the fragmentation region, where wounded-nucleon (x = 0) scaling is expected; k goes low because the
# Landau spread of nMIP per particle widens the per-source distribution, which the fit can only absorb in k.
# Fit range totnMIP >= 200: there 860003's efficiency (official min-bias weight x measured 860003/min-bias ratio)
# is within ~1-2% of its plateau. 10 x 5 x 3 = 150 jobs, 7 x values each.
# RatioChi2Files/ names carry only (npp, k, x, d): run the -multcut checks after this scan, not at the same time.
if ($totnmip) {
  ($nppbin, $nppmin, $nppmax) = (10, 12.0, 39.0);  # step 3
  ($kbin,   $kmin,   $kmax)   = (5, 0.5,  12.5);   # step 3
  ($xbin,   $xmin,   $xmax)   = (7, 0.00, 0.30);   # step 0.05 (inside each job)
  ($effbin, $effmin, $effmax) = (3, 0.00, 0.08);   # step 0.04; eff = 0.98 (1 - d M/540), d = 0 is constant
}
# -totnmip -fine: around the coarse best point (npp, k, x, d) = (24, 3.5, 0.25, 0.08), chi2/ndf 455.6/432 = 1.055
# (Oct 20 2026). The coarse top five lie along an npp-x valley from (24, 0.25) to (30, 0.05), and d sat on its
# upper edge, so: npp 22-31, k 2-6, x 0-0.30, d up to 0.16. 10 x 5 x 5 = 250 jobs, 13 x values each.
if ($totnmip && $fine) {
  ($nppbin, $nppmin, $nppmax) = (10, 22.0, 31.0);  # step 1
  ($kbin,   $kmin,   $kmax)   = (5, 2.0,  6.0);    # step 1
  ($xbin,   $xmin,   $xmax)   = (13, 0.00, 0.30);  # step 0.025 (inside each job)
  ($effbin, $effmin, $effmax) = (5, 0.00, 0.16);   # step 0.04
  $logdir .= "_fine";
}
if ($run && !-d $logdir) { mkdir($logdir) or die "cannot create $logdir\n"; }

my $npp_step = ($nppbin==1) ? 0 : ($nppmax-$nppmin)/($nppbin-1);
my $k_step   = ($kbin==1)   ? 0 : ($kmax-$kmin)/($kbin-1);
my $eff_step = ($effbin==1) ? 0 : ($effmax-$effmin)/($effbin-1);

my $nJobs = $nppbin * $kbin * $effbin;
my $var = $totnmip ? "totnMIP" : "refMult";
print "number of jobs: $nJobs  (each scans $xbin x values), fit $var >= $multcut, logs in $logdir\n";

for (my $i=0; $i<$nppbin; $i++) {
  for (my $j=0; $j<$kbin; $j++) {
    for (my $l=0; $l<$effbin; $l++) {
      my $npp = sprintf("%.3f", $nppmin + $npp_step*$i);
      my $k   = sprintf("%.3f", $kmin   + $k_step*$j);
      my $eff = sprintf("%.3f", $effmin + $eff_step*$l);
      my $extra  = ($multcut == 15) ? "" : " $multcut";   # default jobs keep their old arguments and log names
      $extra = " $multcut $datafile" if $totnmip;          # totnMIP jobs always pass the fit range and the data file
      my $submit = "./submit_condor.pl -s -l $logdir doScanX_OO200.csh $npp $k $eff $xbin $xmin $xmax$extra";
      print "$submit\n";
      system("$submit") if $run;
    }
  }
}
