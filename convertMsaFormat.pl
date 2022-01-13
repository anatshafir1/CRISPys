#!/usr/local/bin/perl -w

use strict;
use FileHandle;
use Bio::AlignIO;

my ($inFile, $outFile, $inFormat, $outFormat) = @ARGV;
die "usage: convertMsaFormat.pl <inFile> <outFile> <inFormat> <outFormat>\n"
    if (!defined $outFormat);

print "inFile = '$inFile'\n";
print "outFile = '$outFile'\n";
print "inFormat = '$inFormat'\n";
print "outFormat = '$outFormat'\n";
my $in  = Bio::AlignIO->new( '-format' => $inFormat , -file => $inFile);
my $out  = Bio::AlignIO->new( '-format' => $outFormat , -file => ">$outFile");

my ($alignObj, $seqStr, $trans);
while ($alignObj = $in->next_aln()) {
	$alignObj->verbose(1);
	# Otherwise, bioperl adds sequence start/stop values, causing problems 
	# with clustal/bali_score
	$alignObj->set_displayname_flat();
	$out->write_aln($alignObj);
}


