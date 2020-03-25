#!/usr/bin/env perl
#
# Purpose: to mask low-coverage regions in reference file and in FASTA.
# Inputs: VCF, BAM, reference
# Outputs: modified VCF, modified reference and/or masked consensus

use Getopt::Long;
use warnings;
use strict;
my $VERSION = "0.1";

for my $exe (qw(samtools bcftools)) {
  my($which) = qx(which $exe 2> /dev/null);
  $which or die("Can not find required '$exe' in PATH");
}

my %opts;
GetOptions(\%opts, qw(bam=s vcf=s reference=s refout=s vcfout=s consout=s depth=i qual=i help version));
usage() if ($opts{help});

if ($opts{version}) { print("vcf_mask_low_coverage v$VERSION\n"); exit;}

my %depth;
my %refhash;
$opts{refout} ||= "$opts{reference}.masked.fasta";
$opts{vcfout} ||= "$opts{vcf}.masked.vcf";

die("No BAM file at $opts{bam}") if (! -f $opts{bam});
open(my $bamd, "-|", "samtools depth -a $opts{bam}") || die "couldn't open bam $opts{bam}";
while(<$bamd>) {
	chomp;
	my @line = split(/\t/);
	$depth{$line[0]}{$line[1]} = $line[2];
}
close($bamd);

open(my $refh, "<", "$opts{reference}") || die("couldn't open reference $opts{reference}");
my $cref;
while(<$refh>) {
	chomp;
	if (/^>(.+?)(\s|$)/) {
		$cref = $1;
		$refhash{$cref} = [];
	}
	elsif ($cref) {
		my @seq = split(//);
		push(@{$refhash{$cref}}, @seq);
	}
	else {
		die("Strange reference file");
	}
}
close($refh);

for my $ref (keys %refhash) {
	for my $i (keys @{$refhash{$ref}}) {
		if (! defined($depth{$ref}{$i+1}) || $depth{$ref}{$i+1} < $opts{depth}) {
			$refhash{$ref}[$i] = "N";
		}
	}
}

my $vcfh;
my $vcfout;
if (defined($opts{qual})) {
	open( $vcfh, "-|", qq(bcftools view -i "QUAL >= $opts{qual}" $opts{vcf}) ) || die "coudn't open vcf $opts{vcf}";
} else {
	open( $vcfh, "-|", qq(bcftools view $opts{vcf}) ) || die "coudn't open vcf $opts{vcf}";
}

open($vcfout, "|-", "bcftools view -Oz - >$opts{vcfout}.gz") || die "trouble opening vcf output location $vcfout";

while(<$vcfh>) {
	if (/^#/) {
		print($vcfout $_);
	} 

	else {
		my @line = split(/\t/);
		my $restore = 0;
		# line[1] coord,  [3]-ref ; [4]-alt ; [0] - chrom
		if (length($line[3]) != length($line[4]) ) {     #need to sub letters back in to refereence for deletion to be applied
			my @sites = split(//, $line[3]);
			for my $i (keys @sites) {
				if ($depth{$line[0]}{$line[1]+$i} >= $opts{depth}) { # depth is 1-based
					$restore = 1;
				}
				if ($restore == 1) {
					$refhash{$line[0]}[$line[1]-1+$i] = $sites[$i]; # -1 since refhash is 0-based, restore the base
				}
			}
		} 
		if ($depth{$line[0]}{$line[1]} >= $opts{depth} || $restore == 1) {
			print($vcfout $_) ;
		}
			
	}
}
close($vcfout);
system("bcftools index $opts{vcfout}.gz");
open(my $refout, ">", "$opts{refout}");
for my $ref (keys %refhash) {
	printf($refout ">%s\n", $ref);
	my $line;
	for my $i (keys @{$refhash{$ref}}) {
		$line .= $refhash{$ref}[$i];
		if ($i % 60 == 59 || $i == $#{$refhash{$ref}} ) {
			printf($refout "%s\n", $line);
			$line = "";
		}
	}
}
close($refout);

if ($opts{consout}) {
	my $conscmd = qq(cat $opts{refout} | bcftools consensus $opts{vcfout}.gz > $opts{consout});
	print(STDERR "Consensus command: $conscmd\n");
	system($conscmd);
}

sub usage {
	print("$0 --bam <bamfile> --vcf <vcffile> --reference <reference.fasta>  [--depth <min-depth>] [--qual <QUAL score> [--refout <masked-ref output>] [--vcfout <masked vcf output>] [--consout <masked consensus output] [--help] [--version]\n");
	exit;
}