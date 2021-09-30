#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Scalar::Util qw(looks_like_number);

# Adds submitted mutations to existing table
# M. Weigand - CDC
# 20210901

&GetOptions(	'indir=s' => \my$indir,
		'now=s' => \my$now,
		'latest=s' => \my$latest,
		'log=s' => \my$log,
		'out=s' => \my$out);

my@new = glob("$indir/[!mut_template]*csv");
open(my $fh, '>>', $out) or die "Could not open file '$out' $!";
open(my $lg, '>>', $log) or die "Could not open file '$log' $!";

foreach my$n (@new){
	open IN, "$n";
	while( my$i = <IN> ){
		unless($i =~ m/^Chrom/){
			chomp$i;
			if( &REQ($i) eq "pass"){
				print $fh "$i\n";
				print $lg basename($n,'.csv').": $i\n";
			}else{
				print "ERROR IN SUBMISSION: \'$i\' \n";
			}
		}
	}
	close IN;

}
close $fh;
close $lg;


####
sub REQ {
	my@mut = split(",", $_[0]);
	if(grep { $_ eq 'NA'} @mut[0..5]){
		return("fail");
	}elsif($mut[-2] =~ /YourName/){
		return("fail");
	}elsif($mut[-1] eq 'YYYYMMDD' ){
		return("fail");
	}elsif(!looks_like_number($mut[2])){
		return("fail");
	}elsif($mut[2] >= 29903){
		return("fail");
	}else{
		return("pass");
	}
}

