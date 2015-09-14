#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib/";
use Getopt::Long;


use PARTIES::Config;

my ($MODE);

GetOptions( 
	    '-mode=s' => \$MODE,
	   );
die "$0 -mode $MODE\n" if(!$MODE);


my %MODES = PARTIES::Config->get_parties_modes();

my $PACKAGE = 'PARTIES::'.$MODE;
my $factory = $PACKAGE->new({	});    


my %all_parameters=(%{$factory->get_root_parameters}, %{$factory->get_parameters});
print "[$MODE]\n";
foreach my $param (sort {$all_parameters{$b}->{MANDATORY} <=> $all_parameters{$a}->{MANDATORY} } keys %all_parameters) {
   next if($param eq 'GENOME' or $param eq 'OUT_DIR');
   my $value = $all_parameters{$param}->{DEFAULT};
   print "#" if(!$all_parameters{$param}->{MANDATORY});
   print lc($param),"=$value\n";
}
print "\n";
