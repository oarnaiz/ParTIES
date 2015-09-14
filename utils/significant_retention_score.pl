#!/usr/bin/perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib/";
use Getopt::Long;
use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;
use PARTIES::Utils;

my ($MIRET,$CONTROL,$TAB,$OUT);
my ($METHOD) = ('Boundaries');

GetOptions( 
	    '-miret=s' => \$MIRET,
	    '-control=s' => \$CONTROL,
	    '-out=s' => \$OUT,
	    '-method=s' => \$METHOD,
	    '-tab' => \$TAB,
	   );
die "$0 -miret $MIRET -control $CONTROL -out $OUT (-method $METHOD -tab) -out [out file] \n" if(!$MIRET or !$CONTROL or !$OUT);

die "method=$METHOD ?" if($METHOD ne 'Boundaries' and $METHOD ne 'IES');

## Load CONTROL
my $control_miret = PARTIES::Utils->read_gff_file_by_id($CONTROL);


## Load MIRET
my $miret = PARTIES::Utils->read_gff_file($MIRET);


# calculate significance
my %significant = PARTIES::Utils->significant_retention_score($miret,$control_miret,{ METHOD => $METHOD,
	   											TAB_FILE => $OUT,
												BIN_DIR => "$Bin/../",
												}
												);       


# write gff3 file
if(!$TAB) {
   
    open(GFF,">$OUT") or die "Can not open $OUT";
    foreach my $seq_id (keys %{$miret}){

       next if(!$miret->{$seq_id});
       foreach my $ies (@{$miret->{$seq_id}}) {
    	  my ($ies_id) = @{$ies->{attributes}->{ID}};
    	  if($significant{$ies_id}) {
    	     foreach my $key (keys %{$significant{$ies_id}}) {
    		# replace old value
    		$ies->{attributes}->{$key} = () if($ies->{attributes}->{$key}); 		
    		push @{$ies->{attributes}->{$key}},$significant{$ies_id}->{$key};
    	     }
    	  }
    	  print GFF PARTIES::Utils->gff_line($ies),"\n";
       }
    }
    close GFF;

}
