#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;


my ($GFF, $VERBOSE,@ATTRIBUTES);
my $DELIM =',';
GetOptions( 
	    '-gff=s' => \$GFF,
	    '-attributes=s' => \@ATTRIBUTES,
	    '-delim=s' => \$DELIM,
	    '-verbose' => \$VERBOSE,
	   );


die "$0 -gff $GFF ( -attributes {multiple, default=all attributes} -delim $DELIM -verbose ) > [tabulated_file] \n" if(!$GFF);


print STDERR "# Read $GFF ... " if($VERBOSE);
open(GFF, $GFF) or die "Cannot open $GFF\n$!\n";
my @FEATURES;
my %KEYS;
while(<GFF>){
   next if($_ =~/^#/);
   chomp;
   my $feat = gff3_parse_feature($_);
   if(!@ATTRIBUTES) {
      foreach my $key (grep { ! /^ID$/ } sort keys %{$feat->{attributes}}) {
         $KEYS{$key}=1;
      }
   }
   push @FEATURES,$feat;
}
print STDERR "Done\n" if($VERBOSE);





print STDERR "# Write tab file ... " if($VERBOSE);
@ATTRIBUTES = grep { ! /^ID$/ } sort keys %KEYS if(!@ATTRIBUTES);

# HEADER
my @HEADER = qw(ID SEQ_ID SOURCE TYPE START END SCORE STRAND PHASE);
if(@ATTRIBUTES) {
   foreach my $header (grep { ! /^ID$/ } map { uc $_ } @ATTRIBUTES) {
      $header.="_ATTR" if(grep(/^$header$/,qw(ID SEQ_ID SOURCE TYPE START END SCORE STRAND PHASE))); 
      push @HEADER,$header;      
   }
}
print join("\t",@HEADER),"\n";



foreach my $feat (@FEATURES) {
   my ($id) = @{$feat->{attributes}->{ID}};
   
   my @line = ($id);
   foreach my $key (qw(seq_id source type start end score strand phase)) {
      push @line, ($feat->{$key} eq '') ? '.' : $feat->{$key};
   }
   
   foreach my $key (@ATTRIBUTES) {
      next if($key eq 'ID');
      my @values = ($feat->{attributes}->{$key}) ? @{$feat->{attributes}->{$key}} : ('NA');
      push @line, join($DELIM,@values);      
   }
   print join("\t",@line),"\n";   
   
}
print STDERR "Done\n" if($VERBOSE);


