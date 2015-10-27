#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;
use FindBin qw($Bin);
use lib "$Bin/../lib/";
use PARTIES::Utils;

my ($VERBOSE);
my @MICA;
GetOptions( 
	    '-mica=s' => \@MICA,
	    	   
	    '-verbose' => \$VERBOSE,
	   );


die "$0 {-mica [file] } (-verbose ) > [tabulated_file] \n" if(!@MICA);

my %IES;
print "##gff-version 3\n";

foreach my $gff (@MICA) {
   print STDERR "# Read $gff ... " if($VERBOSE);
   open(GFF, $gff) or die "Cannot open $gff\n$!\n";
   while(<GFF>){
      next if($_ =~/^##/);
      if($_ =~/^#/) {
         print $_;
         next;
      }
      chomp;
      my $feat = gff3_parse_feature($_);
      push @{$feat->{attributes}->{evidences}},1;
      if(!$IES{$feat->{seq_id}}->{$feat->{start}}) {          
          $IES{$feat->{seq_id}}->{$feat->{start}} = $feat;
      } else {
         my @sequences = @{$feat->{attributes}->{sequence}};
	 my ($bounded_by_TA) = @{$feat->{attributes}->{bounded_by_ta}};
	 push @sequences,@{$feat->{attributes}->{alternative_seq}} if($feat->{attributes}->{alternative_seq});
	 
	 my @last_sequences = @{$IES{$feat->{seq_id}}->{$feat->{start}}->{attributes}->{sequence}};
	 push @last_sequences, @{$IES{$feat->{seq_id}}->{$feat->{start}}->{attributes}->{alternative_seq}} if($IES{$feat->{seq_id}}->{$feat->{start}}->{attributes}->{alternative_seq});	 
	 my ($last_bounded_by_TA) = @{$IES{$feat->{seq_id}}->{$feat->{start}}->{attributes}->{bounded_by_ta}};
	 
	 my $found_same_seq = 0;
	 foreach my $seq (@sequences) {
	    foreach my $last_seq (@last_sequences) {
	       $found_same_seq = 1 if($seq eq $last_seq);
	    }
	 }
	 
	 
	 if(!$found_same_seq) {
	    if($last_bounded_by_TA eq $bounded_by_TA or $last_bounded_by_TA eq 'TRUE') {
	       push @{$IES{$feat->{seq_id}}->{$feat->{start}}->{attributes}->{alternative_seq}},@sequences;
	       my ($evidences) = @{$feat->{attributes}->{evidences}};
	       
	    } else {
	       $IES{$feat->{seq_id}}->{$feat->{start}} = $feat;
	       push @{$IES{$feat->{seq_id}}->{$feat->{start}}->{attributes}->{alternative_seq}},@last_sequences;	       
	    }	 
	 } 
	 
	 my ($last_evidences) = @{$IES{$feat->{seq_id}}->{$feat->{start}}->{attributes}->{evidences}};
	 @{$IES{$feat->{seq_id}}->{$feat->{start}}->{attributes}->{evidences}} = ($last_evidences+1);
	 
	 	      
      }
      
   }
   close GFF;
   print STDERR "Done\n" if($VERBOSE);
}   


foreach my $seq_id (sort keys %IES) {
   foreach my $start (sort {$a<=>$b} keys %{$IES{$seq_id}}) {
      print  "",PARTIES::Utils->gff_line($IES{$seq_id}->{$start}),"\n";
   }
}

