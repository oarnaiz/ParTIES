package PARTIES::Utils;
use strict;


use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;
use Statistics::R;
use File::Basename;
use Bio::SeqIO;
  


#############################################
# PUBLIC FUNCTIONS
#############################################


  
=head2 read_gff_file

 Title   : Read gff file 
 Usage   : Read gff file and store the objects by seq_ids
 Function: read_gff_file($gff_file)
 Returns : Hash of arrays 
 Args    : GFF file

=cut

sub read_gff_file {
   my ($self,$gff_file) = @_;
   
   my %features;
   open(FILE,$gff_file) or die "No gff file : $gff_file";
   while(<FILE>) {
      chomp;
      next if /^#/;
      my $feat = gff3_parse_feature($_);
      if(!$feat->{attributes}->{Parent}) {
         push @{$features{$feat->{seq_id}}},$feat;
      }
   }
   close FILE;
   return \%features;

}
  
=head2 read_gff_file_by_id

 Title   : Read gff file 
 Usage   : Read gff file and store the objects by ID
 Function: read_gff_file_by_id($gff_file)
 Returns : Hash of arrays 
 Args    : GFF file

=cut

sub read_gff_file_by_id {
   my ($self,$gff_file) = @_;
   
   my %features;
   open(FILE,$gff_file) or die "No gff file : $gff_file";
   my ($id);
   while(<FILE>) {
      chomp;
      next if /^#/;
      my $feat = gff3_parse_feature($_);
      if(!$feat->{attributes}->{Parent}) {
         ($id) = @{$feat->{attributes}->{ID}};
      }
      push @{$features{$id}},$feat;
   }
   close FILE;
   return \%features;

}


  
=head2 gff_line

 Title   : write a gff line
 Usage   : gff_line($obj)
 Function: format the object in GFF3 format
 Returns : String
 Args    : Bio::GFF3::LowLevel object

=cut

sub gff_line {
   my ($self,$feat) = @_;
   my ($id) = @{$feat->{attributes}->{ID}};
   my @attributes = ("ID=$id");
   foreach my $key (sort keys %{$feat->{attributes}}) {
      next if($key eq 'ID');
      foreach my $value (@{$feat->{attributes}->{$key}}) {
         push @attributes,"$key=$value";
      }
   }
   my $score = ($feat->{score} eq '') ? '.' : $feat->{score};
   my $strand = ($feat->{strand} eq '') ? '.' : $feat->{strand};
   my $phase = ($feat->{phase} eq '') ? '.' : $feat->{phase};
   return join("\t",($feat->{seq_id},$feat->{source},$feat->{type},$feat->{start},$feat->{end},$score,$strand,$phase,join(";",@attributes)));
}

  
=head2 read_fasta_file

 Title   : Read a fasta file
 Usage   : read_fasta_file($fasta_file)
 Function: 
 Returns : Hash of Bio::Seq objects
 Args    : Fasta file

=cut

sub read_fasta_file {
   my ($self,$fasta_file) = @_;  
   my $seqio = new Bio::SeqIO(-file=> $fasta_file,-format=>'Fasta');
   my %fasta;
   while(my $seq = $seqio->next_seq) {
      my $seq_id = $seq->id;
      $fasta{$seq_id}=$seq;
   }
   return \%fasta;
}
  
=head2 read_config_file

 Title   : Read a ParTIES configuration file
 Usage   : read_config_file($cfg_file)
 Function: 
 Returns : Array
 Args    : ParTIES config text file

=cut

sub read_config_file {
   my ($self,$cfg_file) = @_;   
   open(CFG,$cfg_file) or die $cfg_file;
   my @pipeline;
   my %general;
   my $mode;
   while(<CFG>) {
      chomp;
      next if($_=~/^#/);
      if($_=~/^\[(\S+)\]/) {
         $mode = $1;
	 push @pipeline, { MODE => $mode };
	 foreach my $key (keys %general) {
	   foreach my $value (@{$general{$key}}) {
	      push @{$pipeline[$#pipeline]->{GENERAL_PARAMS}->{$key}},$value;
	   }
	 }
	 
      } elsif($_) {
         my ($key,$value) = split /=/,$_;	 
	 die "Pb syntax $_" if(!$key or $value eq '');
	 $key="-$key" if($key!~/^\-/);
	 if($mode) {	 
	    push @{$pipeline[$#pipeline]->{PARAMS}->{$key}},$value;
	    delete $pipeline[$#pipeline]->{GENERAL_PARAMS}->{$key} if($pipeline[$#pipeline]->{GENERAL_PARAMS}->{$key}); 
	 } else {
	    push @{$general{$key}},$value;
	   # $general{$key} = $value;
	 }
	 
	 
      }
   }
   close CFG;
   return @pipeline;
}
  
=head2 get_InDel_from_aln

 Title   : Get Indel from an alignments
 Usage   : get_InDel_from_aln(...)
 Function: 
 Returns : Array of Indels
 Args    : (Bio::Seq, [target alignment], [query name], [query alignment], [target start], [query start], [query length], [strand], {parameters}
 

=cut

sub get_InDel_from_aln {
   my ($self,$seq,$taln, $qname,$qaln, $begin, $qstart, $qlength, $strand, $parameters) = @_;

   my $MAC_FLANK_SEQ_LENGTH = $parameters->{JUNCTION_FLANK_SEQ_LENGTH};
   die "No mac_flank_seq_length " if(!$MAC_FLANK_SEQ_LENGTH);
   
   my $FH_ERR = $parameters->{ERROR_FILE_HANDLER} if($parameters->{ERROR_FILE_HANDLER});
   my $NOT_BOUNDED_BY_TA = ($parameters->{NOT_BOUNDED_BY_TA} eq 'TRUE') ? 1 : 0;
   
   my $CHECK_BREAK_POINTS = ($parameters->{CHECK_BREAK_POINTS}) ? $parameters->{CHECK_BREAK_POINTS} : 0;
   my $MIRAA_BREAK_POINTS = $parameters->{MIRAA_BREAK_POINTS};

   my $MIN_INDEL_LENGTH = ($parameters->{MIN_INDEL_LENGTH}) ? $parameters->{MIN_INDEL_LENGTH} : 20;
   my $MAX_MISMATCH = (defined $parameters->{MAX_MISMATCH} and $parameters->{MAX_MISMATCH} ne '') ? $parameters->{MAX_MISMATCH} : 0;
   
   my $MIN_READS_SUPPORTING_MAC_INDEL_JUNCTIONS = ($parameters->{MIN_READS_SUPPORTING_MAC_INDEL_JUNCTIONS} ne '') ? $parameters->{MIN_READS_SUPPORTING_MAC_INDEL_JUNCTIONS} : 0;
   my $SAMPLE_SAM = $parameters->{SAMPLE_BAM} if($parameters->{SAMPLE_BAM});
   my $CTL_SAM = $parameters->{CONTROL_BAM} if($parameters->{CONTROL_BAM});

   my @InDels;
   return @InDels if($taln!~/-/);

   # get length of opening gap if any:
   my $open = 0;
   if ($taln =~ /^([-]+)[ACGT]/) { $open = length($1) }
   
   die "ERROR length alignmenet taln_length=".length($taln)." qaln_length=".length($qaln)."\n" if(length($taln) != length($qaln));
   while($taln =~ /[ACGT]([-]+)[ACGT]/g) {
      my $gap = length($1);
      my $pos = pos($taln);
      my $ins = substr($qaln,$pos - $gap - 1,$gap); 
      
      # min gap length
      if($gap >= $MIN_INDEL_LENGTH) {         
	 my $nb_short_indel_in_query_before_current_gap = scalar @{[substr($qaln,0,$pos - $gap - 1) =~ /-/g]};
	 my $nb_short_indel_in_target_before_current_gap = scalar @{[substr($taln,0,$pos - $gap - 1) =~ /-/g]};
	 
	 my ($indel_start,$indel_end) = ($pos - $gap, $pos);
         my $indel_position = $begin + $indel_start - $nb_short_indel_in_target_before_current_gap;
         #print "FIND GAP pos=$indel_position => indel_start=$indel_start ",($indel_start+$qstart)," ",($indel_end+$qstart),"\n";
	 
	 my ($left_mac_flank_seq_length,$right_mac_flank_seq_length) = ($MAC_FLANK_SEQ_LENGTH,$MAC_FLANK_SEQ_LENGTH);
	 $left_mac_flank_seq_length = $indel_start - 1  if( ($indel_start-$left_mac_flank_seq_length) < 1); 
	 
	 
	    
	 my ($taln_left_flank,$taln_right_flank) = (substr($taln,$indel_start - $left_mac_flank_seq_length - 1 ,$left_mac_flank_seq_length),substr($taln,$indel_end - 1 ,$right_mac_flank_seq_length));
	 my ($qaln_left_flank,$qaln_right_flank) = (substr($qaln,$indel_start - $left_mac_flank_seq_length - 1 ,$left_mac_flank_seq_length),substr($qaln,$indel_end - 1 ,$right_mac_flank_seq_length));
	 
	 # no gap in flanks
	 if($taln_left_flank!~/[ATGCN]-+[ATGCN]/ and $taln_left_flank!~/^-+$/ and $taln_right_flank!~/[ATGCN]-+[ATGCN]/ and $taln_right_flank!~/^-+$/ and $qaln_left_flank!~/-/ and $qaln_right_flank!~/-/) {
	    
	    
	    my $max_mismatch = $MAX_MISMATCH;
	    # if flank sequence is shorter because there is close gap
	    if($taln_left_flank=~/^-+([ATGCN]+)$/) {
	       $taln_left_flank = $1;
	       $qaln_left_flank = substr($qaln_left_flank,length($qaln_left_flank)-length($taln_left_flank));
	       $max_mismatch = 0;
	    }
	    if($taln_right_flank=~/^([ATGCN]+)-+$/) {
	       $taln_right_flank = $1;
	       $qaln_right_flank = substr($qaln_right_flank,0,length($taln_right_flank));	       
	       $max_mismatch = 0;
	    }
	    
	    if($taln_right_flank!~/-/ and $taln_left_flank!~/-/
	      and $self->number_of_mismatches_between_sequences($taln_left_flank,$qaln_left_flank) <= $max_mismatch 
	      and $self->number_of_mismatches_between_sequences($taln_right_flank,$qaln_right_flank) <= $max_mismatch) {
	    
	       ($indel_position,$ins,$taln_left_flank,$taln_right_flank,$indel_start) = _adjust_InDel_position($seq,$indel_position,$ins,$taln_left_flank,$taln_right_flank,$indel_start);

	       
	       if(!$CHECK_BREAK_POINTS or _check_break_point($seq,$indel_position,$MIRAA_BREAK_POINTS,$SAMPLE_SAM)) {
	 	  my $InDel ;
		  
	 	  if ($ins =~ /^TA/ and $taln_right_flank=~/^TA/) {
	 	     ($indel_position,$ins,$indel_start) = $self->_adjust_TA_InDel_consensus($seq,$indel_position,$ins,$taln_left_flank,$taln_right_flank,$indel_start);
 
		     my $qpos = ($qstart+$indel_start-$nb_short_indel_in_query_before_current_gap);
		     $qpos = ($qlength - $qpos - length($ins)) if($strand eq '-');
		     #print "$indel_position and ",$qpos,"\n";
		     
		     $InDel = { SEQ_ID => $seq->id, POS => $indel_position, SEQ => $ins, BOUNDED_BY_TA => 'TRUE', QNAME=> $qname , QPOS => $qpos };		     
		     
	 	  } elsif($NOT_BOUNDED_BY_TA and length($taln_left_flank)==$left_mac_flank_seq_length and length($taln_right_flank)==$right_mac_flank_seq_length) {
	 	     ($indel_position,$ins,$indel_start) = _shift_left($seq,$indel_position,$ins,$taln_left_flank,$indel_start);
		     my $qpos = ($qstart+$indel_start-$nb_short_indel_in_query_before_current_gap);
		     $qpos = ($qlength - $qpos - length($ins)) if($strand eq '-');
		     $InDel = { SEQ_ID => $seq->id, POS => $indel_position, SEQ => $ins, BOUNDED_BY_TA => 'FALSE', QNAME=> $qname, QPOS => $qpos };
	 	  }
	    
	     	  if($InDel and _check_InDel_junction($seq,$indel_position,$CTL_SAM,$MIN_READS_SUPPORTING_MAC_INDEL_JUNCTIONS)) {
	 	     push @InDels,$InDel;
	          }
	       }
	    }
	 }
      }
   }
   
   return @InDels;
}
  
=head2 number_of_mismatches_between_sequences

 Title   : Number of mismatches between two sequences
 Usage   : number_of_mismatches_between_sequences($seq1,$seq2)
 Function: 
 Returns : Number of mismatch (integer)
 Args    : (String seq1, String seq2)

=cut

sub number_of_mismatches_between_sequences {
   my ($self,$seq1,$seq2) = @_;
   die "ERROR $seq1 $seq2" if(!$seq1 or !$seq2 or $seq1=~/-/ or $seq2=~/-/ or length($seq1)!=length($seq2));
   return 0 if($seq1 eq $seq2);
   my $number_of_mismatches=0;
   for(my $p=0; $p <length($seq1); $p++) {
      $number_of_mismatches++ if(substr($seq1,$p,1) ne substr($seq2,$p,1) );
   }
   return $number_of_mismatches;
}







  
=head2 significant_retention_score

 Title   : Get significancy of a retention score
 Usage   : 
 Function: 
 Returns : Hash
 Args    : (Hash $results,Hash control_rentention scores, Hash $parameters)

=cut

sub significant_retention_score {
   my ($self,$results,$control_miret,$parameters) = @_;
   
   
   my $method = $parameters->{METHOD};
   my $tab_file =$parameters->{TAB_FILE};
   my $bin_dir =$parameters->{BIN_DIR};


   my %significant;

   ## Writes the tabulated file that will be used by the R procedure
   open(TAB, ">$tab_file") or die "Can not open $tab_file";

   if($method eq 'Boundaries'){
      print TAB join("\t", qw(ID CTL_MAC CTL_LEFT CTL_RIGHT CTL_RS_LEFT CTL_RS_RIGHT CUR_MAC CUR_LEFT CUR_RIGHT CUR_RS_LEFT CUR_RS_RIGHT))."\n";	   
   } elsif($method eq 'IES'){
      print TAB join("\t", qw(ID CTL_MAC CTL_IES CTL_RS CUR_MAC CUR_IES CUR_RS))."\n";
   } else { die "Unknown method $method"; }  

   foreach my $seq_id (keys %{$results}){

      next if(!$results->{$seq_id});
      foreach my $ies (@{$results->{$seq_id}}) {
         my ($ies_id) = @{$ies->{attributes}->{ID}};
         my ($ctl_ies) = @{$control_miret->{$ies_id}};
         die "No control IES " if(!$ctl_ies);
	 die "No valid GFF file for $ies_id" if(!$ctl_ies->{attributes}->{support_mac} or !$ies->{attributes}->{support_mac});
      
         my @line = ($ies_id);
         if($method eq 'Boundaries'){
            push @line,($ctl_ies->{attributes}->{support_mac}[0],$ctl_ies->{attributes}->{support_left}[0],$ctl_ies->{attributes}->{support_right}[0],$ctl_ies->{attributes}->{retention_score_left}[0],$ctl_ies->{attributes}->{retention_score_right}[0],
	 		$ies->{attributes}->{support_mac}[0],$ies->{attributes}->{support_left}[0],$ies->{attributes}->{support_right}[0],$ies->{attributes}->{retention_score_left}[0],$ies->{attributes}->{retention_score_right}[0]);
	 
         } else {
            push @line,($ctl_ies->{attributes}->{support_mac}[0], $ctl_ies->{attributes}->{support_ies}[0],$ctl_ies->{attributes}->{retention_score}[0],
	 		$ies->{attributes}->{support_mac}[0], $ies->{attributes}->{support_ies}[0],$ies->{attributes}->{retention_score}[0]);
         } 
         print TAB join("\t", @line )."\n";
      }
   }
   close TAB;


   ## Calls the R procedures and retrieves the results
   my $R = Statistics::R->new();
#   my $path=$self->{PATH};
   $R->send('in_file="'.$tab_file.'"');
   $R->send('out_file="'.$tab_file.'"');
   $R->send('method="'.$method.'"');   
   my $source_file=$bin_dir."/utils/miret_control_comparison.R";
   $R->send('source("'.$source_file.'")');
   $R->stopR();

   ## read results 
   open(TAB, "$tab_file") or die "Can not open $tab_file";
   my $header = <TAB> ;  # header
   chomp $header;
   my @header = map { lc $_ } split("\t", $header);

   while(<TAB>) {
      chomp;
      my ($ies_id,@columns)=split("\t", $_);
      my @selected_colums = @columns[10..$#columns] ;
      for(my $i=10; $i<=$#columns; $i++) {
	 $significant{$ies_id}->{$header[$i+1]} = $columns[$i];
      }
   }
   close TAB;

   return %significant;
}

=head2 significant_retention_score

 Title   : Calculate the consensus score against the consensus 'T','A','Y','A','G'
 Usage   : consensus_scores($ies_seq)
 Function: consensus_scores
 Returns : @array (score left , score right)
 Args    : (String ies sequence)

=cut


# 
sub consensus_scores {
   my ($self,$ies) = @_;
   
   my $left_seq = substr($ies,0,5);
   my $right_seq = substr($ies,length($ies)-5,5);
   $right_seq = reverse $right_seq;
   $right_seq=~tr/ATGC/TACG/;
   my @consensus = ('T','A','Y','A','G');
   my @scores;
   foreach my $seq ( $left_seq, $right_seq) {
      my @seq = split //,$seq;
   
      my $score=0;
      for(my $i=0;$i<5;$i++) {
      
         if($consensus[$i] eq 'Y') {
            $score++ if($seq[$i] eq 'C' or $seq[$i] eq 'T');
         } else {
            $score++ if($seq[$i] eq $consensus[$i]);
         }
      }
      push @scores,($score / 5);
   }
   return @scores;  
}

  
=head2 generate_IES_id

 Title   : Generated IES identifier
 Usage   : generate_IES_id($prefix,$seq_id,$start,$end,$suffix)
 Function: 
 Returns : String IES id
 Args    : ($prefix,$seq_id,$start,$end,$suffix)

=cut

sub generate_IES_id {
   my ($self,$prefix,$seq_id,$start,$end,$suffix) = @_;
   
   my $seq_id_nb;
   if($seq_id=~/(\d+)$/) { $seq_id_nb = $1; } 
   elsif($seq_id=~/(\d+)_with_IES$/) { $seq_id_nb = $1; } 
   elsif($seq_id=~/(\d+)$suffix$/) { $seq_id_nb = $1; } 
   else { $seq_id_nb = $seq_id; }
   
   #my $id = $PREFIX."_PGM_$scaffold_nb\_$ta";
   my @id = ($prefix,$seq_id_nb,$start);
   push @id,$end if($end);
   return join(".",@id);
}


  
=head2 get_junction_seq

 Title   : Get Junction sequence
 Usage   : get_junction_seq($seq,$pos,$flank_seq_length) 
 Function: 
 Returns : String 
 Args    : ($seq,$pos,$flank_seq_length) 

=cut
sub get_junction_seq {
   my ($self,$seq,$pos,$flank_seq_length) = @_;
   die "$seq,$pos,$flank_seq_length" if(!$seq or !$pos or !$flank_seq_length);
   my $left = ($pos >= ($flank_seq_length+1)) ? lc(substr($seq->seq,$pos - ($flank_seq_length+1),$flank_seq_length)) : lc(substr($seq->seq,0,$flank_seq_length));   
   return $left. uc(substr($seq->seq,$pos-1,2)) .lc(substr($seq->seq,$pos+1,$flank_seq_length));
}





#############################################
# PRIVATE FUNCTIONS
#############################################




# Check break point consistency with mapping of reads
sub _check_break_point {
   my ($seq,$pos,$miraa_break_points,$sam) = @_;
   my $seq_id = $seq->id;
   
   my $range = 100;
   
   my ($start_region,$end_region) = ($pos-$range,$pos+$range); 
   
   if($miraa_break_points and scalar @{$miraa_break_points}) {
      foreach my $bp (@{$miraa_break_points}) {
	 return 1 if($bp->{start} <= $end_region and $bp->{end} >= $start_region);
      }
   }
   return 1 if(!$sam);   
   
   my ($full_match,$partial_match) = (0,0);
   for my $a ($sam->get_features_by_location(-seq_id => $seq_id,
                                                 -start  => $start_region,
                                                 -end    => $end_region)) {
      my $cigar_str = $a->cigar_str;
      if($cigar_str=~/^\d+M$/) {
         $full_match++;
      } else {
         $partial_match++;
      }
   }
   return 1 if($partial_match >= 5 or ($full_match+$partial_match) == 0);
   return 0;   
   
   
   
}

# Adjust Indel position (try to place the insertion between two TA)
sub _adjust_InDel_position {
   my ($seq,$indel_position,$indel_seq,$taln_left_flank,$taln_right_flank,$qstart) = @_;
   if($indel_seq =~ /^TA/ and $taln_right_flank=~/^TA/) {
      return ($indel_position,$indel_seq,$taln_left_flank,$taln_right_flank,$qstart);
   } 
   if ($indel_seq =~ /^([ACGT]*?)(TA[TANGC]+)$/) {
      my ($begin_seq,$rest_seq) = ($1,$2);
      if($begin_seq eq substr($taln_right_flank,0,length($begin_seq))) {
         $indel_position+=length($begin_seq);
	 $qstart+=length($begin_seq);
	 $indel_seq = $rest_seq.$begin_seq;
	 substr($taln_right_flank,0,length($begin_seq))='';
	 $taln_left_flank.=$begin_seq;
	 #die "$indel_seq $taln_left_flank $taln_right_flank";
	 return ($indel_position,$indel_seq,$taln_left_flank,$taln_right_flank,$qstart);
      }
   }   
   # very rare case with muscle alignement
   if ($indel_seq =~ /^([ACNGT]+)(TA[ACGT]?)$/) {
      my ($rest_seq,$end_seq) = ($1,$2);
      if($end_seq eq substr($taln_left_flank,length($taln_left_flank)-length($end_seq))) {
         $indel_position-=length($end_seq);
	 $qstart-=length($end_seq);
	 $indel_seq = $end_seq.$rest_seq;
	 substr($taln_left_flank,length($taln_left_flank)-length($end_seq))='';
	 $taln_right_flank="$end_seq$taln_right_flank";
	 return ($indel_position,$indel_seq,$taln_left_flank,$taln_right_flank,$qstart);
      }

   } 
  return ($indel_position,$indel_seq,$taln_left_flank,$taln_right_flank,$qstart);
}


# Between equivalent postions : prefer these with the better consensus sequence (IES consensus score)
sub _adjust_TA_InDel_consensus {
   my ($self,$seq,$indel_position,$indel_seq,$taln_left_flank,$taln_right_flank,$qstart) = @_;
   
   return ($indel_position,$indel_seq,$qstart) if($taln_right_flank!~/^TATA/);
   
   my ($best_ies_seq,$best_shift,$best_left_score,$best_right_score);
   my $i=0;
   while(substr($taln_right_flank,$i,2) eq 'TA') {
      my $new_ies_seq = substr($indel_seq,$i).substr($indel_seq,0,$i).'TA';      
      my ($left_score,$right_score) = $self->consensus_scores($new_ies_seq);
      if(!$best_ies_seq or ($best_left_score <= $left_score and $best_right_score <= $right_score and $new_ies_seq=~/^TA/)) {
         ($best_ies_seq,$best_shift,$best_left_score,$best_right_score) = (substr($new_ies_seq,0,length($new_ies_seq)-2),$i,$left_score,$right_score);
      }
      $i+=2;
   }
   $indel_position+=$best_shift;
   return ($indel_position,$best_ies_seq,$qstart+$best_shift);
}


# Try to shift at the most left position
sub _shift_left {
   my ($seq,$indel_position,$indel_seq,$taln_left_flank,$qstart) = @_;
   
   # shift to left the most as possible
   my ($shift_left,$shift_right) = (0,0);
   my $i = 1;
   while(substr($taln_left_flank,length($taln_left_flank)-$i,1) eq substr($indel_seq,length($indel_seq)-$i,1)) {
      $shift_left++;
      $i++;
   }
   my $new_indel_seq=substr($taln_left_flank,length($taln_left_flank)-$shift_left).substr($indel_seq,0,length($indel_seq)-$shift_left);
   $indel_position-=$shift_left;
   die $seq->id," $indel_seq VS $new_indel_seq" if($shift_left==0 and $indel_seq ne $new_indel_seq);
   return ( $indel_position,$new_indel_seq,$qstart-$shift_left);
   
}

# Check indel junction ; if there is reads which confirm the mac junction
sub _check_InDel_junction {
   my ($seq,$pos,$ctl_sam,$min_reads_supporting_mac_indel_junctions) = @_;
   
   return 1 if(!$ctl_sam or $min_reads_supporting_mac_indel_junctions==0);
   
   my ($range,$min_overlap) = (100,4);   
   my ($start_region,$end_region) = ($pos-$range,$pos+$range);
   my ($full_perfect_match) = (0);
   for my $a ($ctl_sam->get_features_by_location(-seq_id => $seq->id,
                                                 -start  => $start_region,
                                                 -end    => $end_region)) {
      next if($a->end <= $pos or $a->start >=$pos);
      
      
      if(abs($pos-$a->end) >= $min_overlap and abs($pos-$a->start) >= $min_overlap ){
         my $mismatch = $a->get_tag_values('NM');
      
         my $cigar_str = $a->cigar_str;
         if($cigar_str=~/^\d+M$/ and $mismatch ==0) {
            $full_perfect_match++;
         } 
      }
   }
   return 1 if($full_perfect_match >= $min_reads_supporting_mac_indel_junctions); 
   return 0;
}







1;
