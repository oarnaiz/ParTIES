package PARTIES::MILORD;
use strict;
use base 'PARTIES::Root';

use PARTIES::Config;
use File::Basename;
use Data::Dumper;
use List::Util qw(max min);

=head1 NAME

 ParTIES MILORD module - Method of Identification and Localization of Rare Deletion

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Mapping file of reads on a genome"
				},
			MIN_SIZE => {
				MANDATORY=>1, DEFAULT=>5, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Minimum size for a deletion to be reported"
				},
			MAX_SIZE => {
				MANDATORY=>1, DEFAULT=>1e3, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Maximum size for intra-chromosomic deletion to be reported (Inf for no restriction)"
				},
			JUNCTION_FLANK_SEQ_LENGTH => {
				MANDATORY=>1, DEFAULT=>15, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Length of the flanking sequence to report, around the deletion"
				},
			NOT_BOUNDED_BY_TA => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 2,
				DESCRIPTION=>"Should deletion not bounded by TA dinucleotide be reported"
				},	
			PREFIX => {
				MANDATORY=>0, DEFAULT=>'MILORD', TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Prefix of the ID of detected deletions"
				},	
			USE_INSERT_REFERENCE => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Try to use the result of the Insert module as the reference genome"
				},	
			REPORT_READS => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Should read names be reported when showing a deletion"
				},	
			CONSIDER_TRANS_DELETION => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 2,
				DESCRIPTION=>"Consider deletion between sequences (trans-deletions) or superior to -max_size"
				},			
			CONSIDER_OVERLAPPING => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 2,
				DESCRIPTION=>"Consider overlapping deletions"
				}
		);



my %FILE_EXTENSIONS = ( 
			gff3 => {DESC=> 'GFF3 file', EXT => 'gff3' },	
			stats => {DESC=> 'Stats file', EXT => 'stats' },
			'R1.sam' => {DESC=> 'SAM file Part 1', EXT => 'R1.sam' },	
			'R2.sam' => {DESC=> 'SAM file Part 2', EXT => 'R2.sam' },		
			 );

=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::MILORD->new({....})
 Function: Create a factory object
 Returns : the factory object

=cut
 
sub new {
  my ($class,$ref) = @_;
  
  # Check the class
  $class = ref($class) || $class;

  # Link between the object and the class
  my $self = bless {},$class;
  
  # init parent class
  $self = $self->SUPER::new($ref);
  
  # load parameters
  foreach my $param (sort keys %PARAMETERS) {
     $self->{$param} = ((defined $ref->{"-".lc($param)} and $ref->{"-".lc($param)} ne '') and (ref($ref->{"-".lc($param)}) ne 'ARRAY' or scalar @{$ref->{"-".lc($param)}} !=0)) ? $ref->{"-".lc($param)} : $PARAMETERS{$param}->{DEFAULT}; 
	die "$param should be TRUE or FALSE not : $self->{$param}\n" if($PARAMETERS{$param}->{TYPE} eq 'BOOLEAN' and $self->{$param} ne 'TRUE' and $self->{$param} ne 'FALSE');
  } 
  return $self;
}

#############################################
# PRIVATE FUNCTIONS
#############################################

sub get_parameters {
  my ($self) = @_;  
  return \%PARAMETERS;
}

=head2 _check_mandatory_parameters

 Title   : "Private" function _check_mandatory_parameters
 Usage   : $factory->_check_mandatory_parameters
 Function: Check the mandatory parameters or deduce them with -auto
 Returns : 1 (success) or 0 (error)

=cut

sub _check_mandatory_parameters {
  my ($self) = @_;  

  # IF MODE AUTO
  if($self->auto) {  
     
     if($self->{USE_INSERT_REFERENCE} eq 'TRUE') {
        $self->{GENOME} = $self->{OUT_DIR}."/Insert/Insert.fa" if(-e $self->{OUT_DIR}."/Insert/Insert.fa");
	
        my $out_bam = $self->{OUT_DIR}."/Map/".basename($self->{OUT_DIR}).".Insert.fa.BOWTIE.sorted.bam";
        $self->{BAM} = $out_bam if(-e $out_bam);
     }
	  else{
	     my $bam = $self->{OUT_DIR}."/Map/".basename($self->{OUT_DIR}).".".basename($self->{GENOME}).".BOWTIE.sorted.bam";
	     $self->{BAM} = $bam if(!$self->{BAM} and -e $bam);
	  }
  }
  
  return $self->SUPER::_check_mandatory_parameters;
}


sub get_file_extensions {
  my ($self) = @_;  
  return \%FILE_EXTENSIONS;
}


#############################################
# PUBLIC FUNCTIONS
#############################################



=head2 init

 Usage   : $factory->init
 Function: Initiate processes before multi thread calculation
 Returns : Nothing
 Args    : Nothing

=cut

sub init {
  my ($self) = @_;  
  $self->SUPER::_init;  
  
  
  # load GENOME file
  #if( $self->{CONSIDER_TRANS_DELETION} eq 'TRUE') {
     $self->stderr("Read ".basename($self->{GENOME})." ... " );
     $self->{LOADED_GENOME} = PARTIES::Utils->read_fasta_file($self->{GENOME});
     $self->stderr("Done\n" );
  #}

}


=head2 calculate

 Title   : Calculate function
 Usage   : $factory->calculate($seq);
 Function: 
 Returns : Nothing
 Args    : Bio::Seq object

=cut


sub calculate {
   my ($self,$seq) = @_;  
   my $seq_id = $seq->id;

   $self->stderr("Calculation $seq_id\n");
  
   my $sam = new Bio::DB::Sam( -bam => $self->{BAM}, -fasta => $self->{GENOME}, -expand_flags => 'true');
   
   
   my %stats;
   
   ########
   # Step 1 : Gathering partially mapped reads
   #			=> Parse the BAM file to detect partially mapped reads
   #			=> Hash with partially mapped reads
   #			=> File with unmapped sequences
   #	$self->stderr("Gathering partially mapped reads on $seq_id\n");  
   #$self->stderr("_get_partially_mapped $seq_id\n");
   my %reads = $self->_get_partially_mapped($sam,$seq_id,\%stats);

   
   my %segments;
   
   if(scalar keys %reads) {

       ########
       # Step 2 : Retrieving unmapped sequence localization
       #			=> Create single scaffold reference
       #			=> Index the reference (BOWTIE2)
       #			=> Map the sequences on the new ref (BOWTIE2)
       #	$self->stderr("Remapping reads on $seq_id\n");   
       #$self->stderr("_align_reads $seq_id\n");
       my ($bam_file_for_partially_mapped_reads,$target_reference) = $self->_align_reads(\%reads,$seq);

       ########
       # Step 3 : Reconstituting the read structure
       #			=> Parse the BOWTIE2 mapping
       #			=> Check coherence between initially mapped read part and newly mapped part
       #			=> Determine the theoric difference between the read and the reference
       #	$self->stderr("Read reconstitution on $seq_id\n'");
       #$self->stderr("_get_deletions $seq_id\n");
       my @deletions = $self->_get_deletions(\%reads,$seq,$bam_file_for_partially_mapped_reads,$target_reference,\%stats);
       
       ########
       # Step 4 : Simplifiing Complexity
       #			=> Parse the Hash of reads to look for similar coordinates
       #			=> Merge similar observations
       #	$self->stderr("Building segments from $seq_id\n");
       #$self->stderr("_simplify $seq_id\n");
       %segments = $self->_simplify(\@deletions,$seq_id,$sam);
   }



   if(%segments) {
       $self->stderr("Write GFF $seq_id\n" );
       ## Writes the final GFF file
       my $gff_file = $self->{PATH}."/tmp/".$seq_id.".gff3";

       open(GFF,">$gff_file") or die "Can not open $gff_file";
       foreach my $segment_id (sort {$segments{$a}->{start}<=>$segments{$b}->{start}} keys %segments) {
          my ($start,$end,$strand);
          my @attributes= ("ID=$segment_id","Name=$segment_id");
          foreach my $key (qw(bounded_by_ta size deletion_type is_coherent is_overlapping support_variant support_ref)) {
             push @attributes,"$key=".$segments{$segment_id}->{$key} if($key and $segments{$segment_id}->{$key});
          }
          
          
          if(length($segments{$segment_id}->{sequence})  > $self->{MAX_SIZE}) {
              my $trim_seq= substr($segments{$segment_id}->{sequence},0,$self->{MIN_SIZE}).'XXX'.substr($segments{$segment_id}->{sequence},length($segments{$segment_id}->{sequence})-$self->{MIN_SIZE});
              push @attributes,"sequence=$trim_seq";
          } else {
              push @attributes,"sequence=".$segments{$segment_id}->{sequence};
          }
          
          
          if($segments{$segment_id}->{is_coherent} eq 'TRUE') {
             ($start,$end,$strand) = ($segments{$segment_id}->{start}, $segments{$segment_id}->{end},'.');
          } else {
             ($start,$end) = ($segments{$segment_id}->{start}, $segments{$segment_id}->{start});
             foreach my $key (qw(target_seq_id)) {
                push @attributes,"$key=".$segments{$segment_id}->{$key} if($key and $segments{$segment_id}->{$key});
             }
         
             $strand = ($segments{$segment_id}->{aln_type} eq 'MS') ? '+' : '-';
             my $target_strand = ( ($segments{$segment_id}->{aln_type} eq 'MS' and $segments{$segment_id}->{strand}->[1] < 0) 
                    or ($segments{$segment_id}->{aln_type} eq 'SM' and $segments{$segment_id}->{strand}->[1] > 0)) 
                    ? '+' :'-';
         
         
         
             push @attributes,"target_pos=".$segments{$segment_id}->{end},"target_strand=$target_strand";
          
          }
          if($self->{REPORT_READS} eq 'TRUE') {
              my @read_names = @{$segments{$segment_id}->{read_names}};
              @read_names = ($read_names[0],$read_names[1],$read_names[2],$read_names[3],'...') if(scalar @read_names > 4);
              #die @read_names if(scalar @{$segments{$segment_id}->{read_names}} > 4);
             push @attributes,"read_names=".join(',',@read_names);
          }
          print GFF join("\t", $segments{$segment_id}->{seq_id}, 'MILORD', 'internal_eliminated_sequence',
                      $start, 
                      $end,
                      '.', $strand, '.', join(";",@attributes)),"\n";      
       
       

       }
       close(GFF); 
   }
   my %results;
   $results{$self->get_mode()}->{$seq_id} = \%segments;
   
   
   my $stats_file = $self->{PATH}."/tmp/$seq_id.stats";
   open(STATS,">$stats_file") or die "Can not open $stats_file";
   my @stats_line = ($seq_id);
   my @stats_columns = qw(NB_READs NB_READ_TO_REMAP NB_READ_REMAPPED NB_READ_CONSIDERED NB_READ_SHOWING_DELETION NB_READ_REMAPPED_but_mismatch NB_READ_REMAPPED_but_not_full_match NB_READ_REMAPPED_but_close_to_borders NB_READ_REMAPPED_but_overlapping NB_READ_REMAPPED_but_too_short_or_long);
   print STATS "#",join("\t",'SEQ_ID',@stats_columns),"\n";   
   foreach my $key (@stats_columns){
      push @stats_line, ($stats{$key}) ? $stats{$key} : 0;
   }
   print STATS join("\t",@stats_line),"\n";
   close STATS;
   
   
   system("rm -f ".$self->{PATH}."/tmp/$seq_id*.fa* ".$self->{PATH}."/tmp/*$seq_id*.bam*");
      
   $self->stderr("End calculation $seq_id\n" );
   return %results;
}

=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
   my ($self, $results)=@_;
   my $mode = $self->get_mode();

   $self->SUPER::finish;
   
   if($self->{REPORT_READS} eq 'TRUE') {
      $self->stderr("Report reads\n");
      my $samtools = PARTIES::Config->get_program_path('samtools');
      my $ref = $self->{GENOME};
      foreach my $part (qw(R1 R2)) {
         my $base_sam = $self->{PATH}."/$mode.$part";
         #system("$samtools view -Sb -T $ref $base_sam.sam 2> /dev/null | samtools sort - $base_sam.sorted > /dev/null 2>&1 && $samtools index $base_sam.sorted.bam && rm $base_sam.sam");
         system("$samtools view -Sb -T $ref $base_sam.sam > $base_sam.bam 2> /dev/null && rm $base_sam.sam");
      }
   }
   
}



sub _simplify {
   my ($self,$deletions,$seq_id,$sam) = @_;
   
   my %segments;
   if($self->{REPORT_READS} eq 'TRUE') {
      open(FIRSTSAM,">". $self->{PATH}."/tmp/".$seq_id.".R1.sam") or die $self->{PATH}."/tmp/".$seq_id.".R1.sam";
      open(SCDSAM,">". $self->{PATH}."/tmp/".$seq_id.".R2.sam") or die $self->{PATH}."/tmp/".$seq_id.".R2.sam";
      
   }
   
   my ($samtools) = (PARTIES::Config->get_program_path('samtools'));
   
   foreach my $deletion (@$deletions) {
      next if($self->{CONSIDER_TRANS_DELETION} eq 'FALSE' and $deletion->{deletion_type} eq 'INTER_CHR');
      my $segment_id;
      if($deletion->{is_coherent} eq 'TRUE') {
         $segment_id=PARTIES::Utils->generate_IES_id($self->{PREFIX},$deletion->{seq_id},$deletion->{start},$deletion->{end}); 	      
      } else {
         my $first_part_strand = ($deletion->{aln_type} eq 'MS') ? '+' : '-';
         my $second_part_strand = ( ($deletion->{aln_type} eq 'MS' and $deletion->{strand}->[1] < 0) or ($deletion->{aln_type} eq 'SM' and $deletion->{strand}->[1] > 0)) ? '+' :'-';
      
         $segment_id=PARTIES::Utils->generate_IES_id($self->{PREFIX},$deletion->{seq_id},join('.',$deletion->{start},$first_part_strand)
	 						,join('.',PARTIES::Utils->get_seq_id_number($deletion->{target_seq_id}),$deletion->{end},$second_part_strand));
      }
      
      
      if(!$segments{$segment_id}) {
         $segments{$segment_id} = $deletion;
         $segments{$segment_id}->{support_variant} = 0;
         $segments{$segment_id}->{support_ref} = ($deletion->{is_coherent} eq 'TRUE') ? $self->_get_reference_support($sam, $deletion->{seq_id}, $deletion->{start}) : 'NA';
      } 
      $segments{$segment_id}->{support_variant}++;
      push @{$segments{$segment_id}->{read_names}}, $deletion->{read_name};
      
      
      if($self->{REPORT_READS} eq 'TRUE') {
         my ($fisrt_part_sam_string,$second_part_sam_string) =qw(NA NA);
         if($deletion->{aln_type} eq 'MDM') {
	    my $cmd = join(" ",$samtools,'view',$self->{BAM},join("",$seq_id,':',$deletion->{aln_start},'-',$deletion->{aln_end}),' | grep',$deletion->{read_name});
	    my $sam_line = `$cmd`;
	    chomp $sam_line;
	    my @sam_line = split /\t/,$sam_line;
	    if(@sam_line and $sam_line[5] =~/^(\d+)M(\d+)D(\d+)M$/) {
	       my ($match1,$del,$match2) = ($1,$2,$3);
	       		
	       my @atts = qw(AS:i:0 XN:i:0 XM:i:0 XO:i:0 XG:i:0 NM:i:0 YT:Z:UU);
	       my $flag = ($deletion->{strand}->[0] <0) ? 0 : 16;
	       $fisrt_part_sam_string=join("\t",$sam_line[0],$flag,$sam_line[2],$sam_line[3],$sam_line[4],$match1.'M',qw(* 0 0),substr($sam_line[9],0,$match1),substr($sam_line[10],0,$match1),@atts,'MD:Z:'.$match1);
	       $second_part_sam_string=join("\t",$sam_line[0],$flag,$sam_line[2],$sam_line[3]+$match1,$sam_line[4],$match2.'M',qw(* 0 0),substr($sam_line[9],$match1),substr($sam_line[10],$match1),@atts,'MD:Z:'.$match2);
	     
	    }
	 
	 } else {
            ($fisrt_part_sam_string,$second_part_sam_string) = ($deletion->{fisrt_part_sam_string},$deletion->{second_part_sam_string});
	 }
	 if($fisrt_part_sam_string and $fisrt_part_sam_string ne 'NA' and $second_part_sam_string and $second_part_sam_string ne 'NA') {
            print FIRSTSAM join("\t",$fisrt_part_sam_string,'CO:Z:'.$segment_id),"\n";
            print SCDSAM join("\t",$second_part_sam_string,'CO:Z:'.$segment_id),"\n";
	 	 
	 }

      }       
      
   }   
   if($self->{REPORT_READS} eq 'TRUE') {
      close(FIRSTSAM);
      close(SCDSAM);
      
   }


   return %segments;
}

sub _get_deletions {
   my ($self,$reads,$seq,$bam_file,$target_reference_file,$stats) = @_;
   my $seq_id = $seq->id;
   my @deletions;
   
   if(-s $bam_file and $reads and scalar keys %$reads){
      # Parsing the mapping.
      my $bam = Bio::DB::Sam->new(-bam  =>$bam_file, -fasta=>$target_reference_file)  ;
      my @target_ids = ($self->{CONSIDER_TRANS_DELETION} eq 'FALSE') ? ($seq_id) : $bam->seq_ids;
      

      my %remapped_reads;
      foreach my $target_id (@target_ids) {
         foreach my $aln ($bam->get_features_by_location(-seq_id => $target_id)){
              my $rname=$aln->query->name;
	      if($aln->get_tag_values("NM") != 0) {
	         $stats->{NB_READ_REMAPPED_but_mismatch}++;
	         next;
	      }
	      #next if($aln->strand != 1 and $self->{CONSIDER_TRANS_DELETION} eq 'FALSE');
	      if($aln->cigar_str !~/^\d+M$/) {
	         $stats->{NB_READ_REMAPPED_but_not_full_match}++;
	         next;
	      }
	      
	      # multiple mappe
	      $reads->{$rname}->{skip} = 1 if($remapped_reads{$rname});
	      
	      $remapped_reads{$rname} = $aln;
         }
      }
      
      my ($samtools) = (PARTIES::Config->get_program_path('samtools'));
      my $nb_reads_remapped = `$samtools view -F 4 $bam_file | awk '{ print \$1 }' | sort -u | wc -l`;
      chomp $nb_reads_remapped;
      $stats->{NB_READ_REMAPPED} = $nb_reads_remapped;
      
      my %mapped_sam_strings;
      my %unmapped_sam_strings;
      if($self->{REPORT_READS} eq 'TRUE') {   

         foreach my $line (`$samtools view $bam_file`) {
	    chomp $line;
	    my ($rname) = split /\t/,$line;
	    next if(!$reads->{$rname});
	    push @{$unmapped_sam_strings{$rname}}, $line;
	 }
	 my $mapped_file_bam = $self->{PATH}."/tmp/mapped_$seq_id\_reads_vs_$seq_id.BOWTIE.sorted.bam"; 
	 foreach my $line (`$samtools view $mapped_file_bam`) {
	    chomp $line;
	    my ($rname) = split /\t/,$line;
	    next if(!$reads->{$rname});
	    push @{$mapped_sam_strings{$rname}}, $line;
	 }
	 
      }
      
      
      
      foreach my $rname (keys %$reads) {
         if($reads->{$rname}->{skip}) {
	    $stats->{NB_READ_REMAPPED_but_multiple_match}++;
	    delete $reads->{$rname};
	    next;
	 }
	 next if(!$remapped_reads{$rname} and $reads->{$rname}->{aln_type} ne 'MDM');
	 #next if($reads->{$rname}->{aln_type} ne 'MDM');
	 
	 
	 $stats->{NB_READ_CONSIDERED}++;
	 
	 my ($deletion_start,$deletion_end,$deletion_size,$deletion_seq,$deletion_type,$is_coherent,$is_overlapping,$ref_seq,$target_seq_id);
	 my ($aln_start,$aln_end);
	 my @aln_strands;
	 my ($del_char);
	 

         
	 
	 # MDM
	 ###################
	 if($reads->{$rname}->{aln_type} eq 'MDM') {
	    $deletion_size = $reads->{$rname}->{aln_unmatch};
	    $deletion_start = $reads->{$rname}->{aln}->start + $reads->{$rname}->{aln_match} ;
	    $deletion_end = $deletion_start + $deletion_size - 1; 
	    ($aln_start,$aln_end) = ($reads->{$rname}->{aln}->start,$reads->{$rname}->{aln}->end);
	    $deletion_seq = $reads->{$rname}->{unmapped_seq};
	    $deletion_type = 'INTRA_CHR';
	    $is_coherent = 'TRUE';
	    $is_overlapping = 'FALSE';
	    $ref_seq = $self->get_genome_sequence($seq_id,$aln_start,$aln_end);
	    $target_seq_id = $seq_id;
	    @aln_strands = ($reads->{$rname}->{aln}->strand,$reads->{$rname}->{aln}->strand);
	    $del_char = '-' x $deletion_size;
	    
	 # MS or SM
	 ###################
	 } else {
	    
	    # check if the deletion is coherent on the same scaffold
	    $is_coherent = is_coherent($reads->{$rname}->{aln_type},$reads->{$rname},$remapped_reads{$rname});
        next if($is_coherent eq 'NA');
	    next if($self->{CONSIDER_TRANS_DELETION} eq 'FALSE' and $is_coherent eq 'FALSE');
	    $target_seq_id = $remapped_reads{$rname}->seq_id;

	    if($reads->{$rname}->{aln}->start < $self->{JUNCTION_FLANK_SEQ_LENGTH} 
	       or $remapped_reads{$rname}->start < $self->{JUNCTION_FLANK_SEQ_LENGTH}
	       or $reads->{$rname}->{aln}->end > ($self->{LOADED_GENOME}->{$seq_id}->length - $self->{JUNCTION_FLANK_SEQ_LENGTH}) 
	       or $remapped_reads{$rname}->end > ($self->{LOADED_GENOME}->{$target_seq_id}->length -  $self->{JUNCTION_FLANK_SEQ_LENGTH})) {
	       $stats->{NB_READ_REMAPPED_but_close_to_borders}++;
	       next;
	    }

	    $deletion_type = ($seq_id eq $target_seq_id) ? 'INTRA_CHR' : 'INTER_CHR';
	    
	    $is_overlapping = _overlap($seq_id,$reads->{$rname}->{aln}->start,$reads->{$rname}->{aln}->end,
   				 			$target_seq_id,$remapped_reads{$rname}->start,$remapped_reads{$rname}->end );

       # my @ordered_deletion_positions = sort {$a<=>$b} ($reads->{$rname}->{aln}->start,$reads->{$rname}->{aln}->end,$remapped_reads{$rname}->start,$remapped_reads{$rname}->end);
        
        #my $coherent_but_too_faraway = ($is_coherent eq 'TRUE' and $self->{MAX_SIZE} ne 'Inf' and ($ordered_deletion_positions[$#ordered_deletion_positions] - $ordered_deletion_positions[0]) > $self->{MAX_SIZE}) ? 'TRUE' : 'FALSE';
        #$is_coherent = 'FALSE' if($coherent_but_too_faraway eq 'TRUE');
        


	    #  IS COHERENT
	    if($is_coherent eq 'TRUE') {
	       if($reads->{$rname}->{aln_type} eq 'MS') {	 
	          ($aln_start,$aln_end) = ($reads->{$rname}->{aln}->start,$remapped_reads{$rname}->end);   
	          $deletion_start = $reads->{$rname}->{aln}->end + 1;
	          $deletion_end = $remapped_reads{$rname}->start - 1;

	       
	    
	       } elsif($reads->{$rname}->{aln_type} eq 'SM') {
	          ($aln_start,$aln_end) = ($remapped_reads{$rname}->start,$reads->{$rname}->{aln}->end);
	          $deletion_end = $reads->{$rname}->{aln}->start - 1;
	          $deletion_start = $remapped_reads{$rname}->end + 1;
	       
	       }	    
	    
	    
	       $ref_seq = $self->get_genome_sequence($seq_id,$aln_start,$aln_end);
	       $deletion_seq = $self->get_genome_sequence($seq_id,$deletion_start,$deletion_end);
	       $deletion_size = length($deletion_seq);
	       $del_char = '-' x $deletion_size;
	     #die join("\n",$ref_seq,$deletion_seq) if($rname eq 'NS500446:494:HCHMCAFXY:4:11612:23002:5121');
	   
	    #  NOT COHERENT
	    } else {
	       
	       
	       @aln_strands = ($reads->{$rname}->{aln}->strand,$remapped_reads{$rname}->strand);
	       my $ref_seq_first_part = $self->get_genome_sequence($seq_id,$reads->{$rname}->{aln}->start ,$reads->{$rname}->{aln}->end);
	       my $ref_seq_second_part = $self->get_genome_sequence($target_seq_id,$remapped_reads{$rname}->start ,$remapped_reads{$rname}->end );
	       $ref_seq_second_part = PARTIES::Utils->revcomp($ref_seq_second_part) if($remapped_reads{$rname}->strand < 0);
	       $deletion_size = 'NA';
	       
	       if($reads->{$rname}->{aln_type} eq 'MS') {
              $aln_start = $reads->{$rname}->{aln}->start;
              $deletion_start = $reads->{$rname}->{aln}->end +1;
		  
              my $deletion_first_part = $self->get_genome_sequence($seq_id,$deletion_start,$deletion_start+$self->{MIN_SIZE});
              #$deletion_first_part = PARTIES::Utils->revcomp($deletion_first_part) if($reads->{$rname}->{aln}->strand < 0);
               
              my ($deletion_second_part);
              if($remapped_reads{$rname}->strand > 0) {
                     $aln_end = $remapped_reads{$rname}->end;
                     $deletion_end =$remapped_reads{$rname}->start - 1;
                 $deletion_second_part = $self->get_genome_sequence($target_seq_id,$deletion_end - $self->{MIN_SIZE},$deletion_end);
              
              } else {
                 $aln_end = $remapped_reads{$rname}->start;
                     $deletion_end =$remapped_reads{$rname}->end+1;
                     $deletion_second_part = PARTIES::Utils->revcomp($self->get_genome_sequence($target_seq_id,$deletion_end,$deletion_end+ $self->{MIN_SIZE}));

              
              }
              $deletion_seq = join("",$deletion_first_part,'XXX',$deletion_second_part);
              $ref_seq = join("",$ref_seq_first_part,$deletion_seq,$ref_seq_second_part);
              $del_char = '-' x length($deletion_seq);
              
	       
          } elsif($reads->{$rname}->{aln_type} eq 'SM') {
              $aln_start = $reads->{$rname}->{aln}->end;
              $deletion_start = $reads->{$rname}->{aln}->start -1;
              
              my $deletion_first_part = $self->get_genome_sequence($seq_id,$deletion_start - $self->{MIN_SIZE},$deletion_start);
              #$deletion_first_part = PARTIES::Utils->revcomp($deletion_first_part) if($reads->{$rname}->{aln}->strand < 0);
              
              my ($deletion_second_part);
              if($remapped_reads{$rname}->strand > 0) {		  
                 $aln_end = $remapped_reads{$rname}->start;
                     $deletion_end =$remapped_reads{$rname}->end+1;
                     $deletion_second_part = $self->get_genome_sequence($target_seq_id,$deletion_end,$deletion_end+ $self->{MIN_SIZE});
                 
              
              } else {

                     $aln_end = $remapped_reads{$rname}->end;
                     $deletion_end =$remapped_reads{$rname}->start - 1;
                 $deletion_second_part = PARTIES::Utils->revcomp($self->get_genome_sequence($target_seq_id,$deletion_end - $self->{MIN_SIZE},$deletion_end));
              
              }
              $deletion_seq = join("",$deletion_second_part,'XXX',$deletion_first_part);
              $ref_seq = join("",$ref_seq_second_part,$deletion_seq,$ref_seq_first_part);
              $del_char = '-' x length($deletion_seq);
	       
	       
	       
	       }
	       

	       
	    
	    } 
	    
	    
	    die $rname," ",$reads->{$rname}->{aln_type}," ",$is_coherent,"\n",$reads->{$rname}->{read_aln} if(substr($ref_seq,0,5) ne substr($reads->{$rname}->{read_aln},0,5) or substr($ref_seq,length($ref_seq)-5) ne substr($reads->{$rname}->{read_aln},length($reads->{$rname}->{read_aln})-5));
            die if(!$ref_seq or !$deletion_seq);
	 

	 }
	 
	 if(length($deletion_seq)  < $self->{MIN_SIZE} or ($self->{CONSIDER_TRANS_DELETION} eq 'FALSE' and $self->{MAX_SIZE} ne 'Inf' and length($deletion_seq) > $self->{MAX_SIZE})) {
	    $stats->{NB_READ_REMAPPED_but_too_short_or_long}++;
	    next ;
	 }
	 
	 if($is_overlapping eq 'TRUE' and $self->{CONSIDER_OVERLAPPING} eq 'FALSE') {
	    $stats->{NB_READ_REMAPPED_but_overlapping}++;
	    next ;
	 
	 }

          
	 $reads->{$rname}->{read_aln}=~s/XXX/$del_char/;
	 
	 #my ($shifted_read_aln,$shift_left) = $self->_shift_left($reads->{$rname}->{read_aln},$deletion_seq);
	 

	 
#	 if($reads->{$rname}->{aln_type} ne 'MDM') {
#	 print STDERR join(" ",'>', $rname,$reads->{$rname}->{aln_type},@aln_strands),"\n";
#	 print STDERR join(" ",$seq_id,$reads->{$rname}->{aln}->start,$reads->{$rname}->{aln}->end,$target_seq_id,$remapped_reads{$rname}->start,$remapped_reads{$rname}->end),"\n";
#	 
#	 print STDERR join(" ",'G',$ref_seq,$seq_id,"aln=$aln_start..$aln_end deletion_type=$deletion_type is_coherent=$is_coherent deletion=$deletion_start..$deletion_end"),"\n";
#	 print STDERR join(" ",'R',$reads->{$rname}->{read_aln}),"\n";
#	 #print STDERR join(" ",'R',$shifted_read_aln),"\n";
#	 print STDERR join(" ",'D',$deletion_seq,$deletion_start,$deletion_end,$deletion_size),"\n";	
#	 die if($is_coherent eq 'FALSE' and $seq_id eq $target_seq_id and $reads->{$rname}->{aln_type} eq 'SM' and $aln_strands[0]>0 and $aln_strands[1]>0); 
#         }
#
	 #my $shift_left = $self->_shift_left($reads->{$rname}->{read_aln},$deletion_seq);
	 #$reads->{$rname}->{read_aln} = $shifted_read_aln if($rname eq 'NS500446:494:HCHMCAFXY:4:21504:9913:2005');
	 my ($deletion) = PARTIES::Utils->get_InDel_from_aln($seq,$reads->{$rname}->{read_aln}, $rname,$ref_seq, 0, $deletion_start, $seq->length, '+',
   						{ 
						JUNCTION_FLANK_SEQ_LENGTH => $self->{JUNCTION_FLANK_SEQ_LENGTH},
						ERROR_FILE_HANDLER => $self->{ERR},
						NOT_BOUNDED_BY_TA =>   $self->{NOT_BOUNDED_BY_TA},
						CONSIDER_TRANS_DELETION =>   $self->{CONSIDER_TRANS_DELETION},
						CHECK_BREAK_POINTS => 0,
						MIN_INDEL_LENGTH => $self->{MIN_SIZE},
						MIRAA_BREAK_POINTS => '',
						SAMPLE_BAM => 0,
						CONTROL_BAM =>  0,
						MIN_READS_SUPPORTING_MAC_INDEL_JUNCTIONS => 1,
						IS_COHERENT => $is_coherent
						 } );
						 
         #print STDERR "shift=$shift_left\n";
         
	 if(defined($deletion) and $deletion->{POS} > $self->{JUNCTION_FLANK_SEQ_LENGTH} and $deletion->{POS} < ($seq->length - $self->{JUNCTION_FLANK_SEQ_LENGTH})){
	 
	    my ($mapped_sam_string,$unmapped_sam_string) = qw(NA NA);
	    if($self->{REPORT_READS} eq 'TRUE' and $reads->{$rname}->{aln_type} ne 'MDM') {	       
	       $unmapped_sam_string = _get_closest_sam_line($target_seq_id,$remapped_reads{$rname}->start,$unmapped_sam_strings{$rname});
	       $mapped_sam_string = _get_closest_sam_line($seq_id,$reads->{$rname}->{aln}->start,$mapped_sam_strings{$rname});
	    }
	 
	 
	 
	 
	     my $shift = $reads->{$rname}->{index} - $deletion->{POS};
             #print STDERR  $deletion->{POS}," ",$reads->{$rname}->{index}," shift=$shift\n",$deletion->{SEQ};
	     
	     if($is_coherent eq 'TRUE') {
	        $deletion_start -= $shift;
	        $deletion_end -= $shift;
	     } else {	      
	       $deletion_end = ($remapped_reads{$rname}->strand > 0) ? ($deletion_end - $shift) : ($deletion_end + $shift);
	       $deletion_start = ($reads->{$rname}->{aln_type} eq 'MS') ? ($deletion_start-$shift) : ($deletion_start+$shift)


	     
	     }
	     push @deletions, {
	 			read_name=> $rname,
				aln_type=> $reads->{$rname}->{aln_type},
				seq_id => $seq_id,
				start => $deletion_start,
				end => $deletion_end,
				target_seq_id => $target_seq_id,
				size => $deletion_size,
				sequence => $deletion->{SEQ},
				deletion_type => $deletion_type,
				is_coherent => $is_coherent,
				is_overlapping => $is_overlapping,
				aln_start=>$aln_start,
				aln_end=>$aln_end,
				aln_read => $reads->{$rname}->{read_aln},
				strand => \@aln_strands,
				bounded_by_ta => $deletion->{BOUNDED_BY_TA},
				fisrt_part_sam_string => $mapped_sam_string,
				second_part_sam_string => $unmapped_sam_string
	 	          };	     
	     


            $stats->{NB_READ_SHOWING_DELETION}++;
	

         } else {
            $stats->{NB_READ_REMAPPED_but_close_to_borders}++;
	 }
	 
      }
   }
   return @deletions;
}

sub _shift_left {
   my ($self,$read_aln,$deletion_seq) = @_;
   my $shift = 0;
   my $i=1; 
   
   if($read_aln=~/^(\w+)(\-+)(\w+)$/) {
      my ($left_part_read_aln,$gap,$right_part_read_aln) = ($1,$2,$3);
      
      my $max_shift = ($deletion_seq=~/XXX/) ? $self->{MIN_SIZE} : length($left_part_read_aln);   
      my $shifted_deletion_seq = $deletion_seq;  
      while(substr($left_part_read_aln, length($left_part_read_aln)-$i,1) eq substr($deletion_seq, length($deletion_seq)-$i,1)
				and $i<=$max_shift 
				) {
	 $shift++;
	 $i++;  
	 $shifted_deletion_seq = substr($left_part_read_aln,length($left_part_read_aln)-$shift).substr($deletion_seq,0,length($deletion_seq)-$shift);	
      }
      
      if($shift) {
         my $seq_to_shift = substr($left_part_read_aln,length($left_part_read_aln)-$shift);
         
	 my $new_read_aln = join("",substr($left_part_read_aln,0,length($left_part_read_aln)-$shift),$gap,substr($left_part_read_aln,length($left_part_read_aln)-$shift),$right_part_read_aln);
         return ($new_read_aln,$shift);
      }
   }  
   return ($read_aln,$shift);
}

sub _get_closest_sam_line {
   my ($seq_id,$pos,$sam_strings) = @_;
   
   my $min_dist='';
   my $sam_string = 'NA';
   return $sam_string if(!$sam_strings);
   foreach my $sam_line (@{$sam_strings}) {
      my @sl = split /\t/,$sam_line;
      next if($sl[2] ne $seq_id);
      if($sl[1] == 256) {
         $sl[1]= 0;
	 $sam_line =join("\t",@sl);
      }
      if($min_dist eq '' or $min_dist > abs($sl[3]-$pos)) {
     	 $min_dist = abs($sl[3]-$pos);
     	 $sam_string = $sam_line;
      }
   }
   return $sam_string;
}

sub is_coherent {
   my ($aln_type,$first_part,$second_aln_part) = @_;
   
   return 'TRUE' if($aln_type eq 'MDM');
   
   my $first_aln_part = $first_part->{aln};
   
   #my ($coherent_range_start,$coherent_range_end) = ($first_part->{coherent_range_start},$first_part->{coherent_range_end});
   
   #print STDERR join(" ",$first_aln_part->seq_id," eq ",$second_aln_part->seq_id," and ", $second_aln_part->start," >= ",$first_part->{coherent_range_start}," and ",$second_aln_part->end," <= ",$first_part->{coherent_range_end}),"\n";
   my $is_in_coherent_range = ($first_aln_part->seq_id eq $second_aln_part->seq_id and $second_aln_part->start >= $first_part->{coherent_range_start} and $second_aln_part->end <= $first_part->{coherent_range_end}) ? 'TRUE' : 'FALSE';
   
   # impssible incoherent mapping 
   return 'NA' if($first_part->{must_be_coherent} eq 'TRUE' and $is_in_coherent_range eq 'FALSE');
   
   
   # Trans deletion
   return 'FALSE' if($first_aln_part->seq_id ne $second_aln_part->seq_id);
   
   
   # strand inconsistancy
   return 'FALSE' if($second_aln_part->strand <0);
   
   
   
                #~ coherent_range_start=>$coherent_range_start,
                #~ =>$coherent_range_end,
                #~ must_be_coherent=>$must_be_coherent
   
   # overlap
   #die $first_aln_part->qname if($first_aln_part->start < $second_aln_part->end and $first_aln_part->end > $second_aln_part->start );
   return 'FALSE' if(_overlap($first_aln_part->seq_id ,$first_aln_part->start,$first_aln_part->end,
   				 $second_aln_part->seq_id, $second_aln_part->start,$second_aln_part->end ) eq 'TRUE');
   
   if($aln_type eq 'MS') {
      return 'FALSE' if($first_aln_part->end > $second_aln_part->start);
   
   } elsif($aln_type eq 'SM') {
      return 'FALSE' if($first_aln_part->start < $second_aln_part->end);
   
   } else { die $aln_type };
   
   
   
   return $is_in_coherent_range;
}

sub _overlap {
   my ($seq_id1,$start1,$end1, $seq_id2,$start2,$end2) = @_;
   return 'FALSE' if($seq_id1 ne $seq_id2);
   return 'TRUE' if($start1 <= $end2 and $end1 >= $start2 );   
   return 'FALSE';
}

sub get_genome_sequence {
   my ($self,$seq_id,$start,$end) =  @_;
   return -1 if($start > $end);
   die "$seq_id $start $end ? length=".$self->{LOADED_GENOME}->{$seq_id}->length if($start <1 or $end >$self->{LOADED_GENOME}->{$seq_id}->length);
   return uc( $self->{LOADED_GENOME}->{$seq_id}->subseq($start,$end));
}

sub _align_reads {
   my ($self,$reads,$seq) = @_;
   
   my $seq_id = $seq->id;
   
   my ($bowtie2_build,$bowtie2,$samtools) = (PARTIES::Config->get_program_path('bowtie2-build'),PARTIES::Config->get_program_path('bowtie2'),PARTIES::Config->get_program_path('samtools'));
   
   # write fasta for unmapped
   my $nb_reads_partially_mapped=0;
   my $fasta_file = $self->{PATH}."/tmp/".$seq_id."_tolocate.fa";
   open(FA,">$fasta_file") or die $fasta_file;
   foreach my $rname (keys %$reads) {
      next if($reads->{$rname}->{aln_type} eq 'MDM');
      print FA ">$rname\n".$reads->{$rname}->{unmapped_seq},"\n"; 
      $nb_reads_partially_mapped++;  
   }
   close FA;
   
   
   # mapping


   


   my $unmapped_file_bam = $self->{PATH}."/tmp/unmapped_$seq_id\_reads_vs_ref.BOWTIE.sorted"; 
   
   my $ref = $self->{GENOME};  
   if(-s $fasta_file and $nb_reads_partially_mapped!=0){
       
      system("$samtools faidx $fasta_file");
   
       my $seq_fasta_file = $self->{PATH}."/tmp/".$seq_id."_as_ref.fa";
       my $seqo = new Bio::SeqIO(-file=>">$seq_fasta_file",-format=>'fasta');
       $seqo->write_seq($seq);
       $seqo->close;
       system("$bowtie2_build $seq_fasta_file $seq_fasta_file > /dev/null 2>&1");
       system("$samtools faidx $seq_fasta_file");       
       
      
      
      if($self->{CONSIDER_TRANS_DELETION} eq 'FALSE') {
         $ref = $seq_fasta_file;
      } 
      
      die "ERROR bowtie2 index genome does not exist for ",$ref,"\n" if(!-e $ref.".1.bt2") ;
      system("$bowtie2 --threads 1 --end-to-end --quiet -f  --very-sensitive -k 2 -x ".$ref." -U $fasta_file | $samtools view -F 4 -uS - 2> /dev/null | $samtools sort -o $unmapped_file_bam.bam - > /dev/null 2>&1");
      system("$samtools index $unmapped_file_bam.bam");
      
      
      if($self->{REPORT_READS} eq 'TRUE') {
         # write fasta for mapped reads
         $fasta_file = $self->{PATH}."/tmp/".$seq_id."_mapped.fa";
         open(FA,">$fasta_file") or die $fasta_file;
         foreach my $rname (keys %$reads) {
            next if($reads->{$rname}->{aln_type} eq 'MDM');
            print FA ">$rname\n".$reads->{$rname}->{mapped_seq},"\n";   
         }
         close FA;  
         my $mapped_file_bam = $self->{PATH}."/tmp/mapped_$seq_id\_reads_vs_$seq_id.BOWTIE.sorted"; 
      
         system("$bowtie2 --threads 1 --end-to-end --quiet -f  --very-sensitive -k 2 -x ".$ref." -U $fasta_file | $samtools view -F 4 -uS - 2> /dev/null | $samtools sort -o $mapped_file_bam.bam - > /dev/null 2>&1");
         system("$samtools index $mapped_file_bam.bam");

      }
      
   } else {
      system("touch $unmapped_file_bam.bam");
   }
   
   return ("$unmapped_file_bam.bam",$ref);
}



sub _get_partially_mapped {
   my ($self,$sam,$seq_id,$stats) = @_;
   
   my %putative_pcr_dup;
   my %reads;
   foreach my $pair  ($sam->features(-type   => 'read_pair',-seq_id => $seq_id)) {
      my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
      next if(!$first_mate or !$second_mate);
      
      next if($first_mate->cigar_str !~/^\d+M$/ and $second_mate->cigar_str !~/^\d+M$/);
      next if($first_mate->seq_id ne $second_mate->seq_id);
      
      my $pair_uniquename = join("__",$first_mate->seq_id,$first_mate->start,$first_mate->end,$first_mate->query->dna,
                                    $second_mate->seq_id,$second_mate->start,$second_mate->end,$second_mate->query->dna);
       # PCR duplicates
      next if($putative_pcr_dup{$pair_uniquename});
      $putative_pcr_dup{$pair_uniquename}=1;
      
      
      my @pos = sort {$a<=>$b} ($first_mate->start,$first_mate->end,$second_mate->start,$second_mate->end);
      my ($insert_start,$insert_end) = ($pos[0],$pos[$#pos]);
      
      my $aln = ($first_mate->cigar_str =~/^\d+M$/) ? $second_mate : $first_mate;
      
      my $rname=$aln->query->name;
      
      die "Left right alignement $rname" if($first_mate->start > $second_mate->end);
      my $read_partially_mapped = ($first_mate->cigar_str =~/^\d+M$/) ? 'RIGHT' : 'LEFT';
   
      my ($aln_type, $unmatch, $match) = check_mapping($aln);
      
      next if($aln_type eq "-1" || $unmatch < $self->{MIN_SIZE});
      $stats->{NB_READs}++;
      
      # mismatches
      next if($aln->get_tag_values('XM') > 0);
      
      my $read_seq = $aln->query->dna;
      next if($read_seq=~/N/ or $aln->dna =~/N/);
      if($reads{$rname}) {
         delete $reads{$rname};
         next;
      }
      
      my ($mapped_seq,$unmapped_seq,$read_aln,$index); 
      my ($coherent_range_start,$coherent_range_end,$must_be_coherent)= qw(NA NA FALSE);
      if($aln_type eq 'MDM') {
         $unmapped_seq = substr($aln->dna,$match, $unmatch);
         $read_aln = substr($read_seq, 0, $match).'XXX'.substr($read_seq, $match);
         $index = $match;
      
      } elsif($aln_type eq 'MS'){
         $mapped_seq = substr($read_seq, 0, $match);
         $unmapped_seq = substr($read_seq, $match);
         $read_aln = $mapped_seq.'XXX'.$unmapped_seq;
         $index = $match+1;
         if($read_partially_mapped eq 'RIGHT') {
             #print STDERR "# seq_id $seq_id \n" if(!$self->{LOADED_GENOME}->{$seq_id});
             ($coherent_range_start,$coherent_range_end) = ($aln->end,$self->{LOADED_GENOME}->{$seq_id}->length);
             
         } else {
             
             ($coherent_range_start,$coherent_range_end,$must_be_coherent) = ($first_mate->start,$second_mate->end,'TRUE');
             #print STDERR join(" " ,$first_mate->seq_id,$first_mate->start,$first_mate->end,$first_mate->cigar_str,"\n",$second_mate->seq_id,$second_mate->start,$second_mate->end,$second_mate->cigar_str),"\n";
             #print STDERR "$rname $read_seq $read_partially_mapped $coherent_range_start,$coherent_range_end\n";
         }
      } elsif($aln_type eq 'SM'){
         $unmapped_seq = substr($read_seq,0, $unmatch);
         $mapped_seq = substr($read_seq, $unmatch);
         $read_aln = $unmapped_seq.'XXX'.$mapped_seq;
         $index = $unmatch+1;
         if($read_partially_mapped eq 'LEFT') {
             
             ($coherent_range_start,$coherent_range_end) = (1,$aln->end);
             
         } else {
             ($coherent_range_start,$coherent_range_end,$must_be_coherent) = ($first_mate->start,$second_mate->end,'TRUE');
             #print STDERR join(" " ,$first_mate->seq_id,$first_mate->start,$first_mate->end,$first_mate->cigar_str,"\n",$second_mate->seq_id,$second_mate->start,$second_mate->end,$second_mate->cigar_str),"\n";
             #print STDERR "$rname $read_seq $read_partially_mapped$coherent_range_start,$coherent_range_end\n";
         }
      }
      #die $rname,"\n$read_seq\n".$aln->dna,"\n$unmapped_seq" if($aln->strand <0);
      
      $reads{$rname} = { 
                seq_id => $seq_id,
                read_name => $rname,
                read_seq => $read_seq,
                read_aln => $read_aln,
                cigar => $aln->cigar_str,
                aln_type => $aln_type,
                aln_match => $match,
                aln_unmatch => $unmatch,
                aln_strand => $aln->strand,
                mapped_seq => $mapped_seq,
                unmapped_seq => $unmapped_seq,
                index => $index,
                aln => $aln,
                coherent_range_start=>$coherent_range_start,
                coherent_range_end=>$coherent_range_end,
                must_be_coherent=>$must_be_coherent
        };
        $stats->{NB_READ_TO_REMAP}++;
   }
   return %reads ;


}

# Apply some filters on the mapping to discard reads
# Returns 
#		the type of mapping {MS, SM, MDM}, 
#		the length of the unmapped part,
#		the length of the mapped part,
#		the read sequence.
#
# Return -1 if cigar strand is not MS or SM or MDM
sub check_mapping{
	my ($read)=@_;
	my $cigar=$read->cigar_str;
	
	if($cigar =~ /^\d+M$/ ){return (-1,-1);} # Perfectly mapped read
	if($read->query->dna =~ /N/i ){return (-1,-1);} # N-containing read
	#if($read->get_tag_values('XM') > 0 ){return (-1,-1);} # Mismatch
	if(!$read->paired || !defined($read->mate_seq_id)){ return (-1,-1);} # Unpaired read
	if($read->seq_id ne $read->mate_seq_id){ return (-1,-1);} # Read mapped on different seq_id
	if($read->strand == $read->mstrand){return (-1,-1);} # Strand inconsistency between pair
	
	if(	 $cigar =~ /^(\d+)M(\d+)S$/){			return('MS' ,$2,$1);}
	elsif ($cigar =~ /^(\d+)S(\d+)M$/ ){		return('SM' ,$1,$2);}
	elsif ($cigar =~ /^(\d+)M(\d+)D(\d+)M$/ ){	return('MDM',$2,$1);}
	else{ return (-1,-1);}
}
sub _get_reference_support{
	my ($self, $sam, $seq_id, $pos)=@_;
	my $range = 100;
	my ($start_region,$end_region) = ($pos-$range,$pos+$range);
	my ($full_perfect_match,$partial_match) = (0,0);
	for my $a ($sam->get_features_by_location(-seq_id => $seq_id,
                              	 -start  => $start_region,
                              	 -end    => $end_region)) {
		next if($a->end <= $pos or $a->start >=$pos);
		my $mismatch = $a->get_tag_values('NM');
		my $cigar_str = $a->cigar_str;
		if($cigar_str=~/^\d+M$/ and $mismatch ==0) {$full_perfect_match++;}
		else {$partial_match++;}
	}
	return ($full_perfect_match);
}







1;
