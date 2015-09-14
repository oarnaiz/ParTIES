package PARTIES::MIRAA;
use strict;
use base 'PARTIES::Root';
use PARTIES::Config;

use File::Basename;

=head1 NAME

 ParTIES MIRAA module - Method of Identification by Read Alignment Anomalies

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Mapping file"
				},
			MAX_INSERT_SIZE => {
				MANDATORY=>0, DEFAULT=>1000, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Max sequencing insert size"
				},
			MIN_MATCH_LENGTH => {
				MANDATORY=>0, DEFAULT=>30, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Minimum match length for a read to be used"
				},
			MAX_MISMATCH => {
				MANDATORY=>0, DEFAULT=>2, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Maximum mismatch in the alignment for a read to be used"
				},
			MIN_DIST_BETWEEN_BP => {
				MANDATORY=>0, DEFAULT=>10, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Minimum distance between two breakpoints"
				},
			MIN_DIST_FROM_EXTREMITIES => {
				MANDATORY=>0, DEFAULT=>500, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Minimum distance from seq id ends"
				},
			MIN_BREAK_COVERAGE => {
				MANDATORY=>0, DEFAULT=>10, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Minimum number of partially aligned reads to define a breakpoint"
				},
			MAX_COVERAGE => {
				MANDATORY=>0, DEFAULT=>'1000', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Maximum coverage to consider a breakpoint"
				},
		);


my %FILE_EXTENSIONS = ( 
			gff3 => {DESC=> 'GFF3 file', EXT => 'gff3' },		
			 );
			 
=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::MIRAA->new({....})
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
     # Check if parties mapping exist
     my $bam = $self->{OUT_DIR}."/Map/".basename($self->{OUT_DIR}).".".basename($self->{GENOME}).".BOWTIE.sorted.bam";
     $self->{BAM} = $bam if(!$self->{BAM} and -e $bam);
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
   my $sequence = $seq->seq;
   my $seq_length = $seq->length;
   $self->stderr("Calculation $seq_id\n");
   
   
   my $sam = new Bio::DB::Sam( -bam => $self->{BAM}, -fasta => $self->{GENOME} );
   
   my $gff_file = $self->{PATH}."/tmp/$seq_id.gff3";
   open(GFF,">$gff_file") or die "Can not open $gff_file";
   print GFF "# ",join(" ",$self->_command_line),"\n";
   
   my $log_file = $self->{PATH}."/tmp/$seq_id.log";
   open(LOG,">$log_file") or die "Can not open $log_file";
   
   my $err_file = $self->{PATH}."/tmp/$seq_id.err";
   open(ERR,">$err_file") or die "Can not open $err_file";
   
  

   my %breaks;
   foreach my $aln ($sam->get_features_by_location( -seq_id => $seq_id)) {
      next if ($seq_id ne $aln->mate_seq_id);
      next if ( abs($aln->start - $aln->mate_start) > $self->{MAX_INSERT_SIZE} );
      my $mismatch = $aln->get_tag_values('NM');
      
      my $length = $aln->length;
      next if ($length < $self->{MIN_MATCH_LENGTH});
      # end some filtering
      my $aln_start = $aln->start;
      my $aln_end = $aln->end;
      my $cigar_str = $aln->cigar_str;
      #print  $aln_start," ",$aln_end," ",$mismatch," ",$cigar_str,"\n";
      if($cigar_str=~/^(\d+)M(\d+)I\d+M$/) {
         my ($match_part,$ilength) = ($1,$2);
	 next if($ilength <= 25); # min IES length
         $mismatch -= $ilength;
	 next if($mismatch > $self->{MAX_MISMATCH});
	 my $pos = $aln_start+$match_part;
	 $breaks{$pos}->{CIGARS}->{MIM}++;	 
	 $breaks{$pos}->{COUNTS}+=2;
      } elsif($cigar_str=~/^\d+M\d+[SH]$/) {
         my $pos = $aln_end;
	 $breaks{$pos}->{CIGARS}->{MS}++; 
	 $breaks{$pos}->{COUNTS}++;
      
      } elsif($cigar_str=~/^\d+[SH]\d+M$/) {
         my $pos = $aln_start;
	 $breaks{$pos}->{CIGARS}->{SM}++; 
	 $breaks{$pos}->{COUNTS}++;      
      
      }	
   }
   
   # merge break points
   my %merged_breaks;
   my @last_positions;
   foreach my $pos (sort {$a<=>$b} keys %breaks) {
       next if($pos <= $self->{MIN_DIST_FROM_EXTREMITIES} or $pos >= ($seq_length-$self->{MIN_DIST_FROM_EXTREMITIES}));
       
       if(@last_positions and abs($last_positions[$#last_positions]-$pos) >= $self->{MIN_DIST_BETWEEN_BP}) {
	  my ($best_pos,$break_info) = _create_break(\%breaks,@last_positions);	 
	  #die %$break_info;
	  $merged_breaks{$best_pos}=$break_info if($best_pos); 
          @last_positions =($pos);
       } else {
          push @last_positions,$pos;
          
       }
   }
   my ($best_pos,$break_info) = _create_break(\%breaks,@last_positions);	  
   $merged_breaks{$best_pos}=$break_info if($best_pos); 
   
   my $nb_break_points=0;
   
   # THE HEADER MUST START BY #
   foreach my $pos (sort {$a<=>$b} keys %merged_breaks) {	 
      my ($sam_coverage) = $sam->features(-type=> 'coverage:1',-seq_id => $seq_id,-start  => $pos,-end => $pos); 
      my ($average_coverage) = $sam_coverage->coverage;
					      
      if($merged_breaks{$pos}->{COUNTS} < $self->{MIN_BREAK_COVERAGE} or ($self->{MAX_COVERAGE} and $average_coverage > $self->{MAX_COVERAGE})) {
      
        print ERR "BREAK POINTS not properly identified : $seq_id $pos counts=".$merged_breaks{$pos}->{COUNTS}." average_coverage=$average_coverage\n";
      
      } else {
      
         my @attributes = ("ID=BREAK_POINTS_$seq_id\_$pos");
         foreach my $cigar (keys %{$merged_breaks{$pos}->{CIGARS}}) {
            push @attributes,"cigar=$cigar ".$merged_breaks{$pos}->{CIGARS}->{$cigar};
         }
	 push @attributes,"average_coverage=$average_coverage";
         print GFF join("\t",($seq_id,'MIRAA','segment',$pos,$pos,$merged_breaks{$pos}->{COUNTS},'.','.',join(";",@attributes))),"\n";
         $nb_break_points++;
      }
   }
   print LOG "Number of BREAK POINTS for $seq_id : $nb_break_points\n";
      
   
   close GFF;
   close LOG;
   close ERR;
   
   $self->stderr("End calculation $seq_id\n" );

   
}


sub _create_break {
   my ($breaks,@last_positions) = @_;
   
   my ($best_score,$best_pos);
   my %break_info;
   foreach my $p (@last_positions) {
   
      if(!$best_score or $best_score < $breaks->{$p}->{COUNTS}) {
   	 ($best_score,$best_pos) = ($breaks->{$p}->{COUNTS},$p);
      }
      $break_info{COUNTS}+=$breaks->{$p}->{COUNTS};
      foreach my $cs (keys %{$breaks->{$p}->{CIGARS}}) {
   	 $break_info{CIGARS}->{$cs}+=$breaks->{$p}->{CIGARS}->{$cs};
      }
      
   }
   
   return ($best_pos,\%break_info);
}




1;
