package PARTIES::Assembly;
use strict;
use base 'PARTIES::Root';

use PARTIES::Config;
#use PARTIES::Utils;
use File::Basename;

=head1 NAME

 ParTIES Assembly module - Filter reads and assemble them

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			FASTQ1 => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1.1,
				DESCRIPTION=>"Sequencing file containing reads in Forward strand"
				},
			FASTQ2 => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1.2,
				DESCRIPTION=>"Sequencing file containing reads in Reverse strand"
				},				
			BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Mapping on the reference genome"
				},
			MIRAA => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"MIRAA output used to filter reads"
				},
			MIN_MATCH_LENGTH => {
				MANDATORY=>0, DEFAULT=>30, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Minimum match length in the alignment for a read to be used"
				},
			MAX_MISMATCH => {
				MANDATORY=>0, DEFAULT=>2, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Maximummismatch in the alignment for a read to be used"
				},
			NO_ZIP => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Do not compress FASTQ files"
				},	
			KMER => {
				MANDATORY=>1, DEFAULT=>[], TYPE=>'MULTIPLE', RANK => 1,
				DESCRIPTION=>"Velvet : Assembly Kmer parameter (integer odd) ; may be used multiple time"
				},
			INSERT_SIZE => {
				MANDATORY=>0, DEFAULT=>'500', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Velvet : Estimation of the sequencing insert size"
				},
			MIN_CONTIG_LENGTH => {
				MANDATORY=>0, DEFAULT=>'100', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Velvet : Minimum contig length"
				},
			MIN_COVERAGE => {
				MANDATORY=>0, DEFAULT=>'3', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Velvet : Minimum contig coverage"
				},	
			MAX_COVERAGE => {
				MANDATORY=>0, DEFAULT=>'500', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Velvet : Maximum contig coverage"
				},
			SKIP_VELVET_ASSEMBLY => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Skip velvet assembly (just filter reads)"
				},
				
		);


my %FILE_EXTENSIONS = (  );


=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::Assembly->new({....})
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
  
     # default KMERs
     if(-e $self->{BAM} and ($self->{KMER} eq '' or scalar @{$self->{KMER}} == 0 )) {
        my $max_read_length = 0;
	my ($samtools) = (PARTIES::Config->get_program_path('samtools'));
        foreach my $rl (`$samtools view -f 2 $self->{BAM} | head -10 | awk '{ print length(\$10) }'`) {
           chomp $rl;
	   $max_read_length = $rl if($rl > $max_read_length);
        }
     
        # HiSeq
        if($max_read_length < 120) {
           #$self->{KMER} = [41,45,51,55];
           $self->{KMER} = [51];
        # MiSeq
        } elsif($max_read_length > 131) {
           $self->{KMER} = [131];
        } else { die "KMER ? (read_length=$max_read_length)"; }     
     }
  
     my $miraa = $self->{OUT_DIR}."/MIRAA/MIRAA.gff3";
     $self->{MIRAA} = $miraa if(!$self->{MIRAA} and -e $miraa);
  }
  $self->{KMER} = [0] if($self->{SKIP_VELVET_ASSEMBLY} eq 'TRUE');
  
  return 0 if(!$self->SUPER::_check_mandatory_parameters);
  
  foreach my $fastq ($self->{FASTQ1},$self->{FASTQ2}) {
     next if($fastq=~/\.fastq$/);
     print STDERR "ERROR : fastq file (-fastq $fastq) should be a FASTQ file or does not exist\n" ;
     return 0;
  }
  return 1;
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

  # load MIRAA file
  if(-e $self->{MIRAA}) {
     $self->stderr("Read ".basename($self->{MIRAA})." ... \n" );
     $self->{LOADED_MIRAA} = PARTIES::Utils->read_gff_file($self->{MIRAA});
     $self->stderr("Done\n" );  
  }
  
}

=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
   my ($self,$results) = @_;  
   
   
   $self->stderr("Create filtered FASTQ files\n");
   my $no_mac_junc_suffix = 'no_mac_junc';
   my $at_least_one_no_match_suffix = 'at_least_one_no_match';
   
#   my $telomere_pattern = $self->{TELOMERE_PATTERN};

   my @at_least_one_no_match_fastq;
   my @no_mac_junc_fastq;
   
   # create filtered fastq files
   foreach my $fastq ($self->{FASTQ1},$self->{FASTQ2}) {
      
      # first method : a fastq file where at least one read does not match on the reference
      my $at_least_one_no_match_out = $self->{PATH}."/".basename($fastq);
      $at_least_one_no_match_out=~s/\.fastq$/\.$at_least_one_no_match_suffix\.fastq/;
      open(ATLEAST,">$at_least_one_no_match_out") or die $at_least_one_no_match_out;     
      $self->stderr("Write $at_least_one_no_match_out ... \n" );
      push @at_least_one_no_match_fastq,$at_least_one_no_match_out;   
      
      # exclude reads on the mac juntcions
      if($self->{MIRAA}) {
         my $no_mac_junc_out = $self->{PATH}."/".basename($fastq);
         $no_mac_junc_out=~s/\.fastq$/\.$no_mac_junc_suffix\.fastq/;
         open(NOMAC,">$no_mac_junc_out") or die $no_mac_junc_out;   
         $self->stderr("Write $no_mac_junc_out ... \n" );             
         push @no_mac_junc_fastq, $no_mac_junc_out;            
      }
      
      $self->stderr("Read ".basename($fastq)." ... \n" );    
      
      
      # read fastq file and write the filtered one
      open(FASTQ,$fastq) or die "Can not open $fastq";      
      while(<FASTQ>) {
         chomp;
         if($_=~/^\@/) {
	 
	    my $name_line = $_;
	    my ($id,$atts) = split / /,$name_line;
	    $id=~s/^\@//;
	    $id=~s/\/\d$//;  
     	    my $seq = <FASTQ>;
	    chomp $seq;
     	    my $n = <FASTQ>;
     	    my $qual = <FASTQ>;

	    if(!$results->{MAPPED_READS}->{$id} or $results->{MAPPED_READS}->{$id} != 2) {
	          print ATLEAST "$name_line\n$seq\n$n$qual";	    
	    }
	    
	    next if(!$self->{MIRAA});
	    if(!$results->{NO_MAC_JUNCTION}->{$id}) {
	          print NOMAC "$name_line\n$seq\n$n$qual";
	    } 

         }
      }   
      close FASTQ;
      close NOMAC if($self->{MIRAA});     
      close ATLEAST;    
      $self->stderr("Done\n" );
   }
   $self->stderr("Done\n" );
   
   # Assemble the fastq files with velvet
   if($self->{SKIP_VELVET_ASSEMBLY} ne 'TRUE') {
      $self->stderr("Assemble FASTQ files\n");
      my $insert_size =  $self->{INSERT_SIZE};
      my $max_coverage = $self->{MAX_COVERAGE};
      my $min_coverage = $self->{MIN_COVERAGE};
      my $min_contig_length = $self->{MIN_CONTIG_LENGTH};
      
      my $velveth = PARTIES::Config->get_program_path('velveth');
      my $velvetg = PARTIES::Config->get_program_path('velvetg');
      
      # you can use several kmers to assemble
      foreach my $kmer (@{$self->{KMER}}) {
         $self->stdlog("\nKMER = $kmer\n");
         # method : no filters
         my $method = 'no_filter';
         my $velvet_out_dir = $self->{PATH}."/VELVET_$kmer\_$method";
         $self->stderr("velvet $velvet_out_dir kmer=$kmer ".$self->{FASTQ1}." ".$self->{FASTQ2}."\n");
         $self->stdlog("velveth $velvet_out_dir $kmer -fastq -shortPaired -separate ".$self->{FASTQ1}." ".$self->{FASTQ2}."\n");
         system("$velveth $velvet_out_dir $kmer -fastq -shortPaired -separate ".$self->{FASTQ1}." ".$self->{FASTQ2});
         $self->stdlog("velvetg $velvet_out_dir -very_clean yes -scaffolding no -max_coverage $max_coverage -exp_cov auto -ins_length $insert_size -min_contig_lgth $min_contig_length -cov_cutoff $min_coverage\n\n");
         system("$velvetg $velvet_out_dir -very_clean yes -scaffolding no -max_coverage $max_coverage -exp_cov auto -ins_length $insert_size -min_contig_lgth $min_contig_length -cov_cutoff $min_coverage");
         my $contig_file = "$velvet_out_dir/VELVET_$kmer\_$method\_contigs.fa";
         system("mv $velvet_out_dir/contigs.fa $contig_file");      

         # method : at least one no match
         $method = 'at_least_one_no_match';
         $velvet_out_dir = $self->{PATH}."/VELVET_$kmer\_$method";
         $self->stderr("velvet $velvet_out_dir kmer=$kmer ".$at_least_one_no_match_fastq[0]." ".$at_least_one_no_match_fastq[1]."\n");
         $self->stdlog("velveth $velvet_out_dir $kmer -fastq -shortPaired -separate ".$at_least_one_no_match_fastq[0]." ".$at_least_one_no_match_fastq[1]."\n");
         system("$velveth $velvet_out_dir $kmer -fastq -shortPaired -separate ".$at_least_one_no_match_fastq[0]." ".$at_least_one_no_match_fastq[1]);
         $self->stdlog("velvetg $velvet_out_dir -very_clean yes -scaffolding no -max_coverage $max_coverage -exp_cov auto -ins_length $insert_size -min_contig_lgth $min_contig_length -cov_cutoff $min_coverage\n\n");
         system("$velvetg $velvet_out_dir -very_clean yes -scaffolding no -max_coverage $max_coverage -exp_cov auto -ins_length $insert_size -min_contig_lgth $min_contig_length -cov_cutoff $min_coverage");
         $contig_file = "$velvet_out_dir/VELVET_$kmer\_$method\_contigs.fa";
         system("mv $velvet_out_dir/contigs.fa $contig_file");
       
       
         next if(!$self->{MIRAA});              
         # method : no MAC junctions
         $method = 'no_mac_junctions';
         $velvet_out_dir = $self->{PATH}."/VELVET_$kmer\_$method";
         $self->stderr("velvet $velvet_out_dir kmer=$kmer ".$no_mac_junc_fastq[0]." ".$no_mac_junc_fastq[1]."\n");
         $self->stdlog("velveth $velvet_out_dir $kmer -fastq -shortPaired -separate ".$no_mac_junc_fastq[0]." ".$no_mac_junc_fastq[1]."\n");
         system("$velveth $velvet_out_dir $kmer -fastq -shortPaired -separate ".$no_mac_junc_fastq[0]." ".$no_mac_junc_fastq[1]);
         $self->stdlog("velvetg $velvet_out_dir -very_clean yes -scaffolding no -max_coverage $max_coverage -exp_cov auto -ins_length $insert_size -min_contig_lgth $min_contig_length -cov_cutoff $min_coverage\n\n");
         system("$velvetg $velvet_out_dir -very_clean yes -scaffolding no -max_coverage $max_coverage -exp_cov auto -ins_length $insert_size -min_contig_lgth $min_contig_length -cov_cutoff $min_coverage");
         $contig_file = "$velvet_out_dir/VELVET_$kmer\_$method\_contigs.fa";
         system("mv $velvet_out_dir/contigs.fa $contig_file");
      
      
      }
      $self->stderr("Done\n" );
   }
   
   # gzip the fastq files
   if($self->{NO_ZIP} ne 'TRUE') {
      $self->stderr("Zip files ... ");
      system(PARTIES::Config->get_program_path('gzip')." $self->{PATH}/*.fastq &");
   }
   
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
   
   my $sam = new Bio::DB::Sam( -bam => $self->{BAM}, -fasta => $self->{GENOME} );
   
   
   my $log_file = $self->{PATH}."/tmp/$seq_id.log";
   open(LOG,">$log_file") or die "Can not open $log_file";
   
   my $err_file = $self->{PATH}."/tmp/$seq_id.err";
   open(ERR,">$err_file") or die "Can not open $err_file";

#   my $tab_file = $self->{PATH}."/tmp/$seq_id.tab";
#   open(TAB,">$tab_file") or die "Can not open $tab_file";
   
   my %results; 
   
   return %results if(!$self->{LOADED_MIRAA} or !$self->{LOADED_MIRAA}->{$seq_id});
   
   # No MAC junctions around break points 
   foreach my $bp (@{$self->{LOADED_MIRAA}->{$seq_id}}) {
   
     my $pos = $bp->{start};
     foreach my $a ($sam->get_features_by_location(-seq_id => $seq_id, -start => $pos, -end => $pos+1)) {
   	 my $cigar_str = $a->cigar_str;
   	 my $mismatch = $a->get_tag_values('NM');
   	 my $length = $a->length;
   	 next if ($cigar_str !~ /^\d+M$/ ); # FIXME what about indels of 1 ?
   	 next if ($length < $self->{MIN_MATCH_LENGTH});
   	 next if ($mismatch/$length > $self->{MAX_MISMATCH} /100); # FIXME what about indels of 1 ?
   	 next if (abs($a->start - $pos) < 5 || abs($a->end - $pos) < 5);
   	 my $id = $a->query->name;
	 #print TAB $id,"\n";
	 $results{NO_MAC_JUNCTION}->{$id}++;
      }
   }
   
   # At least one no match in reads
   my %mapped_reads;
   foreach my $a ($sam->get_features_by_location(-seq_id => $seq_id )) {
      my $cigar_str = $a->cigar_str;
      next if ($cigar_str !~ /^\d+M$/ ); 
      my $id = $a->query->name;
      $results{MAPPED_READS}->{$id}++;
   }
   
   
   
#   close TAB;
   close LOG;
   close ERR;
   
   $self->stderr("End calculation $seq_id\n" );
   return %results;
}


1;
