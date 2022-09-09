package PARTIES::Map;
use strict;
use base 'PARTIES::Root';
use PARTIES::Config;
use File::Basename;

=head1 NAME

 ParTIES Map module - Map reads on a reference using bowtie2

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
			MAX_INSERT_SIZE => {
				MANDATORY=>0, DEFAULT=>'1000', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Max sequencing insert size"
				},
			INDEX_GENOME => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Force to index genome with bowtie2"
				},	
			USE_INSERT_REFERENCE => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Try to use the result of the Insert module as the reference genome"
				},	
		);


my %FILE_EXTENSIONS = ( );


=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::Map->new({....})
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
  if($self->{USE_INSERT_REFERENCE}) {
     $self->{GENOME} = $self->{OUT_DIR}."/Insert/Insert.fa" if(-e $self->{OUT_DIR}."/Insert/Insert.fa");
  }
  
  return 0 if(!$self->SUPER::_check_mandatory_parameters);
  
  foreach my $fastq ($self->{FASTQ1},$self->{FASTQ2}) {
     next if($fastq=~/\.fastq$/ and -e $fastq);
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

  
  if($self->{INDEX_GENOME} eq 'TRUE') {
     $self->stderr("Index genome ...\n" );
     $self->stderr("bowtie2-build $self->{GENOME}\n");
     $self->{BT2_INDEX} = $self->{PATH}."/".basename($self->{GENOME});
     my ($bowtie2_build,$samtools) = (PARTIES::Config->get_program_path('bowtie2-build'),PARTIES::Config->get_program_path('samtools'));
     system("$bowtie2_build $self->{GENOME} $self->{BT2_INDEX} > /dev/null 2>&1");   
     system("$samtools faidx $self->{GENOME}");    
     $self->stderr("Done\n" );
  } else {
    $self->{BT2_INDEX} = $self->{GENOME};
  }
  die "ERROR bowtie2 index genome does not exist for ",$self->{GENOME}," (use -index_genome)\n" if(!-e $self->{BT2_INDEX}.".1.bt2") ;
}

=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
  my ($self,$results) = @_;  
  
  my ($bowtie2,$samtools) = (PARTIES::Config->get_program_path('bowtie2'),PARTIES::Config->get_program_path('samtools'));
      
  $self->stderr("Mapping on ".basename($self->{GENOME})." ... \n" );
  my $out_bam = $self->{PATH}."/".basename($self->{OUT_DIR}).".".basename($self->{GENOME}).".BOWTIE.sorted.bam";
  $self->stdlog("$bowtie2 --threads $self->{THREADS} --local -x $self->{BT2_INDEX} -1 $self->{FASTQ1} -2 $self->{FASTQ2} -X $self->{MAX_INSERT_SIZE}\n");
  system("$bowtie2 --threads $self->{THREADS}  --quiet --local -x $self->{BT2_INDEX} -1 $self->{FASTQ1} -2 $self->{FASTQ2} -X $self->{MAX_INSERT_SIZE} | $samtools view -uS - | $samtools sort -o $out_bam - > /dev/null 2>&1");
  system("$samtools index $out_bam");
  $self->stderr("Done\n" );
  
   
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
  return;
}


1;
