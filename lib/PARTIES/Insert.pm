package PARTIES::Insert;
use strict;
use base 'PARTIES::Root';
use PARTIES::Config;
use File::Basename;

=head1 NAME

 ParTIES Insert module - Insert IES within a genome to create an IES containing reference

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			IES => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"IES file, usually MICA output"
				},
			PREFIX => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Prefix for new seq id names"
				},
			SUFFIX => {
				MANDATORY=>0, DEFAULT=>'_with_IES', TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Suffix for new seq id names"
				},
			LOW_CASE => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Should IES sequences be written in lowercase"
				},
		);



my %FILE_EXTENSIONS = ( 
			gff3 => {DESC=> 'GFF3 file', EXT => 'gff3' },
			fa => {DESC=> 'Fasta file', EXT => 'fa' },
			chr_cor => {DESC=> 'Chromosome correspondancies', EXT => 'chr_cor' },					
			 );
 

=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::Insert->new({....})
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
     my $mica = $self->{OUT_DIR}."/MICA/MICA.gff3";
     $self->{IES} = $mica if(!$self->{IES} and -e $mica);
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
  
  # READ FILES
    
  # load IES file
  $self->stderr("Read ".basename($self->{IES})." ...\n" );
  $self->{LOADED_IES} = PARTIES::Utils->read_gff_file($self->{IES});
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
  
  
   my $seq_id = $seq->id;  
   $self->stderr("Calculation $seq_id\n");
   
   
   my $log_file = $self->{PATH}."/tmp/$seq_id.log";
   open(LOG,">$log_file") or die "Can not open $log_file";
   
   my $err_file = $self->{PATH}."/tmp/$seq_id.err";
   open(ERR,">$err_file") or die "Can not open $err_file";
   
   my $gff_file = $self->{PATH}."/tmp/$seq_id.gff3";
   open(GFF,">$gff_file") or die "Can not open $gff_file";
   print GFF "# ",join(" ",$self->_command_line),"\n";

   my $fasta_file = $self->{PATH}."/tmp/$seq_id.fa";
   open(FA,">$fasta_file") or die "Can not open $fasta_file";
   
   my $new_seq_id = join("",$self->{PREFIX},$seq_id,$self->{SUFFIX});
   
   my $chr_corr_file = $self->{PATH}."/tmp/$seq_id.chr_cor";
   open(CORR,">$chr_corr_file") or die "Can not open $chr_corr_file";   
   print CORR join("\t",($seq_id,$new_seq_id)),"\n";
   close CORR;
   
   my $sequence = $seq->seq;
   my $length = $seq->length;

   my $new_sequence ;
   my $last_indel_position=0;
   my $ies_length=0;
   if($self->{LOADED_IES}->{$seq_id}) {
      foreach my $ies (sort {$a->{start}<=>$b->{start}} @{$self->{LOADED_IES}->{$seq_id}}) {
         my $indel_position =  $ies->{start};	  
         
	 my ($indel_sequence) = @{$ies->{attributes}->{sequence}};	 
	 $indel_sequence = lc($indel_sequence) if($self->{LOW_CASE} eq 'TRUE');
	 
	 my ($id) = @{$ies->{attributes}->{ID}};
	 print ERR "ERROR : $seq_id => No consensus seq for $id\n" if(!$indel_sequence);
	 
	 my $subseq = substr($sequence,$last_indel_position,($indel_position-$last_indel_position)-1);
	 $subseq.=$indel_sequence;
	 #die $new_sequence," ",$last_indel_position if($last_indel_position>2000);
	  
      
         my @attributes;
         foreach my $key (sort keys %{$ies->{attributes}}) {
            foreach my $value (@{$ies->{attributes}->{$key}}) {
               push @attributes,"$key=$value";
            }
         }

         print GFF join("\t",($new_seq_id,$ies->{source},$ies->{type},($indel_position+$ies_length),($indel_position+$ies_length+length($indel_sequence)-1),'.','.',".",join(";",@attributes))),"\n";
      
         $ies_length += length($indel_sequence);
	 
         $new_sequence.=$subseq;
      
      
      
      
         $last_indel_position=$indel_position-1;
      }
   }   
   $new_sequence.=substr($sequence,$last_indel_position);
   my $new_length = length($new_sequence);
   $new_sequence=~s/(.{60})/$1\n/g;
   print FA ">$new_seq_id mac_length=$length mic_length=$new_length ies_length=$ies_length\n";
   print FA $new_sequence,"\n";
   
   # check if the length is consistent
   print ERR "ERROR : $seq_id => diff length $seq_id $length != ".($new_length-$ies_length)."\n" if($length != ($new_length-$ies_length));

   
   close GFF;
   close FA;
   close LOG;
   close ERR;
   
   $self->stderr("End calculation $seq_id\n" );

   
}


1;
