package PARTIES::Run;
use strict;
use base 'PARTIES::Root';

use PARTIES::Config;

=head1 NAME

 ParTIES Run module - Run ParTIES using the configuration file

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			
			CONFIG => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Configuration file"
				},					

		);

my %FILE_EXTENSIONS = ( );


=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::Run->new({....})
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
  
 # load all parameters
  
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
  
  
  return 0 if(!$self->SUPER::_check_mandatory_parameters);
  
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
  #$self->SUPER::_init;
  
  
  my @pipeline = PARTIES::Utils->read_config_file($self->{CONFIG});
  
  my ($out_dir,$genome) = ($self->{OUT_DIR},$self->{GENOME});
  foreach my $process (@pipeline) {
     my $mode = $process->{MODE};
     my @params = ();
     
     
     foreach my $t (qw(GENERAL_PARAMS PARAMS)) {
        next if(!$process->{$t});
        foreach my $key (keys %{$process->{$t}}) {
	   foreach my $value (@{$process->{$t}->{$key}}) {
	      next if($value eq 'FALSE');
	      push @params, ($value eq 'TRUE') ? $key : "$key $value";
	   }
	}
     }
     $self->stderr("$0 $mode -out_dir $out_dir -genome $genome ".join(" ",@params)."\n");
     my $r = system("$0 $mode -out_dir $out_dir -genome $genome ".join(" ",@params));  
     exit 1 if($r == 0);
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

}


1;
