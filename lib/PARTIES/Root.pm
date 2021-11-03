package PARTIES::Root;
use strict;
use POSIX qw(strftime);
use Time::HiRes qw(usleep gettimeofday tv_interval);
use FileHandle;
use File::Basename;

=head1 NAME

 ParTIES Root module

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# COMMON PARAMETERS
my %ROOT_PARAMETERS = (
			VERSION => {
				MANDATORY=>1, DEFAULT=>'?', TYPE=>'VALUE', RANK => 5,
				DESCRIPTION=>"Version"
				},
			GENOME => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 0,
				DESCRIPTION=>"Reference genome"
				},
			OUT_DIR => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 0,
				DESCRIPTION=>"Out directory"
				},
			SEQ_ID => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 4,
				DESCRIPTION=>"Calculate only on this seq_id"
				},			
            LIST_OF_SEQ => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 4,
				DESCRIPTION=>"File with a list of sequences to analyze"
				},		
			MIN_SEQ_LENGTH => {
				MANDATORY=>0, DEFAULT=>0, TYPE=>'VALUE', RANK => 4,
				DESCRIPTION=>"Calculate only on sequences longest than this parameter"
				},		
			THREADS => {
				MANDATORY=>0, DEFAULT=>'6', TYPE=>'VALUE', RANK => 4,
				DESCRIPTION=>"Number of threads"
				},
			VERBOSE => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 4,
				DESCRIPTION=>"Verbose"
				},			
			QUIET => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 4,
				DESCRIPTION=>"Quiet mode"
				},
			KEEP_TMP => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 4,
				DESCRIPTION=>"Keep tmp files"
				},			
			FORCE => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 4,
				DESCRIPTION=>"Overrides the results when the command has already been executed"
				},			
			AUTO => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 4,
				DESCRIPTION=>"Enable auto-detection of parameters (usually files generated from other modules)"
				},			
			CHECK_CONFIG => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 4,
				DESCRIPTION=>"Check mandatory parameters (do not take into account auto-detection)"
				},
		);
		

my $STDLOG_FILE = "stdout.log";
my $STDERR_FILE = "stderr.log";

=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::Root->new({....})
 Function: Create a factory object
 Returns : the factory object

=cut

sub new {
  my ($class,$ref) = @_;
  
  # Check the class
  $class = ref($class) || $class;
  
  # Link between the object and the class
  my $self = bless {},$class;
  
  my %root_parameters = %{$self->get_root_parameters};
  foreach my $param (keys %root_parameters) {
     $self->{$param} = ((defined $ref->{"-".lc($param)} and $ref->{"-".lc($param)} ne '') and (ref($ref->{"-".lc($param)}) ne 'ARRAY' or scalar @{$ref->{"-".lc($param)}} !=0)) ? $ref->{"-".lc($param)} : $root_parameters{$param}->{DEFAULT}; 
	die "$param should be TRUE or FALSE not : $self->{$param}\n" if($root_parameters{$param}->{TYPE} eq 'BOOLEAN' and $self->{$param} ne 'TRUE' and $self->{$param} ne 'FALSE');
     
  }
  
  
  return $self;
}




=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
  my ($self) = @_;  
  
  $self->{LOG}->close;
  my $base = $self->get_mode;
  
  foreach my $ext (qw(log err)) {
     next if(!$ext);
     #chomp $ext;
     if($ext eq 'log') {
        $self->stderr("Write file : $self->{PATH}/$STDLOG_FILE\n");
        system("cat $self->{PATH}/tmp/*.$ext 2> /dev/null >> $self->{PATH}/$STDLOG_FILE");
     } elsif($ext eq 'err') {
        $self->stderr("Write file : $self->{PATH}/$STDERR_FILE\n");
        system("cat $self->{PATH}/tmp/*.$ext 2> /dev/null > $self->{PATH}/$STDERR_FILE");     
     }
  } 
  
  my %file_ext = %{$self->get_file_extensions};
  foreach my $key (keys %file_ext) {
     my $ext = ($file_ext{$key}->{EXT}) ? $file_ext{$key}->{EXT} : $key;
     next if(!$ext);
     $self->stderr("Write file : $self->{PATH}/$base.$ext\n");
     system("cat $self->{PATH}/tmp/*.$key 2> /dev/null | grep '^#' | head -1  > $self->{PATH}/$base.$ext");
     system("echo \"\$(echo '##gff-version 3' | cat - $self->{PATH}/$base.$ext)\" > $self->{PATH}/$base.$ext") if($ext=~/gff3$/);
     system("cat $self->{PATH}/tmp/*.$key 2> /dev/null | grep -v '^#' >> $self->{PATH}/$base.$ext");     
     
  }
  $self->remove_temporary_files;   
}


sub remove_temporary_files {
  my ($self) = @_;  

  if($self->{KEEP_TMP} eq 'FALSE') {
      $self->stderr("Remove temporary files : $self->{PATH}/tmp\n");
      system("rm -rf $self->{PATH}/tmp");  
  }    

}

###################################
# PRIVATE FUNCTIONS
###################################

sub _init {
  my ($self,$version) = @_;  
  
  $self->{TIME_POINT_0} = [gettimeofday];
  $self->{LAST_TIME_POINT} = $self->{TIME_POINT_0};
  
  $self->stderr("Initialization\n");
  # create out directory
  my $path = $self->{OUT_DIR}."/".$self->get_mode;
  
  die "ERROR directory already exists $path ; use -force option to skip\n" if($self->{FORCE} eq 'FALSE' and -e $path);
  system("rm -rf $path/tmp") if($self->{FORCE} eq 'TRUE' and -e $path);
  $self->stderr("Create directory $path\n");
  system("mkdir -p $path/tmp");
  $self->{PATH}=$path;

  $self->{LOG} = FileHandle->new(">$self->{PATH}/$STDLOG_FILE");
  $self->stdlog("# VERSION : ".$self->{VERSION}."\n"); 
  $self->stdlog("# COMMAND LINE :\n".join(" \\\n",$self->_command_line)."\n\n");
  $self->stdlog("# PATH OUT DIR : ".$self->{PATH}."\n\n"); 
  
  

  
}



# Display command line
sub _command_line {
   my ($self) = @_;

   my @command = "$0 ".$self->get_mode;

   my %all_parameters=(%{$self->get_root_parameters}, %{$self->get_parameters});
   
   foreach my $param (sort keys %all_parameters) {
      next if($param eq 'VERSION');
      my $type = $all_parameters{$param}->{TYPE} ;   
      if($type eq 'VALUE') {
        my $param_value = $self->{$param};
	$param_value = "'$param_value'" if($param_value=~/ /);
	next if($param_value eq '');
        push @command," -".lc($param)." ".$param_value; 
      } elsif($type eq 'MULTIPLE') {
         foreach my $param_value (@{$self->{$param}}) {
	    $param_value = "'$param_value'" if($param_value=~/ /);
            push @command, " -".lc($param)." ".$param_value;	
	 }
     } elsif($type eq 'BOOLEAN') {
        push @command," -".lc($param) if($self->{$param} eq 'TRUE');    
     }
   }   
   return @command;      
}



=head2 _check_mandatory_parameters

 Title   : "Private" function _check_mandatory_parameters
 Usage   : $factory->_check_mandatory_parameters
 Function: Check the mandatory parameters or deduce them with -auto
 Returns : 1 (success) or 0 (error)

=cut

sub _check_mandatory_parameters {
  my ($self) = @_;  
  my %all_parameters=(%{$self->get_root_parameters}, %{$self->get_parameters});
  foreach my $param (keys %all_parameters) {
     if($all_parameters{$param}->{MANDATORY} and ($self->{$param} eq '' or (ref($self->{$param}) eq 'ARRAY' and scalar @{$self->{$param}} ==0 ))) {
        print STDERR "\n# ERROR -",lc($param)," IS MANDATORY\n\n";
        return 0;
     }
  }
  
  if($self->{GENOME} !~/\.fa$/i and $self->{GENOME} !~/\.fasta$/i)  {
     print STDERR "Genome file (-genome) should be a FASTA file";
     return 0;
  }
    
  return 1;
}


# Print usage
sub _usage {
  my ($self) = @_;
  
  
  my %all_parameters=(%{$self->get_root_parameters}, %{$self->get_parameters});
  print STDERR "USAGE : \n$0 ".$self->get_mode."\n";
  my $last_rank;
  foreach my $param (sort {$all_parameters{$a}->{RANK} <=> $all_parameters{$b}->{RANK} } keys %all_parameters) {
      next if($param eq 'VERSION');
     my $mandatory = '(MANDATORY)' if($all_parameters{$param}->{MANDATORY});
     my $type = $all_parameters{$param}->{TYPE};
     my $desc = $all_parameters{$param}->{DESCRIPTION};
     
     print STDERR "\n" if($last_rank ne $all_parameters{$param}->{RANK});
     if($type eq 'VALUE') {
        my $value = ($self->{$param} ne '') ? $self->{$param} : $all_parameters{$param}->{DEFAULT};
        print STDERR "\t-".lc($param)." [".$value."]\t".$desc." $mandatory\n";
     } elsif($type eq 'MULTIPLE') {
        print STDERR "\t{\n";
	my @values = ($self->{$param}) ? @{$self->{$param}} : @{$all_parameters{$param}->{DEFAULT}};
	print STDERR "\t\t-".lc($param)." []\n" if(scalar @values == 0);
        foreach my $value (@values) {
	   print STDERR "\t\t-".lc($param)." [".$value."]\n";
	}
	print STDERR "\t} (MULTIPLE VALUES) $mandatory\t".$desc."\n";	
     } elsif($type eq 'BOOLEAN') {
        my $value = ($self->{$param}) ? $self->{$param} : $all_parameters{$param}->{DEFAULT};
	print STDERR "\t(-".lc($param)." $mandatory (default $all_parameters{$param}->{DEFAULT}))\t$desc\n";   
     }
     else { die $param,' ',$type; }
     $last_rank = $all_parameters{$param}->{RANK}; 
  }   
  print STDERR "\n";
 
  exit 0;
}

###################################
# PROTOTYPES
###################################


=head2 calculate

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut


sub calculate {
  my ($self,$seq) = @_;  
  die "Not managed for $self";
}



###################################
# GENERIC/UTILS FUNCTIONS
###################################

sub get_root_parameters {
  my ($self) = @_;  
  return \%ROOT_PARAMETERS;
}


sub get_mode {
  my ($self) = @_;  
  my $module_name = ref($self);
  if( $module_name =~/::(\S+)$/) {
     return $1;
  }
  die $module_name;
}

sub verbose {
   my ($self) = @_;
   return ($self->{VERBOSE} eq 'TRUE') ? 1 : 0;
}
sub auto {
   my ($self) = @_;
   return ($self->{AUTO} eq 'TRUE') ? 1 : 0;
}
sub stdlog {
   my ($self,$message) = @_;
   my $fh = $self->{LOG};
   print $fh "$message";
}

sub stderr {
   my ($self,$message) = @_;
   return if($self->{QUIET} eq 'TRUE');
   my $now = [gettimeofday];
   my $duration = sec2human(tv_interval ( $self->{TIME_POINT_0}, $now));
   my $elapsed = sec2human(tv_interval ( $self->{LAST_TIME_POINT}, $now)) ;
   print STDERR "# [".$self->get_mode."] ".strftime("%H:%M:%S",localtime(time))," - $duration (elapsed=$elapsed) : $message" if($self->verbose);
   $self->{LAST_TIME_POINT} = $now;



}

sub sec2human {
    my $secs = shift;
    if    ($secs >= 365*24*60*60) { return sprintf '%.1fy', $secs/(365*24*60*60) }
    elsif ($secs >=     24*60*60) { return sprintf '%.1fd', $secs/(    24*60*60) }
    elsif ($secs >=        60*60) { return sprintf '%.1fh', $secs/(       60*60) }
    elsif ($secs >=           60) { return sprintf '%.1fm', $secs/(          60) }
    else                          { return sprintf '%.1fs', $secs                }
}

1;
