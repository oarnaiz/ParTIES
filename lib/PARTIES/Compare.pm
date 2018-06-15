package PARTIES::Compare;
use strict;
use base 'PARTIES::Root';
use PARTIES::Config;
use List::Util qw(max min);
use FindBin qw($Bin);
use File::Basename;

=head1 NAME

 ParTIES Compare module -  Compare IES/InDel datasets

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			REFERENCE_SET => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Reference gff3 file, usually MICA or MILORD result"
				},			
			CURRENT_SET => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"gff3 to compare to the reference set, usually MICA or MILORD result"
				},							
			MAX_DIST => {
				MANDATORY=>1, DEFAULT=>'0', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Maximum distance between two elements to link them together"
				},
			TAB => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 2,
				DESCRIPTION=>"Write results in tabulated format"
				},
		);



my %FILE_EXTENSIONS = ( 
			IES => {DESC=> 'IES GFF3 file', EXT => 'current.gff3' },
			IES_PLUS => {DESC=> 'IES GFF3 file', EXT => 'reference.gff3' }	
			 );

=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::Compare->new({....})
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
  $self->{BIN}=$Bin;  
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
  my $genome = PARTIES::Utils->read_fasta_file($self->{GENOME});
    
    
    
  # load CURRENT_SET file
  $self->stderr("Read ".basename($self->{CURRENT_SET})." ...\n" );
  $self->{LOADED_CURRENT_SET} = PARTIES::Utils->read_gff_file($self->{CURRENT_SET});
  foreach my $seq_id (keys %{$self->{LOADED_CURRENT_SET}}) {
     next if($genome->{$seq_id});
     die "No seq_id $seq_id in genome=".basename($self->{GENOME})." (file=".basename($self->{CURRENT_SET}).") Wrong genome ?\n";
  }
  $self->stderr("Done\n" );
  
  $self->stderr("Read ".basename($self->{REFERENCE_SET})." ...\n" );
  $self->{LOADED_REFERENCE_SET} = PARTIES::Utils->read_gff_file($self->{REFERENCE_SET});
  foreach my $seq_id (keys %{$self->{LOADED_IES}}) {
     next if($genome->{$seq_id});
     die "No seq_id $seq_id in genome=".basename($self->{GENOME})." (file=".basename($self->{REFERENCE_SET}).") Wrong genome ?\n";
  }
  $self->stderr("Done\n" );  

}

=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
   my ($self,$results) = @_;  
 
   $self->SUPER::finish($results);
   
   if($self->{TAB} eq 'TRUE') {
      $self->stderr("Write tab file ...\n" );
      my $gff2tab=$self->{BIN}.'/utils/gff2tab.pl';
   
      foreach my $key (keys %FILE_EXTENSIONS) {
         my $fext = ($FILE_EXTENSIONS{$key}->{EXT}) ? $FILE_EXTENSIONS{$key}->{EXT} : $key;
         next if($fext!~/\.gff3$/);
         my $gff_file = $self->{PATH}."/".$self->get_mode.'.'.$fext;
         my $tab_file = $gff_file ;
         $tab_file=~s/\.gff3$/\.tab/;
         system("$gff2tab -gff $gff_file > $tab_file");
      }
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
   
   
   my $log_file = $self->{PATH}."/tmp/$seq_id.log";
   open(LOG,">$log_file") or die "Can not open $log_file";
   
   my $err_file = $self->{PATH}."/tmp/$seq_id.err";
   open(ERR,">$err_file") or die "Can not open $err_file";
      
   my $wa_gff_file = $self->{PATH}."/tmp/$seq_id.IES";
   open(WA,">$wa_gff_file") or die "Can not open $wa_gff_file";
   print WA "# ",join(" ",$self->_command_line),"\n";

   my %log; 
   my %links_with_ref_ies;
   if($self->{LOADED_CURRENT_SET}->{$seq_id}) {
      foreach my $ies (sort {$a->{start}<=>$b->{start}} @{$self->{LOADED_CURRENT_SET}->{$seq_id}}) {
          my ($start,$end) = ($ies->{start},$ies->{end});  
	  my ($ies_id) = @{$ies->{attributes}->{ID}};  
	  
	  my @overlaping_features = $self->_find_overlaping_features($ies,$self->{LOADED_REFERENCE_SET}->{$seq_id});
          
	  # get the closest/best ies in reference dataset (by default in in the overlapping ones)
          my ($closest_feature,$boundary_relation,$identical,$min_dist) = $self->_find_closest_feature($ies,$self->{LOADED_REFERENCE_SET}->{$seq_id},@overlaping_features);	  
	  my ($closest_feature_id) = ($closest_feature) ? @{$closest_feature->{attributes}->{ID}} : ('NA');
          
	  # create a data struct with all the information 
	  my %closest_feature_info = (
	    			 found_related_features => (@overlaping_features) ? 'TRUE' : 'FALSE',
	  			 nb_related_features => scalar @overlaping_features,
				 starts_distance => ($closest_feature) ? abs($closest_feature->{start} - $start) : 'NA',
				 ends_distance => ($closest_feature) ? abs($closest_feature->{end} - $end) : 'NA',
				 min_boundaries_distance => ($closest_feature) ? $min_dist : 'NA',
				 boundary_relation => ($closest_feature) ? $boundary_relation : 'NA',
				 closest_feature_id => $closest_feature_id,
				 identical => $identical,
	  			 );	  
	  # stats for report
	  $log{$boundary_relation}++;
	  
	  # add info to gff object
	  foreach my $key (keys %closest_feature_info) { push @{$ies->{attributes}->{$key}},$closest_feature_info{$key}; }
	  
	  print WA PARTIES::Utils->gff_line($ies),"\n";
	  
	  if(@overlaping_features and $closest_feature) {	     
	     $links_with_ref_ies{$closest_feature_id}->{$ies_id} = \%closest_feature_info;	
	     
	  }
      }
      print LOG "For the seq_id $seq_id ";
      foreach my $k (sort {$log{$b}<=>$log{$a}} keys %log) {
         print LOG " $k : $log{$k} ",sprintf('%.2f',($log{$k}*100/scalar @{$self->{LOADED_CURRENT_SET}->{$seq_id}})),"% ";
   
      }
      print LOG "\n";
   } else {
      print ERR "No elements for the seq_id : $seq_id in the file $self->{CURRENT_SET}\n";
   }
   close WA;
   
   

   
   # write ies ref file
   my $wb_gff_file = $self->{PATH}."/tmp/$seq_id.IES_PLUS";
   open(WB,">$wb_gff_file") or die "Can not open $wb_gff_file";
   print WB "# ",join(" ",$self->_command_line),"\n";

   if($self->{LOADED_REFERENCE_SET}->{$seq_id}) {
      foreach my $ies (sort {$a->{start}<=>$b->{start}} @{$self->{LOADED_REFERENCE_SET}->{$seq_id}}) {  
         
	 my ($ies_id) = @{$ies->{attributes}->{ID}};  
	  
	 my %closest_feature_info = (
	 			found_related_features => 'FALSE',
	 			nb_related_features => 0,
	        		min_boundaries_distance => 'NA',
	        		identical => 'FALSE',
	 			);	 
	 if($links_with_ref_ies{$ies_id}) {
	    $closest_feature_info{found_related_features} = 'TRUE';
	    
	    foreach my $i (keys %{$links_with_ref_ies{$ies_id}}) {
	       $closest_feature_info{identical} = 'TRUE' if($links_with_ref_ies{$ies_id}->{$i}->{identical} eq 'TRUE');
	       $closest_feature_info{nb_related_features}++;
	       $closest_feature_info{min_boundaries_distance} = $links_with_ref_ies{$ies_id}->{$i}->{min_boundaries_distance} if($closest_feature_info{min_boundaries_distance} eq 'NA' 
	       																or $links_with_ref_ies{$ies_id}->{$i}->{min_boundaries_distance} < $closest_feature_info{min_boundaries_distance} );	       
	    }
	 
	 }
	 
	 
	 # add info to gff object
	 foreach my $key (keys %closest_feature_info) { push @{$ies->{attributes}->{$key}},$closest_feature_info{$key}; }	
	 
	 print WB PARTIES::Utils->gff_line($ies),"\n";
          
      }
   } else {
      print ERR "No elements for the seq_id : $seq_id in the file $self->{REFERENCE_SET}\n";
   }
   
   close WB;
   close LOG;
   close ERR;
   
   $self->stderr("End calculation $seq_id\n" );
  
}


# find features which overlap an IES
sub _find_overlaping_features {
   my ($self,$ies,$indels) = @_;
   my ($seq_id,$start,$end) = ($ies->{seq_id},$ies->{start},$ies->{end});
   my @overlaping_features;
   return @overlaping_features if(!$indels);
   
   my $min_overlapping = 1;
   foreach my $ref_candidate (@{$indels}) {
      next if($ref_candidate->{end} < $start or $ref_candidate->{start} > $end);
      if(abs($ref_candidate->{end} - $start) >= $min_overlapping or ($ref_candidate->{start} - $end) >= $min_overlapping) {
         push @overlaping_features,$ref_candidate;
      }
   }
   return @overlaping_features;
}

# find the closest element
sub _find_closest_feature {
   my ($self,$ies,$all_features,@overlapping_features) = @_;
   my ($start,$end) = ($ies->{start},$ies->{end});
   my ($ies_id) = @{$ies->{attributes}->{ID}};
   my ($ies_sequence) = @{$ies->{attributes}->{sequence}};
   my ($closest_feature,$min_dist);
   my $boundary_relation = 'no_overlap';
   
   my %boundary_relation_scores = (no_overlap => 0, only_overlap=> 1, internal => 2, external =>2, 
   					right_external => 3 , left_external => 3, right_internal => 3 , left_internal => 3,
					identical => 4);
   
   
   if(@overlapping_features) {
   
      foreach my $feat (@overlapping_features) {
         my ($feat_sequence) = @{$feat->{attributes}->{sequence}};
         my ($feat_id) = @{$feat->{attributes}->{ID}};
         my ($feat_start,$feat_end) = ($feat->{start},$feat->{end});
         my ($starts_distance,$ends_distance) = (abs($feat_start-$start), abs($feat_end-$end));
	     
	 
	 # exactly the same
	 if($starts_distance == 0 and $ends_distance == 0 and $ies_sequence eq $feat_sequence) {
	    return ($feat,'identical','TRUE',0);	 
	 } else {
	    my $dist =min($starts_distance,$ends_distance);
	    
	    my $feat_boundary_relation = $self->_get_boundary_relation($ies,$feat);
	    # get the closest, or in the case where distance is the same : get the best one
	    if($min_dist eq '' or $dist < $min_dist 
	       or ($dist == $min_dist and $boundary_relation_scores{$feat_boundary_relation} > $boundary_relation_scores{$boundary_relation})) {
	       ($closest_feature,$min_dist,$boundary_relation) = ($feat,$dist,$feat_boundary_relation);
	    
	    # multiple mainly left and right
	    } elsif($dist == $min_dist and $boundary_relation_scores{$feat_boundary_relation} == $boundary_relation_scores{$boundary_relation}) {
	       
	       $self->stdlog("# WARN : $ies_id $boundary_relation AND $feat_boundary_relation\n");
	       $boundary_relation = join("_and_",sort($boundary_relation,$feat_boundary_relation));
	       $boundary_relation_scores{$boundary_relation} = $boundary_relation_scores{$feat_boundary_relation};
	       
	    }

         }
      }
   } elsif($all_features) {
      
      foreach my $feat (@{$all_features}) {
         my $dist = $self->_get_features_distance($ies,$feat);
	 my $feat_boundary_relation = $self->_get_boundary_relation($ies,$feat);	 
         if($min_dist eq '' or $dist < $min_dist
	    or ($dist == $min_dist and $boundary_relation_scores{$feat_boundary_relation} > $boundary_relation_scores{$boundary_relation})) {
            ($closest_feature,$min_dist,$boundary_relation) = ($feat,$dist,$self->_get_boundary_relation($ies,$feat));
         } 
      }
   }
   return ($closest_feature,$boundary_relation,'FALSE',$min_dist);  
}


# get distance between elements
sub _get_features_distance {
   my ($self,$ies,$indel) = @_;
   my ($ies_start,$ies_end,$feat_start,$feat_end) = ($ies->{start},$ies->{end},$indel->{start},$indel->{end});
   
   # perfect => negative value
   if($ies_start <= $feat_end and $ies_end >= $feat_start) {
      return -1;      
   # overlap 
   } elsif($ies_start <= $feat_end and $ies_end >= $feat_start) {
      return 0;
   } elsif($ies_end < $feat_start) {
      return abs( $feat_start - $ies_end);   
   }  elsif($ies_start > $feat_end) {
      return abs($ies_start - $feat_end);   
   }
   $self->stderr("# ERROR : _get_features_distance $ies_start,$ies_end,$feat_start,$feat_end \n");  
}



# get the relation type between two elements
sub _get_boundary_relation {
   my ($self,$ies,$feat) = @_;


   my ($start,$end) = ($ies->{start},$ies->{end});
   my ($ies_sequence) = @{$ies->{attributes}->{sequence}};
   
   my ($feat_sequence) = @{$feat->{attributes}->{sequence}};
   
   my ($feat_start,$feat_end) = ($feat->{start},$feat->{end});
   my ($starts_distance,$ends_distance) = (abs($feat_start-$start), abs($feat_end-$end));
  
   my $boundary_relation;      
   # exactly the same
   if(($starts_distance <= $self->{MAX_DIST} and $ends_distance <= $self->{MAX_DIST} and $ies_sequence eq $feat_sequence)
     or (PARTIES::Utils->compare_IES_seq($ies_sequence, $feat_sequence) and $starts_distance==$ends_distance)) {
      ($boundary_relation) = ('identical');
   # same position 
   } elsif($starts_distance <= $self->{MAX_DIST} and $ends_distance <= $self->{MAX_DIST}) {
      ($boundary_relation) = ('identical');
   # same left
   } elsif($starts_distance <= $self->{MAX_DIST}) {
       #($boundary_relation) = ('left');
       if($feat_end < $end) {
          ($boundary_relation) = ('left_external');
       } else {
          ($boundary_relation) = ('left_internal');
       }
   # same right
   } elsif($ends_distance <= $self->{MAX_DIST}) {
       #($boundary_relation) = ('right');
       if($feat_start < $start) {
          ($boundary_relation) = ('right_internal');
       } else {
          ($boundary_relation) = ('right_external');
       }     
   # internal
   } elsif($start < $feat_start and $feat_end < $end) {
       ($boundary_relation) = ('external');
   # external
   } elsif($feat_start < $start and $end < $feat_end) {
       ($boundary_relation) = ('internal');
   # simple overlap
   } elsif($start <= $feat_end and $end >= $feat_start ) {
       ($boundary_relation) = ('only_overlap');
   } else {
       ($boundary_relation) = ('no_overlap');  
   }

   return $boundary_relation;

}


1;


