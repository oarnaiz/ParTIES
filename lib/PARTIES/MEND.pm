
package PARTIES::MEND;
use strict;
use base 'PARTIES::Root';

use PARTIES::Config;
use Data::Dumper;
use File::Basename;
use FindBin qw($Bin);

=head1 NAME

 ParTIES MIRET module - Method of Ies RETention

=head1 AUTHORS

2021, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Mapping file of reads on the genome"
				},
			GERMLINE_BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Mapping file on the germline genome"
				},
			IES => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"IES file, usually MICA output"
				},
			GERMLINE_IES => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"IES file, with coordinates on the germline genome (i.e given by Insert)"
				},
			GERMLINE_GENOME => {
				MANDATORY=>1, DEFAULT=>[], TYPE=>'MULTIPLE', RANK => 1,
				DESCRIPTION=>"Reference genome with IES (i.e output of Insert)"
				},		
			EXCISION_ERRORS => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"IES excision error file"
				},
			MAX_MISMATCH => {
				MANDATORY=>0, DEFAULT=> 1, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Number of mismatch allowed in a read."
				},
			MIN_NON_AMBIGOUS_DISTANCE => {
				MANDATORY=>1, DEFAULT=> 15, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Minimum distance to be non ambigous mapped read"
				},		
			MAX_INSERT_SIZE => {
				MANDATORY=>0, DEFAULT=>1000, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Max sequencing insert size"
				},						
			MAX_DIST => {
				MANDATORY=>1, DEFAULT=>2, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Max distance to TA"
				},
		);



my %FILE_EXTENSIONS = ( 
			 );
 

=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::MIRET->new({....})
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
		
	  my $germline_ies = $self->{OUT_DIR}."/Insert/Insert.gff3";
	  $self->{GERMLINE_IES} = $germline_ies if(!$self->{GERMLINE_IES} and -e $germline_ies);
	  
		
     my $bam = $self->{OUT_DIR}."/Map/".basename($self->{OUT_DIR}).".".basename($self->{GENOME}).".BOWTIE.sorted.bam";
     $self->{BAM} = $bam if(!$self->{BAM} and -e $bam);
    
     my $germline_bam = $self->{OUT_DIR}."/Map/".basename($self->{OUT_DIR}).".".basename(${$self->{GERMLINE_GENOME}}[0]).".BOWTIE.sorted.bam";
     $self->{GERMLINE_BAM} = $germline_bam if(!$self->{GERMLINE_BAM} and -e $germline_bam);
     
     my $excision_error_file = $self->{OUT_DIR}."/Compare/Compare.current.gff3";
     $self->{EXCISION_ERRORS} = $excision_error_file if(!$self->{EXCISION_ERRORS} and -e $excision_error_file);
     
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


    # load IES excision error file
    $self->stderr("Read ".basename($self->{EXCISION_ERRORS})." ...\n" );
    $self->{LOADED_EXCISION_ERRORS} = PARTIES::Utils->read_gff_file($self->{EXCISION_ERRORS});
    $self->stderr("Done\n" ); 

    # Loading IES file with germline genome coordinates
    $self->stderr("Read ".basename($self->{GERMLINE_IES})." ...\n" );
    $self->{LOADED_GERMLINE_IES} = PARTIES::Utils->read_gff_file_by_id($self->{GERMLINE_IES});
    $self->stderr("Done\n" );

    # load GENOME file
    $self->stderr("Read ".basename($self->{GENOME})." ... " );
    $self->{LOADED_GENOME} = PARTIES::Utils->read_fasta_file($self->{GENOME});
    $self->stderr("Done\n" );

    $self->stderr("Read ".basename($self->{GERMLINE_GENOME}[0])." ... " );
    $self->{LOADED_GERMLINE_GENOME} = PARTIES::Utils->read_fasta_file($self->{GERMLINE_GENOME}[0]);
    $self->stderr("Done\n" );
        
    # bin directory
    $self->{BIN}=$Bin;

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
    my %results;
    return %results if(!$self->{LOADED_IES}->{$seq_id});
    
    my $min_non_ambigous_dist=$self->{MIN_NON_AMBIGOUS_DISTANCE};
    my $window_around_ta= $self->{MAX_DIST};
    my $max_insert_size=$self->{MAX_INSERT_SIZE};
    my $max_mismatch=$self->{MAX_MISMATCH};
    my $TELO_PATTERN = '(GGG(G|T)TT){3,}|(AA(A|C)CCC){3,}';

    my $mac_sam = new Bio::DB::Sam( -bam => $self->{BAM}, -fasta => $self->{GENOME}, -expand_flags => 'true');
    my $mac_ies_sam = new Bio::DB::Sam( -bam => $self->{GERMLINE_BAM}, -fasta => $self->{GERMLINE_GENOME}[0], -expand_flags => 'true');
    
    # error file
    my $err_file = $self->{PATH}."/tmp/$seq_id.err";
    open(ERR,">$err_file") or die "Can not open $err_file";
    print ERR join("\t",qw(RNAME REARRANGEMENT CIGAR SEQ PART_LENGTH)),"\n";
    
    my %excision_errors;
    my %unknown;
    # for each IESs
    foreach my $ies_mac (sort {$a->{start}<=>$b->{start}} @{$self->{LOADED_IES}->{$seq_id}}) {
        my ($ies_id)=@{$ies_mac->{attributes}->{ID}};

        my ($ta_start,$ta_end) = ($ies_mac->{start},$ies_mac->{end});
        # define ROI
        my ($mac_region_start,$mac_region_end) = ($ta_start-$max_insert_size,$ta_end+$max_insert_size);
        $mac_region_start=1 if($mac_region_start <1);
        $mac_region_end=$self->{LOADED_GENOME}->{$seq_id}->length if($mac_region_end > $self->{LOADED_GENOME}->{$seq_id}->length);
       
        
        my %status;
        
        
        my %putative_pcr_dup;
        # get reads by pairs on the MAC genome 
        foreach my $pair  ($mac_sam->features(-type   => 'read_pair',-seq_id => $seq_id, -start=>$mac_region_start, -end => $mac_region_end)) {
            my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
            next if(!$first_mate or !$second_mate);
      
            next if($first_mate->seq_id ne $second_mate->seq_id);
            
            my @pos = sort {$a<=>$b} ($first_mate->start,$first_mate->end,$second_mate->start,$second_mate->end);
            my ($insert_start,$insert_end) = ($pos[0],$pos[$#pos]);
            next if(($insert_end-$insert_start+1) > $max_insert_size);
            
            # uniquename of the paair
            my $pair_uniquename = join("__",$first_mate->seq_id,$first_mate->start,$first_mate->end,$first_mate->query->dna,
                                        $second_mate->seq_id,$second_mate->start,$second_mate->end,$second_mate->query->dna);
            
            my $qname = $first_mate->qname;
            my ($first_cigar,$second_cigar) = ($first_mate->cigar_str,$second_mate->cigar_str);
            
            if($first_cigar =~/^\d+M$/ and $second_cigar =~/^\d+M$/) {
                next if(_mismatch_count($first_mate) > $max_mismatch or _mismatch_count($second_mate) > $max_mismatch);
                # left mac_junction
                if($first_mate->end < ($ta_start-$window_around_ta) and $second_mate->start <= ($ta_start-$window_around_ta) and $second_mate->end > ($ta_end+$min_non_ambigous_dist)) {
                    if(!$putative_pcr_dup{$pair_uniquename}) {
                        $status{LEFT}->{mac_junction}->{$qname}=1;
                    } else {
                        $status{SKIP}->{mac_junction}->{$qname}=1;
                    }
                   
                    
                # right mac_junction
                } elsif($second_mate->start > ($ta_end+$window_around_ta) and $first_mate->end >= ($ta_end+$window_around_ta) and $first_mate->start < ($ta_start-$min_non_ambigous_dist)) {
                    if(!$putative_pcr_dup{$pair_uniquename}) {
                        $status{RIGHT}->{mac_junction}->{$qname}=1;
                    } else {
                        $status{SKIP}->{mac_junction}->{$qname}=1;
                    }
                    
                }
            # special case deletion through the gap
            } elsif( ($first_cigar =~/^\d+M\d+D\d+M$/ and $second_cigar =~/^\d+M$/) or ($second_cigar =~/^\d+M\d+D\d+M$/ and $first_cigar =~/^\d+M$/) ) {
                
                my ($del_mate,$fm_mate,$nb_del,$del_pos);
                if($first_cigar =~/^(\d+)M(\d+)D\d+M$/) {
                    $del_pos= $first_mate->start+$1;
                    $nb_del = $2;
                    ($del_mate,$fm_mate)=($first_mate,$second_mate);
                } elsif($second_cigar =~/^(\d+)M(\d+)D\d+M$/) {
                    $del_pos= $second_mate->start+$1;
                    $nb_del = $2;
                    ($del_mate,$fm_mate)=($second_mate,$first_mate);
                }
                
                next if(_mismatch_count($del_mate) > $nb_del or _mismatch_count($fm_mate) > $max_mismatch);
                
                my $dist_ta = (abs($ta_start-$del_pos) < abs($ta_end-$del_pos)) ? abs($ta_start-$del_pos) : abs($ta_end-$del_pos);
                
                # left
                if($fm_mate->end < ($ta_start-$window_around_ta) and $del_mate->start <= ($ta_start-$window_around_ta) and $del_mate->end > ($ta_end+$min_non_ambigous_dist) and $dist_ta <=$window_around_ta) {
                    if(!$putative_pcr_dup{$pair_uniquename}) {
                        $status{LEFT}->{mac_junction_with_del}->{$qname}=1;
                    } else {
                        $status{SKIP}->{mac_junction_with_del}->{$qname}=1;
                    }
                # right
                } elsif($fm_mate->start > ($ta_end+$window_around_ta) and $del_mate->start <= ($ta_start-$window_around_ta) and $del_mate->end > ($ta_end+$min_non_ambigous_dist) and $dist_ta <=$window_around_ta) {
                                
                    if(!$putative_pcr_dup{$pair_uniquename}) {
                        $status{RIGHT}->{mac_junction_with_del}->{$qname}=1;
                    } else {
                        $status{SKIP}->{mac_junction_with_del}->{$qname}=1;
                    }
                }
                
                
            }
            
            # PCR duplicates
            $putative_pcr_dup{$pair_uniquename}=1;
            
            
        }

        
        # Check MAC+IES genome 
        my $ies = $self->{LOADED_GERMLINE_IES}->{$ies_id}->[0];
        my ($mac_ies_seq_id,$ies_start,$ies_end)=($ies->{seq_id},$ies->{start},$ies->{end}+2);
        
        %excision_errors = $self->_load_excision_errors($mac_ies_seq_id) if(!%excision_errors);
        
        my ($region_start,$region_end) = ($ies_start-$max_insert_size,$ies_end+$max_insert_size);
        
        $region_start=1 if($region_start <1);
        $region_end=$self->{LOADED_GERMLINE_GENOME}->{$mac_ies_seq_id}->length if($region_end > $self->{LOADED_GERMLINE_GENOME}->{$mac_ies_seq_id}->length);
        foreach my $pair  ($mac_ies_sam->features(-type   => 'read_pair',-seq_id => $mac_ies_seq_id, -start=>$region_start, -end => $region_end)) {
            my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
            next if(!$first_mate or !$second_mate);
            my $qname = $first_mate->qname;
            next if(_has_status($qname,\%status));
            
            next if(_mismatch_count($first_mate) > $max_mismatch or _mismatch_count($second_mate) > $max_mismatch);
                
            my ($first_cigar,$second_cigar) = ($first_mate->cigar_str,$second_mate->cigar_str);
            next if(!_valid_cigar($first_cigar) or  !_valid_cigar($second_cigar));
            next if($first_mate->seq_id ne $second_mate->seq_id);
            my @pos = sort {$a<=>$b} ($first_mate->start,$first_mate->end,$second_mate->start,$second_mate->end);
            my ($insert_start,$insert_end) = ($pos[0],$pos[$#pos]);
            next if(($insert_end-$insert_start+1) > $max_insert_size);
            
            my $pair_uniquename = join("__",$first_mate->seq_id,$first_mate->start,$first_mate->end,$first_mate->query->dna,
                                        $second_mate->seq_id,$second_mate->start,$second_mate->end,$second_mate->query->dna);
            # next if PCR duplicates
            next if($putative_pcr_dup{$pair_uniquename});
            $putative_pcr_dup{$pair_uniquename}=1;
            
            if($insert_start <=$ies_end and $insert_end>=$ies_start) {
                
                    
                # left  boundary
                if($first_mate->end <= ($ies_start-$window_around_ta) and  $second_mate->start <=($ies_start-$window_around_ta)) {
                    my $boundary='LEFT';
                    if($second_cigar=~/^\d+M$/ ) {
                        # extremeties if
                        if($second_mate->end >= ($ies_start+$min_non_ambigous_dist)) {
                            $status{$boundary}->{ies}->{$qname}=1;
                        } elsif($second_mate->end >=  ($ies_start-$window_around_ta) and $second_mate->end <=  ($ies_start+$window_around_ta) ) {
                            $status{$boundary}->{extremities}->{$qname}=1;
                        } 
                    } elsif($first_cigar=~/^\d+M$/ and $second_cigar=~/^\d+M(\d+)S$/ and ($second_mate->end >=  ($ies_start-$window_around_ta) and $second_mate->end <=  ($ies_start+$window_around_ta)) ) {
                        my $space_length=$1;
                        
                        if($space_length >= $min_non_ambigous_dist) {
                            if($excision_errors{$qname}) {
                                 $status{$boundary}->{illegitimate_rearrangement}->{$qname}=1;
                                
                            } else {
                                if(_unmapped_seq($second_mate) =~/$TELO_PATTERN/) {
                                    $status{$boundary}->{telomerization}->{$qname}=1;
                                } else {
                                    $unknown{$ies_id}->{$boundary}->{$qname}={match=>$second_mate,unmapped_seq=>_unmapped_seq($second_mate)};
                                }
                                
                            
                            }
                        }
                    }
                    
                    
                # right boundary
                } elsif($second_mate->start >= ($ies_end+$window_around_ta) and $first_mate->end >= ($ies_end+$window_around_ta)) {
                    my $boundary='RIGHT';
                    if($first_cigar=~/^\d+M$/) {
                        if($first_mate->start <= ($ies_end-$min_non_ambigous_dist)  ) {
                            $status{$boundary}->{ies}->{$qname}=1;
                        } elsif($first_mate->start >=  ($ies_end-$window_around_ta) and $first_mate->start <=  ($ies_end+$window_around_ta) ) {
                            $status{$boundary}->{extremities}->{$qname}=1;
                        } 
                       
                        
                    } elsif($second_cigar=~/^\d+M$/ and $first_cigar=~/^(\d+)S\d+M$/ and  ($first_mate->start >=  ($ies_end-$window_around_ta) and $first_mate->start <=  ($ies_end+$window_around_ta)) ) {
                        my $space_length=$1;
                        
                        if($space_length >= $min_non_ambigous_dist) {
                            if($excision_errors{$qname}) {
                                $status{$boundary}->{illegitimate_rearrangement}->{$qname}=1;
                            } else {
                                if(_unmapped_seq($first_mate) =~/$TELO_PATTERN/) {
                                    $status{$boundary}->{telomerization}->{$qname}=1;
                                } else {
                                    #$status{$boundary}->{unknown}->{$qname}=1;
                                    $unknown{$ies_id}->{$boundary}->{$qname}={match=>$first_mate,unmapped_seq=>_unmapped_seq($first_mate)};
                                }
                            
                            }
                        }
                    }
                }
            }
               
        }

        
        foreach my $boundary (qw(LEFT RIGHT)) {
            foreach my $state (qw(mac_junction mac_junction_with_del ies extremities illegitimate_rearrangement telomerization unknown)) {
                my $nb= ($status{$boundary}->{$state}) ? scalar keys %{$status{$boundary}->{$state}}: 0;
                $results{$self->get_mode()}->{$seq_id}->{$ies_id}->{$boundary}->{$state}=$nb;
            }
            
            
        }
        
        
    }
    
    if(%unknown) {
        my ($bowtie2_build,$bowtie2,$samtools) = (PARTIES::Config->get_program_path('bowtie2-build'),PARTIES::Config->get_program_path('bowtie2'),PARTIES::Config->get_program_path('samtools'));
       
        # write fasta for unmapped
        my $unknown_fasta_file = $self->{PATH}."/tmp/".$seq_id."_tolocate.fa";
        open(FA,">$unknown_fasta_file") or die $unknown_fasta_file;
        foreach my $ies_id (keys %unknown) {
            foreach my $boundary (keys %{$unknown{$ies_id}}) {
                foreach my $qname (keys %{$unknown{$ies_id}->{$boundary}}) {
                    my $name = join("__",$qname,$ies_id,$boundary);
                    print FA ">$name\n",$unknown{$ies_id}->{$boundary}->{$qname}->{unmapped_seq},"\n";
                }
            }
        }
        close FA;
        

        # mapping
        my $unmapped_file_bam = $self->{PATH}."/tmp/unmapped_$seq_id\_reads_vs_ref.BOWTIE.sorted"; 
        my $target_ref = $self->{GERMLINE_GENOME}[0];
        system("$bowtie2 --threads 1 --end-to-end --quiet -f  --very-sensitive -k 2 -x ".$target_ref." -U $unknown_fasta_file | $samtools view -F 4 -uS - 2> /dev/null | $samtools sort -o $unmapped_file_bam.bam - > /dev/null 2>&1");
        system("$samtools index $unmapped_file_bam.bam");
        my %remapped;
        foreach my $line (`$samtools view $unmapped_file_bam.bam`) {
            chomp $line;
            if($line=~/NM:i:(\d)/) {
                if($1 <=$max_mismatch) {
                    my ($rname,$flag,$target_seq_id,$pos,$mapq,$cigar) = split /\t/,$line;
                    if($cigar=~/^(\d+)M$/) {
                        my ($target_start,$target_end) = ($pos,$pos+$1);
                        my ($qname,$ies_id,$boundary) = split /__/,$rname;
                        $remapped{$ies_id}->{$boundary}->{$qname}={seq_id=>$target_seq_id,start=>$target_start,end=>$target_end};
                    }
                }
            }
            
            
        }
        
        foreach my $ies_id (keys %unknown) {
            my $ies = $self->{LOADED_GERMLINE_IES}->{$ies_id}->[0];
            my ($ies_start,$ies_end)=($ies->{start},$ies->{end}+2);
        
            foreach my $boundary (keys %{$unknown{$ies_id}}) {
                foreach my $qname (keys %{$unknown{$ies_id}->{$boundary}}) {
                    if($remapped{$ies_id}->{$boundary}->{$qname}) {
                        my ($target_start,$target_end) = ($remapped{$ies_id}->{$boundary}->{$qname}->{start},$remapped{$ies_id}->{$boundary}->{$qname}->{end});
                        my ($match_start,$match_end) = ($unknown{$ies_id}->{$boundary}->{$qname}->{match}->start,$unknown{$ies_id}->{$boundary}->{$qname}->{match}->end);
                        
                        if($remapped{$ies_id}->{$boundary}->{$qname}->{seq_id} eq $unknown{$ies_id}->{$boundary}->{$qname}->{match}->seq_id and $target_start <= $match_end and $target_end >= $match_start) {
                            $self->stderr( "OVERLAP : $target_start,$target_end <=> $match_start,$match_end\n");
                        } else {
                            
                            
                            # missing MAC junction
                            if($remapped{$ies_id}->{$boundary}->{$qname}->{seq_id} eq $unknown{$ies_id}->{$boundary}->{$qname}->{match}->seq_id
                                and (
                                    ($boundary eq 'LEFT' and abs($target_start-$ies_end) <= $window_around_ta)
                                    or ($boundary eq 'RIGHT' and abs($target_end-$ies_start) <= $window_around_ta))) {
                                    
                               
                                $results{$self->get_mode()}->{$seq_id}->{$ies_id}->{$boundary}->{mac_junction}++; 
                            } else {
                                
                                $results{$self->get_mode()}->{$seq_id}->{$ies_id}->{$boundary}->{illegitimate_rearrangement}++; 
                            }
                            
                            
                            
                            
                            
                        }
                    } else {
                        $results{$self->get_mode()}->{$seq_id}->{$ies_id}->{$boundary}->{unknown}++;
                    }
                }
            }
        }
        system("rm -f $unmapped_file_bam.bam $unmapped_file_bam.bam.bai $unknown_fasta_file");
        
    }
    close ERR;
    


    return %results;
}



=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
           Details :
 	      - compare the current retention scores to a control dataset if provided (CONTROL)
 	      - calls R source to compute statistical tests and returns significance (T/F)
 	      - write the final GFF file

 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
	my ($self, $results)=@_;
	$results=$results->{$self->get_mode()};
    
    
    ## Writes the final GFF file
	my $tab_file = $self->{PATH}."/".uc($self->get_mode).".tab";
  	open(TAB,">$tab_file") or die "Can not open $tab_file";   
	print TAB "## ",join(" ",$self->_command_line),"\n";
    
    
    my @boundaries=qw(LEFT RIGHT);
    my @states=qw(mac_junction mac_junction_with_del ies extremities illegitimate_rearrangement telomerization unknown);
    
    my @header=qw(ID);
    foreach my $boundary (@boundaries) {
        foreach my $state (@states) {
            push @header,$boundary."_".$state;
        }
    }
    print TAB join("\t",@header),"\n";
    
    foreach my $seq_id (sort keys %{$results}){
        foreach my $ies_id (sort keys %{$results->{$seq_id}}) {
            my @line =($ies_id);
            foreach my $boundary (@boundaries) {
                foreach my $state (@states) {
                    my $n = ($results->{$seq_id}->{$ies_id}->{$boundary}->{$state}) ? $results->{$seq_id}->{$ies_id}->{$boundary}->{$state} : 0;
                    push @line,$n;
                }
            }
            print TAB join("\t",@line),"\n";
            
        }
        
    }
    close TAB;
    
	$self->remove_temporary_files;

}

sub _has_status {
    my ($qname,$status) = @_;
    return 0 if(!$status);
    foreach my $boundary (keys %$status) {
        foreach my $st (keys %{$status->{$boundary}}) {
            return 1 if($status->{$boundary}->{$st}->{$qname});
        }
    }
    return 0;
}  
sub _valid_cigar {
    my ($cigar)=@_;
    if($cigar=~/^\d+M$/ or $cigar=~/^\d+M\d+S$/ or $cigar=~/^\d+S\d+M$/) {
        return 1;
    }
    return 0;
}

sub _mismatch_count {
    my ($aln)=@_;
    my $nm=$aln->get_tag_values("NM");
    $nm=0 if(!$nm);
    return $nm;
}
sub _unmapped_seq {
    my ($aln)=@_;
    my $cigar = $aln->cigar_str;
    my $unmap_seq='';
    if($cigar=~/^\d+M(\d+)S$/) {
        $unmap_seq= substr($aln->query->dna,length($aln->query->dna)-$1);
        die if(length($unmap_seq)!=$1);
    } elsif($cigar=~/^(\d+)S\d+M$/) {
        $unmap_seq= substr($aln->query->dna,0,$1);
        die if(length($unmap_seq)!=$1);
        
        
    }
    
    
    return $unmap_seq;
}



sub _load_excision_errors {
    my ($self,$seq_id)=@_;
    my %errors;
    if($self->{LOADED_EXCISION_ERRORS}->{$seq_id}) {
        foreach my $error (@{$self->{LOADED_EXCISION_ERRORS}->{$seq_id}}) {
            my ($boundary_relation)=@{$error->{attributes}->{boundary_relation}};
            next if($boundary_relation eq 'identical');
            foreach my $read_name (@{$error->{attributes}->{read_names}}) {
                $errors{$read_name}=$boundary_relation;
            }
        }
    }
    
    return %errors;
}


1;

