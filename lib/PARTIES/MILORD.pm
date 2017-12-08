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
				MANDATORY=>1, DEFAULT=>1e4, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Maximum size for a deletion to be reported (The maximum is 1e4, high values slows the calculation)"
				},
			JUNCTION_FLANK_SEQ_LENGTH => {
				MANDATORY=>1, DEFAULT=>15, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Length of the flanking sequence to report, arround the deletion"
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
			REPORT_READ_NAMES => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Should read names be reported when showing a deletion"
				},	

		);



my %FILE_EXTENSIONS = ( 
			gff3 => {DESC=> 'GFF3 file', EXT => 'gff3' },		
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
  my $file2locate=$self->{PATH}."/tmp/".$seq_id."_tolocate.fa";
  my $current_seq= $seq->seq;
  
  
########
# Step 1 : Gathering partially mapped reads
#			=> Parse the BAM file to detect partially mapped reads
#			=> Hash with partially mapped reads
#			=> File with unmapped sequences
#	$self->stderr("Gathering partially mapped reads on $seq_id\n");
   my %reads=$self->_gather_partially_mapped(\$sam, $seq_id, $file2locate);
   $self->{$seq_id}->{PARTIAL_MAP}+=scalar(keys %reads);
   return 1 if(scalar(keys %reads)==0);

	

########
# Step 2 : Retrieving unmapped sequence localization
#			=> Create single scaffold reference
#			=> Index the reference (BOWTIE2)
#			=> Map the sequences on the new ref (BOWTIE2)
#	$self->stderr("Remapping reads on $seq_id\n");
   my ($new_mapping, $new_ref)=$self->_bowtie2_mapping($seq_id, \$current_seq, $file2locate);


########
# Step 3 : Reconstituting the read structure
#			=> Parse the BOWTIE2 mapping
#			=> Check coherence between initially mapped read part and newly mapped part
#			=> Determine the theoric difference between the read and the reference
#	$self->stderr("Read reconstitution on $seq_id\n");
   $self->_reconstitute_read_structure(\%reads, $new_mapping, $new_ref, $seq_id, \$current_seq);

########
# Step 4 : Simplifiing Complexity
#			=> Parse the Hash of reads to look for similar coordinates
#			=> Merge similar observations
#	$self->stderr("Building segments from $seq_id\n");
   my %segments;
   $self->_get_segment_from_reads(\%reads, \%{$segments{$seq_id}}, $seq_id, \$current_seq);

   $self->{$seq_id}->{COHERENT_SEGMENTS}+=scalar(keys %{$segments{$seq_id}});

########
# Step 5 : Detecting TA bounded deletions
#			=> Run muscle alignment
#			=> Use the get_IES_from_aln
#   $self->stderr("Refining segments on $seq_id\n");
    $self->_refine_segments(\%{$segments{$seq_id}}, $seq_id, \$current_seq, $sam);

#   $self->stderr("Merging results on $seq_id\n");
    my %results;
    foreach my $i (sort {$segments{$seq_id}->{$a}->{start}<=>$segments{$seq_id}->{$b}->{start}} keys %{$segments{$seq_id}}){
    	    next if($segments{$seq_id}->{$i}->{tobereported}==0);
    	    my $indel_id=PARTIES::Utils->generate_IES_id($self->{PREFIX},$seq_id,$segments{$seq_id}->{$i}->{start},$segments{$seq_id}->{$i}->{end});
    	    if(defined($results{$indel_id})){
    		    $results{$indel_id}->{support_ref}+=$segments{$seq_id}->{$i}->{support_ref};
    		    $results{$indel_id}->{support_variant}+=$segments{$seq_id}->{$i}->{support_variant};
    		    $results{$indel_id}->{read_names}=join(',',$results{$indel_id}->{read_names}, @{$segments{$seq_id}->{$i}->{read_names}});
    	    }
    	    else{
    		    $results{$indel_id}={
    			    start => $segments{$seq_id}->{$i}->{start},
    			    end => $segments{$seq_id}->{$i}->{end},
    			    ID => $indel_id,
    			    Name => $indel_id,
    			    junction_seq => $segments{$seq_id}->{$i}->{junction_seq},
    			    sequence => $segments{$seq_id}->{$i}->{sequence},
    			    bounded_by_ta => $segments{$seq_id}->{$i}->{bounded_by_ta},
    			    support_ref => $segments{$seq_id}->{$i}->{support_ref}, 
    			    support_variant => $segments{$seq_id}->{$i}->{support_variant},
    		    };
    		    $results{$indel_id}->{read_names} = join(',', @{$segments{$seq_id}->{$i}->{read_names}}); 
    	    }
    }
    $self->stdlog(join("\t",$seq_id, $self->{$seq_id}->{PARTIAL_MAP}, $self->{$seq_id}->{COHERENT_SEGMENTS}, scalar(keys %results))."\n");
#   $self->stderr("Writting results from on $seq_id\n");

    ## Writes the final GFF file
    my $gff_file = $self->{PATH}."/tmp/".$seq_id.".gff3";
    open(GFF,">$gff_file") or die "Can not open $gff_file";
    foreach my $indel (sort {$results{$a}->{start}<=>$results{$b}->{start}} keys %results){
    	    print GFF join("\t", $seq_id, 'MILORD', 'internal_eliminated_sequence',
    			    $results{$indel}->{start}, 
    			    $results{$indel}->{end},
    			    '.', '.', '.', '');
    	    foreach my $att (sort keys %{$results{$indel}}){
    		    next if($att eq "read_names" &&  $self->{REPORT_READ_NAMES} eq 'FALSE');
    		    print GFF $att."=".$results{$indel}->{$att}.";";
    	    }
    	    print GFF "\n";
    }
    close(GFF);
    my %toreturn;
    $toreturn{$self->get_mode()}->{$seq_id}=\%results;
    return %toreturn;
}	

=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
	my ($self, $results)=@_;
	$results=$results->{$self->get_mode()};
	my %stats;
 	foreach my $seq_id (keys %{$results}){
		foreach my $i (keys %{$results->{$seq_id}}){
			my $type=$results->{$seq_id}->{$i}->{type};
			$stats{$type}->{obs}+=$results->{$seq_id}->{$i}->{support_variant};
			$stats{$type}->{nb}++;
		}
		$stats{partial_map}->{counts}+=$self->{PARTIAL_MAP}->{$seq_id};
		$stats{coherent_segments}->{counts}+=$self->{COHERENT_SEGMENTS}->{$seq_id};
	}
	foreach my $k (keys %stats){
		foreach my $st (keys %{$stats{$k}}){
			$self->stdlog(" $k - $st = ".$stats{$k}->{$st}."\n");
		}
	}
	$self->SUPER::finish;
}




sub _gather_partially_mapped{
	my ($self, $sam, $seq_id, $file2locate)= @_;
	my %reads;
	open(SEQ, ">$file2locate") or die $file2locate;
	my ($skipped_reads, $kept_reads, $mdm, $ms, $sm)=(0,0,0,0,0);
	my @align=$$sam->get_features_by_location(-seq_id=>$seq_id);
	foreach my $aln (@align){
		my ($type, $unmatch, $match, $seq2analyse)=check_mapping($aln);
		$skipped_reads++ if($type eq "-1");
		next if($type eq "-1" || $unmatch < $self->{MIN_SIZE});
		$kept_reads++;
		my ($size, $start, $end, $seq)=(0,0,0,'NA');
		my $rname=$aln->query->name;
		$start=$aln->start;
			
		if($type eq 'MDM'){
			$mdm++;
			$size=$unmatch;
			$reads{$rname}={
				read_seq_id => $seq_id, 
				read_name => $rname,
				read_type => "MDM",
				cigar => $aln->cigar_str,
				remapped_status => "1",
				indel_index => $match,
				indel_index_desc => abs((length($seq2analyse)/2)-$match),
				indel_start => $aln->start+$match-1,
				indel_end => $aln->start+$match+$size-2,
				indel_size => $size,
				indel_seq => substr($aln->dna,$match, $size),
				read_seq => $seq2analyse,
				f_seq => substr($seq2analyse, 0, $match).'DELETION'.substr($seq2analyse, $match),
				read_start => $aln->start,
				read_end => $aln->end,
				read_length => length($seq2analyse),
				tabound => "none",
				read_cigar => $type,
				read_strand => $aln->strand 
			};
		}
		else{
			$reads{$rname}={
				read_seq_id => $seq_id,		
				read_name => $rname,
				read_type => "UN",		
				cigar => $aln->cigar_str,
				remapped_status => "not_retrieved",
				read_strand => $aln->strand,		
				indel_start => $start-1,		
				read_start => $start,
				read_end => $aln->end,		
				read_length => length($seq2analyse),
				read_seq => $aln->query->dna,
				unmatch => $unmatch,
				read_cigar => $type,
				mate_start => $aln->mate_start,
				mate_end => $aln->mate_end,
				tabound => "none"
			};
			my $seq2locate;
			if($type eq 'MS'){
				$ms++;
				$seq2locate=substr($aln->query->dna, $match);
				$reads{$rname}->{read_type}="MS";
				$reads{$rname}->{indel_index}=$match+1;
				$reads{$rname}->{indel_index_desc}=abs((length($seq2analyse)/2)-($match+1));
				$reads{$rname}->{indel_start}=$start+$match;
				$reads{$rname}->{seq2locate}=$seq2locate;
				$reads{$rname}->{f_seq}=substr($aln->query->dna, 0, $match).'DELETION'.substr($aln->query->dna, $match);
			}
			elsif($type eq 'SM'){
				$sm++;
				$seq2locate=substr($aln->query->dna,0, $unmatch);
				$reads{$rname}->{read_type}="SM";
				$reads{$rname}->{indel_index}= $unmatch+1;
				$reads{$rname}->{indel_index_desc}=abs((length($seq2analyse)/2)-($unmatch+1));
				$reads{$rname}->{indel_end}=$start+1;#+$match;
				$reads{$rname}->{seq2locate}=$seq2locate;
				$reads{$rname}->{f_seq} = substr($aln->query->dna, 0, $unmatch).'DELETION'.substr($aln->query->dna, $unmatch);
			}
			$seq2locate=~s/(.{60})/$1\n/g;
			print SEQ ">$rname\n$seq2locate\n";
		}	
	}
	close(SEQ);
	return %reads ;
}



sub _bowtie2_mapping{
	my ($self, $seq_id, $ref_seq, $file2locate)=@_;
	my $file_ref=$self->{PATH}."/tmp/".$seq_id."_as_ref.fa";
	my $file_bam = $self->{PATH}."/tmp/seq2.BT2.".$seq_id."_as_ref";
	
	my ($bowtie2_build,$bowtie2,$samtools) = (PARTIES::Config->get_program_path('bowtie2-build'),PARTIES::Config->get_program_path('bowtie2'),PARTIES::Config->get_program_path('samtools'));
	if(-s $file2locate){ 
		my $reformat_seq=$$ref_seq;
		$reformat_seq=~s/(.{60})/$1\n/g;
		open(REF, ">$file_ref");print REF ">$seq_id\n$reformat_seq\n";close(REF);
		system("$bowtie2_build $file_ref $file_ref > /dev/null 2>&1");
		system("$bowtie2 --threads 1 --end-to-end --quiet -f  --very-sensitive -k 2 -x $file_ref -U $file2locate | $samtools view -F 4 -uS - 2> /dev/null | $samtools sort - $file_bam  > /dev/null 2>&1");
		system("$samtools index $file_bam.bam");
	}else{ system("touch $file_bam.bam");}

	return("$file_bam.bam",$file_ref);
}	


	
sub _reconstitute_read_structure{
	my ($self, $reads, $mapping, $ref, $seq_id, $seq_id_seq)=@_;

	if(-s $mapping){
		# Parsing the mapping.
		my $bam = Bio::DB::Sam->new(-bam  =>$mapping,-fasta=>$ref)  ;
		my %remapped;
		my @align = $bam->get_features_by_location(-seq_id => $seq_id);
		foreach my $aln (@align){
			my $name=$aln->query->name;
			my $cur_read=$reads->{$name};
			$remapped{$name}->{remapped_status}='multiple_hit' if(defined($remapped{$name}));
			$remapped{$name}->{remapped_status}='strand_insconsistencies' if( $aln->strand != 1);
			$remapped{$name}->{remapped_status}='mismatch' if($aln->get_tag_values("NM") != 0);
			$remapped{$name}->{name}=$name;
			next if(defined($remapped{$name}->{remapped_status})); 
			$remapped{$name}->{remapped_status}='1';
			$remapped{$name}->{remapped_strand}=$aln->strand;
			$remapped{$name}->{remapped_start}=$aln->start;
			$remapped{$name}->{remapped_end}=$aln->end;
		}
		foreach my $name (keys %remapped){
			my $cur_read=$reads->{$name};
			$cur_read->{remapped_status}=$remapped{$name}->{remapped_status};
			next if($remapped{$name}->{remapped_status} ne '1');
			$cur_read->{remapped_strand}=$remapped{$name}->{remapped_strand};
			$cur_read->{remapped_start}=$remapped{$name}->{remapped_start};
			$cur_read->{remapped_end}=$remapped{$name}->{remapped_end};
			

			#### CHECK COHERENCE WITH ORIGINALLY MAPPED READS ###;
			my $coherent=1;
			if($cur_read->{read_type} eq 'MS'){
				if($cur_read->{read_strand}==1){
					if($cur_read->{indel_start} >= $cur_read->{remapped_start} || $cur_read->{indel_start} >= $cur_read->{remapped_end}){
						$cur_read->{remapped_status}="position_inconsistenties";						
						$coherent=0;
					}
					if($cur_read->{mate_end} <= $cur_read->{remapped_start}){
						$cur_read->{remapped_status}="mate_inconsistenties";
						$coherent=0;
					}
					if($coherent==1){
						$cur_read->{read_end}=$cur_read->{remapped_end}; # La position de fin de read devient celle du remapped
						$cur_read->{indel_end}=$cur_read->{remapped_start}-1; # La position de fin de segment devient la position de start du remapped
						$cur_read->{indel_seq}=substr($$seq_id_seq, $cur_read->{indel_start}-1, $cur_read->{indel_end}-$cur_read->{indel_start}+1);
						$cur_read->{indel_size}=length($cur_read->{indel_seq});
						$cur_read->{remapped_status}="1";
					}
				}
				elsif($cur_read->{read_strand}==-1){
					if($cur_read->{indel_start} >= $cur_read->{remapped_start} || $cur_read->{indel_start} >= $cur_read->{remapped_end}){
						$cur_read->{remapped_status}="position_inconsistenties";						
						$coherent=0;
					}
					if($coherent==1){
						$cur_read->{read_end}=$cur_read->{remapped_end}; # La position de fin de read devient celle du remapped
						$cur_read->{indel_end}=$cur_read->{remapped_start}-1;# La position de fin de segment devient la position de start du remapped -1
						$cur_read->{indel_seq}=substr($$seq_id_seq, $cur_read->{indel_start}-1, $cur_read->{indel_end}-$cur_read->{indel_start}+1);
						$cur_read->{indel_size}=length($cur_read->{indel_seq});
						$cur_read->{remapped_status}="1";
					}
				}
			}

			elsif($cur_read->{read_type}  eq 'SM'){
				if($cur_read->{read_strand}==1)	 {
					if($cur_read->{indel_start} <= $cur_read->{remapped_start} || $cur_read->{indel_start} <= $cur_read->{remapped_end}){
						$cur_read->{remapped_status}="position_inconsistenties";						
						$coherent=0;
					}
					if($coherent==1){
						$cur_read->{indel_end}=$cur_read->{read_start}-1; # La position de fin de segment devient la position de start du read initial
						$cur_read->{read_start}=$cur_read->{remapped_start}; # La position de start du read est celle du remapped
						$cur_read->{indel_start}=$cur_read->{remapped_end}+1; # la position de start du segment devient celle de la fin de remapped +1
						$cur_read->{indel_seq}=substr($$seq_id_seq, $cur_read->{indel_start}-1, $cur_read->{indel_end}-$cur_read->{indel_start}+1);
						$cur_read->{indel_size}=length($cur_read->{indel_seq});
						$cur_read->{remapped_status}="1";
					}
				}
				elsif($cur_read->{read_strand}==-1){
					if($cur_read->{indel_start} <= $cur_read->{remapped_start} || $cur_read->{indel_start} <= $cur_read->{remapped_end}){
						$cur_read->{remapped_status}="position_inconsistenties";						
						$coherent=0;
					}
					if($cur_read->{mate_start} >= $cur_read->{remapped_start}){
						$cur_read->{remapped_status}="mate_inconsistenties";
						$coherent=0;
					}
					if($coherent==1){
						$cur_read->{indel_end}=$cur_read->{read_start}-1; # La position de fin de segment devient la position de start du read initial
						$cur_read->{read_start}=$cur_read->{remapped_start}; # La position de start du read est celle du remapped
						$cur_read->{indel_start}=$cur_read->{remapped_end}+1; # la position de start du segment devient celle de la fin de remapped +1
						$cur_read->{indel_seq}=substr($$seq_id_seq, $cur_read->{indel_start}-1, $cur_read->{indel_end}-$cur_read->{indel_start}+1);
						$cur_read->{indel_size}=length($cur_read->{indel_seq});
						$cur_read->{remapped_status}="1";
					}
				}
			}
		}
	}
}

sub _get_segment_from_reads{
	my ($self, $reads, $segment, $seq_id, $current_seq) = @_;
	foreach my $rname (keys %{$reads}){
		next if($reads->{$rname}->{remapped_status} ne '1');
		next if($reads->{$rname}->{indel_size} < $self->{MIN_SIZE} );
		next if($reads->{$rname}->{indel_size} >= $self->{MAX_SIZE} );
		my $seg_id=PARTIES::Utils->generate_IES_id($self->{PREFIX},$seq_id,$reads->{$rname}->{indel_start},$reads->{$rname}->{indel_end});
		if(!defined($segment->{$seg_id}) || $segment->{$seg_id}->{indel_index_desc} gt $reads->{$rname}->{indel_index_desc}){
			$segment->{$seg_id}={
				seq_id => $seq_id,
				id => $seg_id,
				start => $reads->{$rname}->{indel_start},
				end => $reads->{$rname}->{indel_end},
				size => $reads->{$rname}->{indel_size},
				sequence => uc($reads->{$rname}->{indel_seq}),
				bounded_by_ta => 'FALSE', 
				indel_index => $reads->{$rname}->{indel_index},
				indel_index_desc => $reads->{$rname}->{indel_index_desc},
				read_seq => $reads->{$rname}->{read_seq},
				read_name => $rname,
				read_start => $reads->{$rname}->{read_start},
				read_end => $reads->{$rname}->{read_end},
				read_length => $reads->{$rname}->{read_length},
				cigar => $reads->{$rname}->{cigar},
				f_seq => $reads->{$rname}->{f_seq},
				ref_seq =>  uc(substr($$current_seq, $reads->{$rname}->{read_start}-1, $reads->{$rname}->{read_end} - $reads->{$rname}->{read_start}+1))
			};
			
		} 
		$segment->{$seg_id}->{support_variant}+=1;
		push @{$segment->{$seg_id}->{read_names}}, $rname;
	}
}

sub _refine_segments{
	my ($self, $segments, $seq_id, $current_seq, $sam) = @_;
	my $a=0;my $ndif=0;my $non_ies=0;
	foreach my $seg_id (keys %{$segments}){
		my $tname = $segments->{$seg_id}->{read_name};
		my $seq = Bio::Seq->new(-seq => $segments->{$seg_id}->{read_seq}, -id => $segments->{$seg_id}->{read_name});		
	   my $qname = $seq_id;
		my $contigseq = $$current_seq;

		my $taln=$segments->{$seg_id}->{f_seq};
		my $qaln=$segments->{$seg_id}->{ref_seq};
		my $del_char="-" x $segments->{$seg_id}->{size};
		$taln=~/^(\w+)DELETION(\w+)$/;
		my ($taln_left_flank, $taln_right_flank)=($1, $2);
		
		my ($mv_lft, $i)=(0,1);
		while(substr($taln_left_flank, length($taln_left_flank)-$i,1) eq 
				substr($segments->{$seg_id}->{sequence}, length($segments->{$seg_id}->{sequence})-$i,1)
				&& $i<=length($taln_left_flank) ){
				$mv_lft++;
				$i++;
		}
		$taln_left_flank=~/(\w+)(\w{$mv_lft})$/;
		$taln_left_flank=$1; my $to_move=$2;

		$taln=$taln_left_flank."DELETION".$to_move.$taln_right_flank;
		$taln=~s/DELETION/$del_char/;

		my $qstart_tmp = $segments->{$seg_id}->{start};
		my $tstart_tmp = $segments->{$seg_id}->{indel_index};
		
#   	if($qaln =~/N/){$self->stdlog("$seg_id\n");next;}
   	next if($qaln =~/N/);
		die "No seq align seq1=$taln seq2=$qaln" if(!$taln or !$qaln);
		next if(length($taln) != length($qaln));
		my ($mic_sam, $ctl_sam);
		my @ies = PARTIES::Utils->get_InDel_from_aln($seq,$taln, $qname,$qaln, 0, $qstart_tmp, length($contigseq), '+',
   						{ 
						JUNCTION_FLANK_SEQ_LENGTH => $self->{JUNCTION_FLANK_SEQ_LENGTH},
						ERROR_FILE_HANDLER => $self->{ERR},
						NOT_BOUNDED_BY_TA =>   $self->{NOT_BOUNDED_BY_TA},
						CHECK_BREAK_POINTS => 0,
						MIN_INDEL_LENGTH => $self->{MIN_SIZE},
						MIRAA_BREAK_POINTS => ($self->{BREAK_POINTS}) ? $self->{BREAK_POINTS}->{$tname} : '',
						SAMPLE_BAM => ($mic_sam) ? $mic_sam : 0,
						CONTROL_BAM => ($ctl_sam) ? $ctl_sam : 0,
						MIN_READS_SUPPORTING_MAC_INDEL_JUNCTIONS => 1,
						TNAME=>$tname,
						 } );

		if(defined($ies[0]->{POS}) && $ies[0]->{POS} > $self->{JUNCTION_FLANK_SEQ_LENGTH} && $ies[0]->{POS} < (length($contigseq) - $self->{JUNCTION_FLANK_SEQ_LENGTH})){
			my $shift = $tstart_tmp - $ies[0]->{POS};
			$segments->{$seg_id}->{bounded_by_ta}= $ies[0]->{BOUNDED_BY_TA};
			$segments->{$seg_id}->{start} = $segments->{$seg_id}->{start} - $shift;
			$segments->{$seg_id}->{end} = $segments->{$seg_id}->{end} - $shift;
			$segments->{$seg_id}->{sequence} = $ies[0]->{SEQ};
			push @{$segments->{$seg_id}->{IES_found}}, @ies;
			$segments->{$seg_id}->{junction_seq} = PARTIES::Utils->get_junction_seq($seq, $ies[0]->{POS}, $self->{JUNCTION_FLANK_SEQ_LENGTH} );
			$segments->{$seg_id}->{id} = PARTIES::Utils->generate_IES_id($self->{PREFIX},$seq_id,$segments->{$seg_id}->{start}, $segments->{$seg_id}->{end});
			$segments->{$seg_id}->{support_ref}=$self->_get_reference_support($sam, $seq_id, $segments->{$seg_id}->{start});		
#			$self->_compare_segment_to_ies($seq_id, $segments->{$seg_id});
			$segments->{$seg_id}->{tobereported}=1;
		}
		else{
			$segments->{$seg_id}->{tobereported}=0;
		}
		
		$a++;
	}

}


sub get_min_dist{ # return the shortest distance of a TA indel to the given IES
	my ($s1, $e1, $s2, $e2)=@_;
	my @res = (abs($s2-$s1), abs($e2-$e1));
	return( min(@res) );
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
	if($read->get_tag_values('XM') > 0 ){return (-1,-1);} # N-containing read
	if(!$read->paired || !defined($read->mate_seq_id)){ return (-1,-1);} # Unpaired read
	if($read->seq_id ne $read->mate_seq_id){ return (-1,-1);} # Read mapped on different seq_id
	if($read->strand == $read->mstrand){return (-1,-1);} # Strand inconsistency between pair
	
	if(	 $cigar =~ /^(\d+)M(\d+)S$/){			return('MS' ,$2,$1, $read->query->dna);}
	elsif ($cigar =~ /^(\d+)S(\d+)M$/ ){		return('SM' ,$1,$2, $read->query->dna);}
	elsif ($cigar =~ /^(\d+)M(\d+)D(\d+)M$/ ){	return('MDM',$2,$1, $read->query->dna);}
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
