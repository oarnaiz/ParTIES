package PARTIES::MIRET;
use strict;
use base 'PARTIES::Root';

use PARTIES::Config;
use Data::Dumper;
use File::Basename;
use FindBin qw($Bin);

=head1 NAME

 ParTIES MIRET module - Method of Ies RETention

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Mapping file of reads on the genome"
				},
			GERMLINE_BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Description parameter"
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
#			SNPFILE => {
#				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 3,
#				DESCRIPTION=>"File of SNPs on the Macronuclear genome with IES"
#				},
			MAX_MISMATCH => {
				MANDATORY=>0, DEFAULT=> 1, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Number of mismatch allowed in a read."
				},
			SCORE_METHOD => {
				MANDATORY=>1, DEFAULT=>"Boundaries", TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"This indicate how the retention score should be calculated, either for each Boundaries or for IES [Boundaries;IES]"
				},
			CONTROL => {
				MANDATORY=>0, DEFAULT=>"", TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"IES retention in a control sample. It will be used, if provided to test the upper retention in the current sample compare to the control sample (GFF3 MIRET output)"
				},
			TOPHAT => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"RNAseq mapped with tophat"
				},			
			TAB => {
				MANDATORY=>0, DEFAULT=>'TRUE', TYPE=>'BOOLEAN', RANK => 2,
				DESCRIPTION=>"Write results in tabulated format"
				},

		);



my %FILE_EXTENSIONS = ( 
			gff3 => {DESC=> 'GFF3 file', EXT => 'gff3' },		
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
    
     push @{$self->{GERMLINE_GENOME}},$self->{OUT_DIR}."/Insert/Insert.fa" if(!@{$self->{GERMLINE_GENOME}} and -e $self->{OUT_DIR}."/Insert/Insert.fa");

     my $germline_bam = $self->{OUT_DIR}."/Map/".basename($self->{OUT_DIR}).".".basename(${$self->{GERMLINE_GENOME}}[0]).".BOWTIE.sorted.bam";
     $self->{GERMLINE_BAM} = $germline_bam if(!$self->{GERMLINE_BAM} and -e $germline_bam);
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

  # Loading IES file with germline genome coordinates
  $self->stderr("Read ".basename($self->{GERMLINE_IES})." ...\n" );
  $self->{LOADED_GERMLINE_IES} = PARTIES::Utils->read_gff_file_by_id($self->{GERMLINE_IES});
  $self->stderr("Done\n" );
	
  # bin directory
  $self->{BIN}=$Bin;
  
 #Loading the SNP tab file
#	if($self->{'SNPFILE'}){
#		my $snpfile=\$self->{'SNPFILE'};
#		open(SNP, "<".$$snpfile) or die $!;
#		while(<SNP>){
#			chomp;
#			my ($seq_id, $pos, $refs, $covers, $pileups, $quals, $majorbases, $freqs)=split("\t", $_);
#			$self->{SNP}->{$seq_id}->{$pos}=lc($majorbases);
#		}
#		close(SNP);
#	}

  

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
   
   my $log_file = $self->{PATH}."/tmp/$seq_id.log";
   open(LOG,">$log_file") or die "Can not open $log_file";
   print LOG "LOG MESSAGE\n";

   my $err_file = $self->{PATH}."/tmp/$seq_id.err";
   open(ERR,">$err_file") or die "Can not open $err_file";
   if(defined($self->{LOADED_IES}->{$seq_id})){

   	my $sam = new Bio::DB::Sam( -bam => $self->{BAM}, -fasta => $self->{GENOME}, -expand_flags => 'true');
		my $sam_ies = new Bio::DB::Sam( -bam => $self->{GERMLINE_BAM}, -fasta => $self->{GERMLINE_GENOME}[0], , -expand_flags=> 'true' );
		my %read_cat;
		my $shift=0;
		my $min_dist=4;
		foreach my $ies (sort {$a->{start}<=>$b->{start}} @{$self->{LOADED_IES}->{$seq_id}}) {
			my $ies_id=\$ies->{attributes}->{ID}[0];
			my $ies_seq=\$ies->{attributes}->{sequence}[0];
			my $junction_seq=\$ies->{attributes}->{junction_seq}[0];
			my $ies_size=length($$ies_seq);
			my $position = $ies->{start};
			my ($twd_left, $twd_right)=_ambigous_ies_mac_sequence(\$ies_seq, \$junction_seq);

			my ($ies_start, $ies_end) = ($self->{LOADED_GERMLINE_IES}->{$$ies_id}[0]->{start}, $self->{LOADED_GERMLINE_IES}->{$$ies_id}[0]->{end}+1);
			my $seq_id_ies=$self->{LOADED_GERMLINE_IES}->{$$ies_id}[0]->{seq_id};
			$twd_left= $min_dist if ($twd_left<$min_dist);
			$twd_right= $min_dist if ($twd_right<$min_dist);

			$read_cat{$position}={};
		# Getting mac counts
			my ($mac, $mapping_issued, $total, $mac_read_names);
			($mac, $mapping_issued, $total, $mac_read_names,
				$read_cat{$position} )=_get_counts(\$self,
																$position, 
																\$sam->segment($seq_id, $position-$twd_left, $position+$twd_right), 
																$twd_left, 
																$twd_right, 
																{"MAC"=>1}, 
																{"RIGHT"=>0},
																\%{$read_cat{$position}});

	# Getting Left Boundaries IES+ counts
			my ($L_ret, $L_mapping_issues, $L_total, $L_read_names);
			($L_ret, $L_mapping_issues, $L_total, $L_read_names,
				$read_cat{$position})=_get_counts(	\$self,
																$ies_start, 
																\$sam_ies->segment($seq_id_ies, $ies_start-$twd_left, $ies_start+$twd_right),
																$twd_left, 
																$twd_right, 
																$mac_read_names,
																{"RIGHT"=>1},
																\%{$read_cat{$position}});

	# Getting Right Boundaries IES+ counts
			my ($R_ret, $R_mapping_issues, $R_total, $R_read_names);
			($R_ret, $R_mapping_issues, $R_total, $R_read_names,
				$read_cat{$position}) =_get_counts(	\$self,
																$ies_end, 
																\$sam_ies->segment($seq_id_ies, $ies_end-$twd_left, $ies_end+$twd_right),
																$twd_left, 
																$twd_right, 
																$mac_read_names,
																$L_read_names,
																\%{$read_cat{$position}});

			my $ies_ret=$L_ret+$R_ret;
			my $L_bound_ret=$L_ret;
			my $R_bound_ret=$R_mapping_issues;

			if($self->{SCORE_METHOD} eq "Boundaries"){
				my ($L_score,$R_score)=(0,0);
				$L_score=sprintf('%.4f',($L_bound_ret/($L_bound_ret+$mac))) if($L_bound_ret+$mac != 0);
				$R_score=sprintf('%.4f',($R_bound_ret/($R_bound_ret+$mac))) if($R_bound_ret+$mac != 0);
				$ies->{attributes}->{support_mac}[0]=$mac;
				$ies->{attributes}->{support_left}[0]=$L_bound_ret;
				$ies->{attributes}->{support_right}[0]=$R_bound_ret;
				$ies->{attributes}->{retention_score_left}[0]=$L_score;
				$ies->{attributes}->{retention_score_right}[0]=$R_score;
				$ies->{attributes}->{total_counts_left}[0]=$L_total;
				$ies->{attributes}->{total_counts_right}[0]=$R_total;
			}
			elsif($self->{SCORE_METHOD} eq "IES"){
				my $Score=0; $Score=sprintf('%.4f',($ies_ret/($ies_ret+$mac))) if($ies_ret+$mac != 0);
				$ies->{attributes}->{support_mac}[0]=$mac;
				$ies->{attributes}->{support_ies}[0]=$ies_ret;
				$ies->{attributes}->{retention_score}[0]=$Score;
				$ies->{attributes}->{total_counts}[0]=$L_total+$R_total;
			}
			push @{$results{$self->get_mode()}->{$seq_id}},$ies;
		}
		
   }
   print LOG "END calculation for $seq_id\n";
   close LOG;
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
	
	
	my %significant;
	if($self->{CONTROL} ne '' and -e $self->{CONTROL}) {
	   
	   ## Load CONTROL
	   my $control_miret = PARTIES::Utils->read_gff_file_by_id($self->{CONTROL});
        
	   my $method = $self->{SCORE_METHOD};
	   my $tmp_file = $self->{PATH}."/tmp/".uc($self->get_mode).".tab";
	   
	   %significant = PARTIES::Utils->significant_retention_score($results,$control_miret,{ METHOD => $method,
	   											TMP_FILE => $tmp_file,
												BIN_DIR => $self->{BIN},
												}
												);
	   
	  	   
	
	}
 	
	
        ## Writes the final GFF file
	my $gff_file = $self->{PATH}."/".uc($self->get_mode).".gff3";
  	open(GFF,">$gff_file") or die "Can not open $gff_file";   
	print GFF "##gff-version 3\n# ",join(" ",$self->_command_line),"\n";

	foreach my $seq_id (keys %{$results}){
	
	   next if(!$results->{$seq_id});
	   foreach my $ies (@{$results->{$seq_id}}) {
	      my ($ies_id) = @{$ies->{attributes}->{ID}};
	      if($significant{$ies_id}) {
	         foreach my $key (keys %{$significant{$ies_id}}) {
		    # replace old value
		    $ies->{attributes}->{$key} = () if($ies->{attributes}->{$key});		    
		    push @{$ies->{attributes}->{$key}},$significant{$ies_id}->{$key};
		 }
	      }
	      print GFF PARTIES::Utils->gff_line($ies),"\n";
	   }
	}
	close GFF;
	
       if($self->{TAB} eq 'TRUE') {
          $self->stderr("Write tab file ...\n" );
          my $gff2tab=$self->{BIN}.'/utils/gff2tab.pl';
   
          foreach my $key (keys %FILE_EXTENSIONS) {
             my $fext = ($FILE_EXTENSIONS{$key}->{EXT}) ? $FILE_EXTENSIONS{$key}->{EXT} : $key;
             next if($fext!~/gff3$/);
             my $gff_file = $self->{PATH}."/".$self->get_mode.'.'.$fext;
             my $tab_file = $gff_file ;
             $tab_file=~s/\.gff3$/\.tab/;
             system("$gff2tab -gff $gff_file > $tab_file");
          }
       } 	
	
	
	
	
	$self->remove_temporary_files;

}




sub _ambigous_ies_mac_sequence {
# Gives the length for with the sequence between IES and MAC juction are different
	my ($ies_seq, $junction_seq)=@_;
	my ($toward_left, $toward_right) = (0,0);
	my ($Ldif, $Rdif)=(0,0);
	my $junction_side_length=(length($junction_seq)-2)/2;

	while( ($toward_left < $junction_side_length) and $Ldif==0){
		#substring junction_seq from just before TA to before TA + toward_left+1
		my $sub_junc = lc(substr($junction_seq, ($junction_side_length-1)-($toward_left), $toward_left+1));
		my $sub_ies = lc(substr($ies_seq, (length($ies_seq)-3) - ($toward_left), $toward_left+1));

		$Ldif=1 if($sub_junc ne $sub_ies);
		$toward_left++;
	}
	
	while( ($toward_right <= $junction_side_length) and $Rdif==0){
		#substring junction_seq from just after TA to after TA + toward_right+1
		my $sub_junc = lc(substr($junction_seq, ($junction_side_length+2), $toward_right+1));
		my $sub_ies = lc(substr($ies_seq, 2, $toward_right+1));

		$Rdif=1 if($sub_junc ne $sub_ies);
		$toward_right++;
	}
	return ($toward_left, $toward_right);
}




sub _get_counts {
# Require the position of the T from the TA
#			the segment object from a bio::db:sam
#			the min distances required to avoid ambiguities
#			the maximum % of mismatch in the read.
	my (	$self, $position, $segment,
			$min_left_dist, $min_right_dist,
			$mac_reads, $left_reads, $hash_reads) = @_ ;
	my ($mac, $nw_mapped, $total)=(0,0,0);
	my %read_names; 


	my %read_cat=%{$hash_reads};	# Read debugging
	my $step;
	$step="MAC" if(defined($mac_reads->{"MAC"}));
	$step="LEFT" if(!defined($mac_reads->{"MAC"}) && defined($left_reads->{"RIGHT"}) );
	$step="RIGHT" if(!defined($mac_reads->{"MAC"}) && !defined($left_reads->{"RIGHT"}) );
        
	my $TOPHAT_MAPPING = ($$self->{TOPHAT} eq 'TRUE') ? 1 : 0;
	foreach my $aln ($$segment->features() ){
		$total++;
		my ($rname, $rcigar, $rstart, $rend)=($aln->name, $aln->cigar_str, $aln->start, $aln->end);
		
		if(!defined($read_cat{$rname})){$read_cat{$rname}={ "MAC" => 'ns', "LEFT" => 'ns', "RIGHT" => 'ns'};}

# Filtering if the current read has been flag as Mac read
		if(defined($mac_reads->{$rname}) && $mac_reads->{$rname}==$aln->get_tag_values('FIRST_MATE')){ 
			$read_cat{$rname}->{$step}="flagged_as_mac"; $total--; next;}

# If RNAseq option is set : transform all \dM\dN\dM reads into \M reads.
		if($TOPHAT_MAPPING && $rcigar=~/(\d+)M(\d+)N(\d+)M/){
			my $fm=$1+$2+$3;
			$rcigar="$fm\M";
			if($step ne "MAC" &&
				$1+$rstart <= $position +$min_left_dist && 
				$1+$rstart+$2 >= $position+1 - $min_right_dist ){
					
					$read_cat{$rname}->{$step}="too_close_to_TA";
					next;
					
			}
		}

# Filtering on alignement profil (cigar), positions of the read compare to the TA
		if (_allowed_cigar($rcigar)==0){ 
			$read_cat{$rname}->{$step}="mapping_not_allowed";next;}
		if($rstart >= $position - $min_left_dist){ 
			$read_cat{$rname}->{$step}="too_close_to_TA";next;}
		if($rend <= $position+1 + $min_right_dist){ 
			$read_cat{$rname}->{$step}="too_close_to_TA"; next;}


# Filtering using the MD tag indicating where are the differences between read and reference
# Conting the number of differences, exluding those reported in the SNP file.
		my $mismatchs_found=_mismatch_count($self,\$aln);
		if($mismatchs_found > $$self->{'MAX_MISMATCH'}){ 
			$read_cat{$rname}->{$step}="too_many_mismatchs"; next;}
	

	# Recording read name
		if($step eq "MAC"){
			$read_names{$rname}=$aln->get_tag_values('FIRST_MATE');$mac++;
			$read_cat{$rname}->{$step}="MAC";
		}
		elsif($step eq "LEFT"){
			if(!defined($mac_reads->{$rname})){
				$read_names{$rname}=$aln->get_tag_values('FIRST_MATE');$mac++;
				$read_cat{$rname}->{$step}="MIC";
			}
		}
		else{ # Counting on the right Boundaries of the IES
			if(!defined($mac_reads->{$rname}) && !defined($left_reads->{$rname})){
				$read_names{$rname}=$aln->get_tag_values('FIRST_MATE');$mac++;
			}
			$read_cat{$rname}->{$step}="MIC";
			$nw_mapped++;
#			if($$self->{SCORE_METHOD} eq "Boundaries" && !defined($mac_reads->{$rname})){
#				$read_names{$rname}=$aln->get_tag_values('FIRST_MATE');$mac++;
#			}else{
#				$read_cat{$rname}->{$step}="MIC";
#				$nw_mapped++;
#				}
		}
	}
	return ($mac,$nw_mapped,$total,\%read_names, \%read_cat);
}

sub _allowed_cigar {
# Asses whether the cigar strand is compatible with further analysis
#	Allowed is full match, single nucleotide indel.
#	Returns 	0/1 for compatibility
	my ($cigar)=@_;
	if( $cigar =~ /^(\d+)M$/){return 1;}
	elsif( $cigar =~ /^(\d+)M$(\d+)I(\d+)M$/){return 1;}
	elsif( $cigar =~ /^(\d+)M$\d+D(\d+)M$/){return 1;}
	else{return 0; }
}


sub _mismatch_count {
	my ($self, $a)=@_;
	my $md=$$a->get_tag_values("MD");
	return(0) if $md=~/^\d+$/;
	
	#retrieving mismatch positions from MD tag and read start
	my $read_start=$$a->start;
	my $read_based_index=0;
	my $ref_based_index=$read_start;
	my $mismatch_count=0;
	while($md=~/([0-9]+)([A-Z]|\^[A-Z]+)/g)
	{ 
		my $new_base;
		my @letters=split(//,$2);
		
		$read_based_index+=$1;
		$ref_based_index+=$1;
		
		if(scalar(@letters)==1){
			$new_base=substr($$a->query->dna, $read_based_index, 1);
#			$mismatch_count++ if(lc($$self->{SNP}->{$$a->seq_id}) ne lc($new_base));
			$mismatch_count++;# if(lc($$self->{SNP}->{$$a->seq_id}) ne lc($new_base));
		# Incrementing both read-based and ref-based index
			$read_based_index++;
			$ref_based_index++;
		}
		else{
			foreach my $nt (@letters){
				next if($nt eq "^"); # the ^ char indicate a insertion in the reference compare to the read
				if($letters[0] eq "^"){ 
#					$mismatch_count++ if(lc($$self->{SNP}->{$$a->seq_id}) ne '*');
					$mismatch_count++;# if(lc($$self->{SNP}->{$$a->seq_id}) ne '*');
					#Incrementing only the ref-based index since there is no read sequence here
					$ref_based_index++;
				}
				else{
					$new_base=substr($$a->query->dna, $read_based_index, 1);
#					$mismatch_count++ if(lc($$self->{SNP}->{$$a->seq_id}) ne lc($new_base));
					$mismatch_count++;# if(lc($$self->{SNP}->{$$a->seq_id}) ne lc($new_base));
				# Incrementing both read-based and ref-based index
					$read_based_index++;
					$ref_based_index++;
				}
			}
		} 
	}
	return $mismatch_count;
}

1;
