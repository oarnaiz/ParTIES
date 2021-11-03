package PARTIES::MICA;
use strict;
use base 'PARTIES::Root';
use PARTIES::Config;
use File::Basename;
use List::Util qw(max min shuffle);

=head1 NAME

 ParTIES MICA module - Method of Identification by Comparison of Assemblies

=head1 AUTHORS

2015, I2BC - CNRS , GNU GPL3

=cut

# ALL PARAMETERS FOR THIS MODULE
my %PARAMETERS = (
			BAM => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Mapping on the reference genome"
				},
			CONTROL_BAM => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Mapping of somatic DNA-sequencing on the somatic genome"
				},
			GERMLINE_GENOME => {
				MANDATORY=>1, DEFAULT=>[], TYPE=>'MULTIPLE', RANK => 1,
				DESCRIPTION=>"Assembly file in which IES will be searched (may be used multiple time)"
				},
			GERMLINE_BLAT => {
				MANDATORY=>0, DEFAULT=>[], TYPE=>'MULTIPLE', RANK => 2,
				DESCRIPTION=>"Blat of a given germline genome on the reference genome (may be used multiple time; it will be computed if not given)"
				},
			MIRAA => {
				MANDATORY=>0, DEFAULT=>'', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"MIRAA output used to filter reads"
				},			
			INSERT_SIZE => {
				MANDATORY=>0, DEFAULT=>'500', TYPE=>'VALUE', RANK => 1,
				DESCRIPTION=>"Estimation of the sequencing insert size"
				},
			SEQUENCES_PER_BLAT_BATCH => {
				MANDATORY=>0, DEFAULT=>50, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Number of sequences per blat batch"
				},				
			PREFIX => {
				MANDATORY=>0, DEFAULT=>'MICA', TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Prefix of the ID of detected deletions"
				},			
			JUNCTION_FLANK_SEQ_LENGTH => {
				MANDATORY=>0, DEFAULT=>15, TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"Length of the flanking sequence to report, arround the IES"
				},
			NOT_BOUNDED_BY_TA => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 2,
				DESCRIPTION=>"Should insertion not bounded by TA dinucleotide be reported"
				},
			REPEAT_MASKER_PARAMETERS => {
				MANDATORY=>1, DEFAULT=>'-nolow  -x -species "paramecium tetraurelia"', TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"RepeatMasker parameters"
				},
			SKIP_REPEAT_MASKER => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Skip RepeatMasker step"
				},
		);


my %FILE_EXTENSIONS = ( 
			IES => {DESC=> 'IES GFF3 file', EXT => 'gff3' },
			IES_PLUS => {DESC=> 'IES GFF3 file', EXT => 'ies.gff3' }
			 );

my $SHIFT = 1000;
my $LIMIT = 25000;
#my $MAC_FLANK_SEQ_LENGTH = 15;

=head2 new

 Title   : new
 Usage   : my $factory = PARTIE::MICA->new({....})
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
  
     my $miraa = $self->{OUT_DIR}."/MIRAA/MIRAA.gff3";
     $self->{MIRAA} = $miraa if(!$self->{MIRAA} and -e $miraa);

     if(!scalar @{$self->{GERMLINE_GENOME}}) {
        my $assembly_dir = $self->{OUT_DIR}."/Assembly/";
        if(-e $assembly_dir ) {
           foreach my $assembly (`ls $assembly_dir/VELVET_*/*fa`) {
              chomp $assembly;
              push @{$self->{GERMLINE_GENOME}},$assembly;
           }
        }
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
  
  # READ FILES
    
  # load MIRAA file
  if( $self->{MIRAA}) {
     $self->stderr("Read ".basename($self->{MIRAA})." ... " );
     $self->{LOADED_MIRAA} = PARTIES::Utils->read_gff_file($self->{MIRAA});
     $self->stderr("Done\n" );
  }
  
  
  #die "Need same number of blat alignment than assemblies" if($self->{BLAT} and scalar @{$self->{GERMLINE_GENOME}} != scalar @{$self->{BLAT}});
  
  # Mask assemblies
  $self->stderr("Mask assemblies ... " );
  my ($RepeatMasker) = (PARTIES::Config->get_program_path('RepeatMasker'));
  my $repeat_masker_parameters = $self->{REPEAT_MASKER_PARAMETERS};
  for(my $aln_number=0;$aln_number<@{$self->{GERMLINE_GENOME}};$aln_number++) {
     my $assembly = $self->{GERMLINE_GENOME}->[$aln_number];  
     my $base_name_assembly = basename($assembly);
     my $mask_dir = $self->{PATH}."/tmp/mask";
     mkdir $mask_dir;
     next if(-e "$assembly.masked");
     $self->stderr( $base_name_assembly);
     #$self->stderr("Masking ".$base_name_assembly." ... " );
     if($self->{SKIP_REPEAT_MASKER} eq 'TRUE') {
        system("ln -s $base_name_assembly $assembly.masked");
     } else {
        system("$RepeatMasker -pa $self->{THREADS} -dir $mask_dir/ $repeat_masker_parameters $assembly > /dev/null");
        # if no repetitive sequences
        if(!-e "$mask_dir/$base_name_assembly.masked") {
           my $rm_out = `head -1 $mask_dir/$base_name_assembly.out`;
	   system("ln -s $base_name_assembly $assembly.masked") if($rm_out=~/^There were no repetitive sequences detected in /);
        } else {
           die "No file $mask_dir/$base_name_assembly.masked" if(!-e "$mask_dir/$base_name_assembly.masked");
           system("mv $mask_dir/$base_name_assembly.masked $assembly.masked");  
	}
     }
     
     #$self->stderr("Done\n" );
  }
  $self->stderr("Done\n" );
   
  # multi thread blat calculation
  $self->stderr("Align assemblies on the reference\n" );
  #my $pm = new Parallel::ForkManager($self->{THREADS});
  my @blat_files = @{$self->{GERMLINE_BLAT}} if($self->{GERMLINE_BLAT});
  for(my $aln_number=0;$aln_number<@{$self->{GERMLINE_GENOME}};$aln_number++) {
     my $assembly = $self->{GERMLINE_GENOME}->[$aln_number].".masked";
     my $blat_file = $blat_files[$aln_number];
     next if(defined $blat_file and $blat_file ne '' and -e $blat_file);
     #my $pid = $pm->start and next; 
     $self->_calculate_blat_file($assembly,$blat_file);
     #$pm->finish(0);
  }
  #$pm->wait_all_children; 
  
  
  @blat_files = @{$self->{GERMLINE_BLAT}} if($self->{GERMLINE_BLAT});
  for(my $aln_number=0;$aln_number<@{$self->{GERMLINE_GENOME}};$aln_number++) {
     my $assembly = $self->{GERMLINE_GENOME}->[$aln_number];
     my $base_name_assembly = basename($assembly);
     my $blat_file = ($blat_files[$aln_number]) ? $blat_files[$aln_number] : $self->{PATH}."/".$base_name_assembly.".masked_VS_".basename($self->{GENOME}).".psl";

     # load assembly fasta sequences
     $self->stderr("Read ".$base_name_assembly." ... " );
     my $seqio = new Bio::SeqIO(-file=> $assembly,-format=>'Fasta');
     while(my $seq = $seqio->next_seq) {
        my $seq_id = $seq->id;
        $self->{CONTIGS}->{$aln_number}->{$seq_id}=$seq->seq;
     } 
     $seqio->close; 
     $self->stderr("Done\n" );

     
  
     $self->stderr("Read blat alignment : ".basename($blat_file)." ... " );
     open(BLAT,$blat_file) or die $blat_file;
     my %blat;
     while(<BLAT>) {
        chomp;
        next if($_!~/^\d/);
        my ($match,$mismatch,$repmatch,$n,
	  $qgap_count,$qgap_bases,$hgap_count,$hgap_bases,$strand,
	  $qname,$qsize,$qstart,$qend,
	  $tname,$tsize,$tstart,$tend,
	  $block_count,$blockSizes,$qStarts,$tStarts) = split /\t/,$_;
        die $_ if(!$self->{CONTIGS}->{$aln_number}->{$qname});
        my $score = $match - $mismatch;
        #next if($strand eq '-');
        # FIND THE BEST MATCHES
        if(!$blat{$qname} or $score > $blat{$qname}->{score}) {
           $blat{$qname}->{score}=$score;
           my @lines = ($_);
           $blat{$qname}->{line}=\@lines;
        } elsif($score == $blat{$qname}->{score}) {
           push @{$blat{$qname}->{line}},$_;
        }
      }
      close BLAT;
      foreach my $qname (keys %blat) {
         foreach my $line (@{$blat{$qname}->{line}}) {
           my ($match,$mismatch,$repmatch,$n,
           $qgap_count,$qgap_bases,$hgap_count,$hgap_bases,$strand,
           $qname,$qsize,$qstart,$qend,
           $tname,$tsize,$tstart,$tend,
           $block_count,$blockSizes,$qStarts,$tStarts) = split /\t/, $line;
	      push @{$self->{ALIGNMENTS}->{$aln_number}->{$tname}},$line;
         }
      }
      $self->stderr("Done\n" );
      
   }    
   
}


# do a blat alignment
sub _calculate_blat_file {
  my ($self,$assembly,$blat_file) = @_;  
  
  my $base_name_assembly = basename($assembly);
  # BLAT ALIGNMENT  
  if(!$blat_file or !-e $blat_file) {


     $self->stderr("Alignment of ".$base_name_assembly." ON ".basename($self->{GENOME})." ... "  );
     $blat_file=$self->{PATH}."/".$base_name_assembly."_VS_".basename($self->{GENOME}).".psl";
     if(-e $blat_file) {
     	$self->stderr(" already calculated ... ");
     } else {
     	my $insert_size = $self->{INSERT_SIZE};
	
        my $qio = new Bio::SeqIO(-file =>$assembly ,-format => 'Fasta');
	my @qseqs;
        while (my $q = $qio->next_seq) {
           push @qseqs , $q;
        }
	$qio->close;
	
	mkdir $self->{PATH}."/tmp/".$base_name_assembly;
	# write fasta
	my $number_of_sequence_by_batch = $self->{SEQUENCES_PER_BLAT_BATCH};
	my @qseq_ids = shuffle(@qseqs);
	
	my $out;
        my $n=0;
        my @base_files;
	foreach my $q (@qseq_ids) {
	
           if($n == 0 or $n >= $number_of_sequence_by_batch) {
               $out->close if($out);
	       my $base_filename = File::Temp->new(DIR => $self->{PATH}."/tmp/".$base_name_assembly."/")->filename;
	       while(-e "$base_filename.fa") {
	          $base_filename = File::Temp->new(DIR => $self->{PATH}."/tmp/".$base_name_assembly."/")->filename;
	       }
               $out = new Bio::SeqIO(-file =>">$base_filename.fa" ,-format => 'Fasta');
               push @base_files,$base_filename;
               $n=0;
           }
           $out->write_seq($q);
           $n++;
        }
        $out->close if($out);
	
	# execute blat
	my $pm = new Parallel::ForkManager($self->{THREADS});
        foreach my $file (@base_files) {
           my $pid = $pm->start and next; 
	   system(PARTIES::Config->get_program_path('blat')." ".$self->{GENOME}." $file.fa -noHead -t=dna -q=dna  -minScore=$insert_size $file.blat > /dev/null");  
	   unlink  "$file.fa";
           # Terminates the child process
           $pm->finish(0);
        }
        $pm->wait_all_children; 

	# compile all results
	system("rm -f $blat_file && touch $blat_file");
	foreach my $base_file (@base_files) {
	   system("cat $base_file.blat >> $blat_file");
	}
        system("rm -rf ".$self->{PATH}."/tmp/".$base_name_assembly);

	
        #system(PARTIES::Config->get_program_path('blat')." ".$self->{GENOME}." ".$assembly." -t=dna -q=dna -noHead -minScore=$insert_size $blat_file");
	
	
     }
     $self->stderr("Done\n" );
  
  }
  return $blat_file;
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
   
   my $mic_sam = new Bio::DB::Sam(-bam  =>$self->{BAM},-fasta=>$self->{GENOME}) if($self->{BAM});
   my $ctl_sam = new Bio::DB::Sam(-bam  =>$self->{CONTROL_BAM},-fasta=>$self->{GENOME}) if($self->{CONTROL_BAM});   

   # THE HEADER MUST START BY "#"
   
   
   my %IES;
   # for each alignment search for insertion
   foreach my $aln_number (sort {$a<=>$b} keys %{$self->{ALIGNMENTS}}) {

      my %ies_in_seq_id;
      
      foreach my $blat_line (@{$self->{ALIGNMENTS}->{$aln_number}->{$seq_id}}) {
      
         my ($match,$mismatch,$repmatch,$n,
           $qgap_count,$qgap_bases,$hgap_count,$hgap_bases,$strand,
           $qname,$qsize,$qstart,$qend,
           $tname,$tsize,$tstart,$tend,
           $block_count,$blockSizes,$qStarts,$tStarts) = split /\t/, $blat_line;
	   
          my @blockSizes = split /,/,$blockSizes;
          my @qStarts = split /,/,$qStarts;
          my @tStarts = split /,/,$tStarts;
    
          $qstart = $qStarts[0]; # BUG IN THE BLAT OUTPUT
          $qend = $qStarts[$#qStarts]+$blockSizes[$#blockSizes];      
          
	  # Global alignment
          # if the segment is small, we can do a muscle alignment with all the segment
          if(($qend-$qstart) < $LIMIT and $hgap_bases < 10000) {
	     # find ies in the segment
             my @ies_in_segments = $self->_find_ies_in_segment($aln_number,$seq,$tstart,$tend,$qname,$qstart,$qend,$strand,'GLOBAL',$mic_sam,$ctl_sam,\*LOG,\*ERR) ;
             foreach my $ies (@ies_in_segments) {
                $ies_in_seq_id{$ies->{POS}}=$ies;
             }
	  # analyse block by block   
          } else {
             
	     my $npart = @qStarts;
	     
	     # create blocks
	     my @blocks;
             for(my $p=0; $p<$npart; $p++) {
                my($qstart_part,$tstart_part,$len_part)= ($qStarts[$p], $tStarts[$p], $blockSizes[$p]);
                my $qend_part= $qstart_part+$len_part;
                my $tend_part= $tstart_part+$len_part;
                $qstart_part++; $tstart_part++; # move to 1-origin

	        push @blocks , { QNAME => $qname, QSTART => $qstart_part, QEND =>$qend_part, TNAME => $tname, TSTART => $tstart_part, TEND => $tend_part};
		
	     }
	     
	     my @ordered_bloks = sort {$a->{TSTART}<=>$b->{TSTART}} @blocks; 
	     
	     for(my $b=0; $b<@ordered_bloks-1; $b++) {
	        my ($tname,$qname,$tpos) = ($ordered_bloks[$b]->{TNAME},$ordered_bloks[$b]->{QNAME},$ordered_bloks[$b]->{TEND});
		
		my ($left_shift,$right_shift) = ($SHIFT,$SHIFT);
		if($tpos < $SHIFT or ($tpos - $ordered_bloks[$b]->{TSTART}) < $SHIFT) {
		   $left_shift = min( $tpos,($tpos - $ordered_bloks[$b]->{TSTART}));
		}
		if(($tpos+$SHIFT) > $seq->length or ($ordered_bloks[$b+1] and ($ordered_bloks[$b+1]->{TEND} - $tpos) < $SHIFT) ) {
		   $right_shift = min(($seq->length-$tpos), ($ordered_bloks[$b+1]->{TEND} - $tpos));
		}
		
		my ($tstart_tmp,$tend_tmp) = ($tpos-$left_shift ,$tpos+$right_shift);
		my ($qstart_tmp,$qend_tmp) = ($ordered_bloks[$b]->{QEND}-$left_shift ,$ordered_bloks[$b+1]->{QSTART}+$right_shift);
		
		
		# find ies in the segment
		my @ies_in_segments = $self->_find_ies_in_segment($aln_number,$seq,$tstart_tmp,$tend_tmp,$qname,$qstart_tmp,$qend_tmp,$strand,'LOCAL',$mic_sam,$ctl_sam,\*LOG,\*ERR);
	        foreach my $ies (@ies_in_segments) {
                   $ies_in_seq_id{$ies->{POS}}=$ies;
                }
	     }
	     
          }      
      }
      
      # Merge IESs
      my $N=0;
      foreach my $id (sort {$a <=> $b} keys %ies_in_seq_id) {
         my $ies = $ies_in_seq_id{$id};
         next if($ies->{POS} < 1 or $ies->{POS} > $seq->length);
	 
         my ($pos,$ies_seq) = ($ies->{POS},$ies->{SEQ});
	 
         if($IES{($pos-1)}) {
            $pos--;
         } elsif($IES{($pos+1)}) {
            $pos++;
         }
      
         $IES{$pos}->{$ies_seq}->{OBJECT} = $ies;
         $IES{$pos}->{$ies_seq}->{EVIDENCES}++;
         $N++;
     
     
      }
      print LOG "NUMBER OF IES in $seq_id for ".basename($self->{GERMLINE_GENOME}->[$aln_number])." : $N\n";
      
      
   }
   my $WRITE_IES_GFF3 = (scalar @{$self->{GERMLINE_GENOME}} == 1) ? 1 : 0;
   
   
   my $gff_file = $self->{PATH}."/tmp/$seq_id.IES";
   open(GFF,">$gff_file") or die "Can not open $gff_file";
   print GFF "# ",join(" ",$self->_command_line),"\n";
   
   if($WRITE_IES_GFF3) {
      my $gff_file = $self->{PATH}."/tmp/$seq_id.IES_PLUS";
      open(GFFIES,">$gff_file") or die "Can not open $gff_file";
      print GFFIES "# ",join(" ",$self->_command_line),"\n";
   
   }

   my $N=0;
   my $mac_flank_seq_length = $self->{JUNCTION_FLANK_SEQ_LENGTH};
   foreach my $pos (sort {$a<=>$b} keys %IES) {
      next if($pos < ($mac_flank_seq_length+1) or $pos > ($seq->length-$mac_flank_seq_length-1));
      my $id = PARTIES::Utils->generate_IES_id($self->{PREFIX},$seq_id,$pos);
      
      # get best IES
      my ($best_ies,$evidences);
      foreach my $key (keys %{$IES{$pos}}) {
         if(!$best_ies or $IES{$pos}->{$key}->{EVIDENCES} > $evidences) {
	    $best_ies = $IES{$pos}->{$key}->{OBJECT} ;
	    $evidences = $IES{$pos}->{$key}->{EVIDENCES};
	 }
      }
      
      my $ies_seq = $best_ies->{SEQ};
      my $junction_seq = PARTIES::Utils->get_junction_seq($seq,$best_ies->{POS},$mac_flank_seq_length);
      my $bounded_by_TA = $best_ies->{BOUNDED_BY_TA};
      
      if(!$ies_seq or !$junction_seq or ($bounded_by_TA eq 'TRUE' and $junction_seq!~/TA/)) {
         print ERR "ERROR id=$id ies_seq=$ies_seq junction_seq=$junction_seq bounded_by_TA=$bounded_by_TA" ;
	 next;
      }
      
      my @attributes = ("ID=$id","Name=$id","junction_seq=$junction_seq","sequence=$ies_seq","bounded_by_ta=$bounded_by_TA");
      

      foreach my $alt_ies_seq (keys %{$IES{$pos}}) {
         next if($alt_ies_seq eq $ies_seq);
	 push @attributes,"alternative_seq=$alt_ies_seq";
      }
      

      print GFF join("\t",($seq_id,'MICA','internal_eliminated_sequence',$best_ies->{POS},($best_ies->{POS}+1),'.','.','.',join(";",@attributes))),"\n";
      
      if($WRITE_IES_GFF3) {
         print GFFIES join("\t",($best_ies->{QNAME},'MICA','internal_eliminated_sequence',$best_ies->{QPOS},($best_ies->{QPOS}+length($ies_seq)-1),'.','.','.',join(";",@attributes))),"\n";
      
      }
      
   }
   

   
   close GFF;
   close GFFIES if($WRITE_IES_GFF3);
   close LOG;
   close ERR;
   
   $self->stderr("End calculation $seq_id\n" );

   
}

# The most important method 
# find inserted sequences in the region, based on th alignment
sub _find_ies_in_segment {

   my ($self,$aln_number,$seq,$tstart_tmp,$tend_tmp,$qname,$qstart_tmp,$qend_tmp,$strand,$aln_type,$mic_sam,$ctl_sam,$fh_log,$fh_err) = @_;

   my $tmp_fname = File::Temp->new(DIR => $self->{PATH}."/tmp/")->filename;
   #$tmp_fname=$self->{PATH}."/tmp/N3ufHJWuE";
   
   while(-e "$tmp_fname.fa") {
      $tmp_fname = File::Temp->new(DIR => $self->{PATH}."/tmp/")->filename;
   }
   
   my $tname = $seq->id;
   

   my $scafseq = $seq->seq;
   my $contigseq = $self->{CONTIGS}->{$aln_number}->{$qname};
   if ($strand eq '-') { $contigseq =~ tr/ACGTacgt/TGCAtgca/; $contigseq = reverse $contigseq; }

   my $tseq = substr $scafseq,$tstart_tmp, ($tend_tmp-$tstart_tmp);
   my $ftseq = $tseq;
   $ftseq=~s/(.{60})/$1\n/g;
   my $qseq = substr $contigseq,$qstart_tmp, ($qend_tmp-$qstart_tmp);
   #if ($strand eq '-') { $qseq =~ tr/ACGT/TGCA/; $qseq = reverse $qseq; }
   $qseq=~s/(.{60})/$1\n/g;
   open(TMP,">$tmp_fname.fa") or die "Can not open $tmp_fname.fa $tname VS $qname \n $!\n";
   print TMP ">$tname:$tstart_tmp-$tend_tmp\n$ftseq\n";
   print TMP ">$qname:$qstart_tmp-$qend_tmp\n$qseq\n";
   close TMP;
   die "No seq\n$qname and $tname" if(!$ftseq or !$qseq);
   system(PARTIES::Config->get_program_path('muscle')." -in $tmp_fname.fa -out $tmp_fname.aln -quiet -maxiters 1 -diags");
   unlink "$tmp_fname.fa";
   
   my @ies;
   if(-e "$tmp_fname.aln") {    
   
   	#print STDERR "muscle -in $tmp_fname.fa -out $tmp_fname.aln -quiet -maxiters 1 -diags -clwstrict\n";
   	#system("muscle -in $tmp_fname.fa -out $tmp_fname.aln -quiet -maxiters 1 -diags -clwstrict");
   	open(ALN,"$tmp_fname.aln") or die "No file $tmp_fname.aln";
   	my %aln;
   	my $nb_aln=0;
   	while(<ALN>) {
   		chomp;
   		next if(!$_);
   		if($_=~/^>/) {
   		   $nb_aln++;
   		} else {
   		   $aln{$nb_aln}.=$_;
   		}
   	}
   	close ALN;
   	my $taln = $aln{1};
   	my $qaln = $aln{2};     
   	die "No seq align seq1=$taln seq2=$qaln in $tmp_fname.aln " if(!$taln or !$qaln);


   

   	# common procedure to find ies/indel in an alignment
   	@ies = PARTIES::Utils->get_InDel_from_aln($seq,$taln, $qname,$qaln, $tstart_tmp, $qstart_tmp, length($contigseq), $strand,
   						{ 
						JUNCTION_FLANK_SEQ_LENGTH => $self->{JUNCTION_FLANK_SEQ_LENGTH},
						ERROR_FILE_HANDLER => $fh_err,
						NOT_BOUNDED_BY_TA =>   $self->{NOT_BOUNDED_BY_TA},
						CHECK_BREAK_POINTS => 1,
						MIN_INDEL_LENGTH => 20,
						MIRAA_BREAK_POINTS => ($self->{LOADED_MIRAA}) ? $self->{LOADED_MIRAA}->{$tname} : '',
						SAMPLE_BAM => ($mic_sam) ? $mic_sam : 0,
						CONTROL_BAM => ($ctl_sam) ? $ctl_sam : 0,
						MIN_READS_SUPPORTING_MAC_INDEL_JUNCTIONS => 1,
						 } );
   } else {
   	$self->stderr("Warning no file $tmp_fname.aln\n" );
   }
   unlink "$tmp_fname.aln";
   return @ies;

}   






1;
