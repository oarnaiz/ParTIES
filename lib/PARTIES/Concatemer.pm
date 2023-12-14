package PARTIES::Concatemer;
use strict;
use base 'PARTIES::Root';
use PARTIES::Config;
use File::Basename;
use List::MoreUtils qw(uniq); 

=head1 NAME

 ParTIES Concatemer module - Find trace of concatemere excision product

=head1 AUTHORS

2020, I2BC - CNRS , GNU GPL3

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
			IES => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 2,
				DESCRIPTION=>"IES file, usually MICA output"
				},
			BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 2.1,
				DESCRIPTION=>"Mapping on the reference genome"
				},
			GERMLINE_BAM => {
				MANDATORY=>1, DEFAULT=>'', TYPE=>'VALUE', RANK => 2.2,
				DESCRIPTION=>"Mapping file on the germline genome"
				},
			MAX_MISMATCH => {
				MANDATORY=>0, DEFAULT=>0, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Maximum mismatch in the alignment for a read to be used"
				},   
			MIN_MATCH_LENGTH => {
				MANDATORY=>0, DEFAULT=>6, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Minimum match length for a read to be used"
				},          
			MIN_MAPPING_QUALITY => {
				MANDATORY=>0, DEFAULT=>10, TYPE=>'VALUE', RANK => 3,
				DESCRIPTION=>"Minimum mapping quality"
				},          
			NO_REMOVE_PCR_DUPLICATES => {
				MANDATORY=>0, DEFAULT=>'FALSE', TYPE=>'BOOLEAN', RANK => 3,
				DESCRIPTION=>"Skip PCR duplicates removal procedure"
				},
		);



my %FILE_EXTENSIONS = ( 
			gff3 => {DESC=> 'GFF3 file', EXT => 'gff3' }
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
  foreach my $fastq ($self->{FASTQ1},$self->{FASTQ2}) {
     next if($fastq=~/\.fastq$/ or $fastq=~/.fastq.gz/);
     print STDERR "ERROR : fastq file (-fastq $fastq) should be a FASTQ file or does not exist\n" ;
     return 0;
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
    my $IES = PARTIES::Utils->read_gff_file($self->{IES});
    $self->stderr("Done\n" );

    my ($threads,$bowtie2_build,$bowtie2,$samtools) = ($self->{THREADS},PARTIES::Config->get_program_path('bowtie2-build'),PARTIES::Config->get_program_path('bowtie2'),PARTIES::Config->get_program_path('samtools'));
    my ($bwa) =PARTIES::Config->get_program_path('bwa');
    $self->stderr("Write IES fasta file ...\n" );
    my $ies_fasta_file = $self->{PATH}."/tmp/ies.fa";
    $self->{IES_FASTA_FILE} = $ies_fasta_file;
    
    my $min_ies_length = '';
    my $seqout = new Bio::SeqIO(-file => ">$ies_fasta_file", -format => 'fasta');
    foreach my $seq_id (sort keys %{$IES}) {
        foreach my $ies (sort {$a->{start}<=>$b->{start}} @{$IES->{$seq_id}}) {
            my ($id) = @{$ies->{attributes}->{ID}};
            my ($sequence) = ($ies->{attributes}->{sequence}) ? @{$ies->{attributes}->{sequence}} : "";
            $sequence.='TA';
            $seqout->write_seq(new Bio::Seq(-id=>$id,-seq=>$sequence));
            $min_ies_length=length($sequence) if($min_ies_length eq '' or $min_ies_length > length($sequence));
            $self->{IES_SEQ}->{$id}=$sequence;
        }
    }
    $seqout->close;
    
    # Find repeated boundaries
    my %boundaries;
    foreach my $ies_id (keys %{$self->{IES_SEQ}}) {
        my $ies_seq = $self->{IES_SEQ}->{$ies_id};
        for(my $size=$self->{MIN_MATCH_LENGTH}; $size <=$min_ies_length ; $size++) {
            foreach my $word (substr($ies_seq,0,$size),substr($ies_seq,length($ies_seq)-$size)) {
                $boundaries{$size}->{$word}++;
                $boundaries{$size}->{PARTIES::Utils->revcomp($word)}++;
            }
        }
    }
    foreach my $size (keys %boundaries) {
        foreach my $word (keys %{$boundaries{$size}}) {
            if($boundaries{$size}->{$word} == 1) {
                delete $boundaries{$size}->{$word};
            }
        }
    }
    $self->{REPEATED_BOUNDARIES}=\%boundaries;
    
    
    system("$bowtie2_build $ies_fasta_file $ies_fasta_file > /dev/null 2>&1 && $samtools faidx $ies_fasta_file");
    system("$bwa index $ies_fasta_file 2> /dev/null && $samtools faidx $ies_fasta_file");
    $self->stderr("Done\n" );



    $self->stderr("Read BAM files ...\n" );
    my %mapped_on_references;
    foreach my $bam_file ($self->{BAM},$self->{GERMLINE_BAM}) {
        $self->stderr("Read $bam_file ...\n" );
        if(!-e $bam_file) {
            print STDERR "\n# WARNING $bam_file does not exist ...";
            exit;
        } else {
            foreach my $seq_id (`samtools view -H $bam_file | grep '\@SQ' | awk '{ print \$2 }' | perl -p -e 's/^SN://'`) {
                foreach my $line (`samtools view -F 4 $bam_file $seq_id  `) {
                    #my ($qname,$rname,$cigar)= split /\t/,$line;
                    my ($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$tlen,$seq,$qual,@opt)= split /\t/,$line;
                    #die "$qname,$rname,$cigar";
                    next if($cigar!~/^\d+M$/);
                    $mapped_on_references{$qname}++ if($rname ne '*');
                }   
            }
        }
    }  
    $self->stderr("Done\n" );  
    
    my @fastq_files;
    # create filtered fastq files
    my $file_idx=1;
    my %reads;
    foreach my $fastq ($self->{FASTQ1},$self->{FASTQ2}) {
        my $file_basename=basename($fastq);
        $file_basename=~s/\.fastq$//;
        $file_basename=~s/\.fastq.gz$//;
        $self->stderr("Read ".basename($fastq)." ...\n" );
        my $outfile= $self->{PATH}."/$file_basename.not_well_mapped.fa";
        my $cmd = ($fastq=~/.fastq.gz/) ? "zcat $fastq |" : "cat $fastq |";
	
        open(FASTQ,$cmd) or die "Can not use $cmd";
        while(<FASTQ>) {
            chomp;
            if($_=~/^\@/) {
                my $name_line = $_;
                my ($read_name,$atts) = split / /,$name_line;
                $read_name=~s/\/\d$//;
                $read_name=~s/^\@//;

                my $seq = <FASTQ>;
                chomp $seq;
                my $n = <FASTQ>;
                my $qual = <FASTQ>;

                if(!$mapped_on_references{$read_name} and $seq!~/N/) {
                    #print OUT ">$read_name\n$seq\n";
                    $reads{$read_name}->{$file_basename}=$seq;
                    
                }
            }
        }
        close FASTQ;
        

        
        $self->stderr("Done\n" );
        push @fastq_files,{input_file=>$outfile,basename=>$file_basename,round=>1, file_idx=>$file_idx};
        $file_idx++;
    }
    
    
    # Remove PCR Duplicates
    if($self->{NO_REMOVE_PCR_DUPLICATES} eq 'FALSE') {
        $self->stderr("Remove PCR duplicates ...\n" );
        my %read_counts;
        my @rnames=keys %reads;
        foreach my $rname (@rnames) {
            my ($read_seq1, $read_seq2) = ($reads{$rname}->{$fastq_files[0]->{basename}},$reads{$rname}->{$fastq_files[1]->{basename}});
            
            if($read_seq1  and $read_seq2) {
                my $unique_seq = join('__',$read_seq1, $read_seq2);
                
                # remove this read
                if($read_counts{$unique_seq} ne '') {
                    $reads{$rname}=undef;
                    delete $reads{$rname};
                    
                }
                $read_counts{$unique_seq} = 1;
                
            }
            
        }   
        $self->stderr("Done\n" );   
    }
    
    foreach my $file_handler (@fastq_files) {
        my $file_basename = $file_handler->{basename};
        my $input_file = $file_handler->{input_file};
        
        $self->stderr("Write $input_file ...\n" );
        open(OUT,">$input_file") or die $input_file;
        foreach my $read_name (keys %reads) {
            next if(!$reads{$read_name} or !$reads{$read_name}->{$file_basename});
            print OUT ">$read_name\n",$reads{$read_name}->{$file_basename},"\n";
        }
        close OUT;
        $self->stderr("Done\n" ); 
    }
    
    
    $self->{INPUT_FILEs} = \@fastq_files;
    

    


}


sub _mapping {
    my ($self,$file_handler,$seq_fasta_file) = @_;
    my $round = $file_handler->{round};
    my $file_idx = $file_handler->{file_idx};
    
    my ($threads,$bowtie2_build,$bowtie2,$samtools) = ($self->{THREADS},PARTIES::Config->get_program_path('bowtie2-build'),PARTIES::Config->get_program_path('bowtie2'),PARTIES::Config->get_program_path('samtools'));
    my ($bwa) =PARTIES::Config->get_program_path('bwa');
  
    my $ies_fasta_file = $self->{IES_FASTA_FILE} ;    
    
    my $file_basename = $file_handler->{basename};
    
    my $bam_file = $self->{PATH}."/$file_basename.round$round.bam";
    my $nb_seq= `grep -c '^>' $seq_fasta_file`;
    chomp $nb_seq;
    $self->stderr("Process ".basename($seq_fasta_file)." (N=$nb_seq) ; ".basename($bam_file)." ...\n" );
    
    
    my $bam_bowtie_file = $self->{PATH}."/tmp/$file_basename.round$round.BOWTIE.bam";
    my $bowtie_unmapped_fasta = $self->{PATH}."/tmp/$file_basename.round$round.unmapped.fa";
    
    #~ #print STDERR "$bowtie2 -f --threads $threads --very-sensitive-local -L 10 --quiet  -x ".$ies_fasta_file." -U $seq_fasta_file | $samtools view -uS - 2> /dev/null | $samtools sort -o $bam_bowtie_file - > /dev/null 2>&1 && $samtools index $bam_bowtie_file \n";
    #~ #print STDERR join('', "$bowtie2 -f --threads $threads --very-sensitive-local -L 10 --quiet  -x ".$ies_fasta_file." -U $seq_fasta_file | $samtools view -uS - 2> /dev/null | $samtools sort -o $bam_bowtie_file - > /dev/null 2>&1 && $samtools index $bam_bowtie_file");
    
    system("$bowtie2 -f --threads $threads --very-sensitive-local -L 10 --quiet  -x ".$ies_fasta_file." -U $seq_fasta_file | $samtools view -uS - 2> /dev/null | $samtools sort -o $bam_bowtie_file - > /dev/null 2>&1 && $samtools index $bam_bowtie_file");
    #print STDERR "$samtools view -f 4 $bam_bowtie_file | awk '{ print \">\" \$1 \"\\n\" \$10 }' > $bowtie_unmapped_fasta\n";
    system("$samtools view -f 4 $bam_bowtie_file | awk '{ print \">\" \$1 \"\\n\" \$10 }' > $bowtie_unmapped_fasta");
    #print STDERR "$samtools view -b -F 4 $bam_bowtie_file  > $bam_bowtie_file.tmp && mv $bam_bowtie_file.tmp $bam_bowtie_file\n";
    system("$samtools view -b -F 4 $bam_bowtie_file  > $bam_bowtie_file.tmp && mv $bam_bowtie_file.tmp $bam_bowtie_file");
    
    my $bam_bwa_file =  $self->{PATH}."/tmp/$file_basename.round$round.BWA.bam";
    #print STDERR join('',"$bwa aln -t $threads  $ies_fasta_file $bowtie_unmapped_fasta > $bowtie_unmapped_fasta.sai 2> /dev/null && $bwa samse $ies_fasta_file $bowtie_unmapped_fasta.sai $bowtie_unmapped_fasta 2> /dev/null | $samtools view -uS - 2> /dev/null | $samtools sort -o $bam_bwa_file - > /dev/null 2>&1");
    system("$bwa aln -t $threads $ies_fasta_file $bowtie_unmapped_fasta > $bowtie_unmapped_fasta.sai 2> /dev/null && $bwa samse $ies_fasta_file $bowtie_unmapped_fasta.sai $bowtie_unmapped_fasta 2> /dev/null | $samtools view -uS - 2> /dev/null | $samtools sort -o $bam_bwa_file - > /dev/null 2>&1");
    
    #print STDERR "$samtools merge -f $bam_file $bam_bowtie_file $bam_bwa_file";
    system("$samtools merge -f $bam_file $bam_bowtie_file $bam_bwa_file");
    #print STDERR "$samtools index $bam_file\n";     
    
   # system("$samtools view -F 4 -q 10 -b $bam_file > t.bam && mv t.bam $bam_file");  
    
    system("$samtools index $bam_file");     
    #~ #system("rm $bam_bowtie_file $bam_bwa_file $bowtie_unmapped_fasta.sai $bowtie_unmapped_fasta");
    $self->stderr("Done\n" );
    
    my $nb_read_mapped=`samtools view -F 4 -c $bam_file`;
    chomp $nb_read_mapped;
    $self->stderr("Read $bam_file (N=$nb_read_mapped) ... \n" );    
    my $bam = Bio::DB::Sam->new(-bam  =>$bam_file, -fasta=>$ies_fasta_file)  ;
    my @unmapped_parts;
    foreach my $aln  ($bam->features(-type=>'match')) {
        #die $aln->seq_id;
        next if(!$aln);
        my $ies_id = $aln->seq_id;
        next if($ies_id eq '');
        next if($aln->get_tag_values('XM') > $self->{MAX_MISMATCH});
        # at least at one edge of the IES
        next if($aln->start !=1 and $aln->end != length($self->{IES_SEQ}->{$ies_id}));
        my $read_seq = $aln->query->dna;
        next if($read_seq=~/N/ or $aln->dna =~/N/);  
        #print STDERR $aln->qual,"\n";
        next if($aln->qual < $self->{MIN_MAPPING_QUALITY});
        my $rname=$aln->query->name;
        my ($part_start,$part_end) = (1,length($read_seq));
        if($rname=~/_part_(\d+)_(\d+)/) {
            ($part_start,$part_end) = ($1,$2);
            $rname=~s/_part_\S+//;
        } 
        
        my $cigar = $aln->cigar_str;  

        
        $self->{ALIGNMENTs}->{$rname}->{$file_idx} = {name =>$rname, seq=>$read_seq, has_a_partial_full_match=>0, aln_length=>0, ambiguous_matches=>'FALSE'} if(!$self->{ALIGNMENTs}->{$rname} or !$self->{ALIGNMENTs}->{$rname}->{$file_idx});
        
        my ($qstart,$qend) = ($aln->query->start+$part_start-1,$aln->query->end+$part_start-1);
        
        my $full_ies_match = ($aln->start==1 and $aln->end==length($self->{IES_SEQ}->{$ies_id})) ? 'TRUE' : 'FALSE';
        print STDERR join(" ","R$round",$aln->query->name,$ies_id,'l='.length($self->{IES_SEQ}->{$ies_id}),$aln->start, $aln->end,$cigar,$aln->strand,"full_ies_match=$full_ies_match", $qstart,$qend,'part',$part_start,$part_end),"\n"
if($rname eq 'NS500446:630:HLMV2AFXY:1:11201:10814:11874');
        # possible alignment types
        if($cigar =~ /^(\d+)M$/) {
            if($1 > $self->{MIN_MATCH_LENGTH}) {
                my $part = { type=>'M', tname =>$ies_id, tstart=>$aln->start, tend=>$aln->end , strand=>$aln->strand, qstart=> $qstart,qend=> $qend,seq=>$read_seq, full_ies_match=>$full_ies_match};
                
                my $ambiguous_match = (!$self->{REPEATED_BOUNDARIES}->{length($part->{seq})} or !$self->{REPEATED_BOUNDARIES}->{length($part->{seq})}->{$part->{seq}}) ? 'FALSE' : 'TRUE';
                $part->{ambiguous_match} = $ambiguous_match;
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{ambiguous_matches} = 'TRUE' if($ambiguous_match eq 'TRUE');
                
                push @{$self->{ALIGNMENTs}->{$rname}->{$file_idx}->{parts}}, $part;
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{has_a_partial_full_match} = 1;
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{aln_length} += length($read_seq);
                
                
            }
            
        } elsif($cigar =~ /^(\d+)M(\d+)S$/ and $aln->end == length($self->{IES_SEQ}->{$ies_id})){
            
            if($1 > $self->{MIN_MATCH_LENGTH}) {
                
                my $part = { type=>'MS', tname =>$ies_id, tstart=>$aln->start, tend=>$aln->end , strand=>$aln->strand , qstart=> $qstart,qend=> $qend, seq => substr($read_seq,0,$1-2) , full_ies_match=>$full_ies_match};
                
                my $ambiguous_match = (!$self->{REPEATED_BOUNDARIES}->{length($part->{seq})} or !$self->{REPEATED_BOUNDARIES}->{length($part->{seq})}->{$part->{seq}}) ? 'FALSE' : 'TRUE';
                $part->{ambiguous_match} = $ambiguous_match;
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{ambiguous_matches} = 'TRUE' if($ambiguous_match eq 'TRUE');
                
                push @{$self->{ALIGNMENTs}->{$rname}->{$file_idx}->{parts}}, $part;
                
                push @unmapped_parts , {id=>join("_",$rname,'part',$qend+1,length(substr($read_seq,$1-2))+$qend), seq => substr($read_seq,$1-2) };
                
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{aln_length} += length($part->{seq});
            }

        } elsif ($cigar =~ /^(\d+)S(\d+)M$/ and $aln->start ==1 ) {
                
            if($2 > $self->{MIN_MATCH_LENGTH}) {
                my $part = { type=>'SM', tname =>$ies_id, tstart=>$aln->start, tend=>$aln->end, strand=>$aln->strand, qstart=> $qstart,qend=> $qend, seq => substr($read_seq,$1+2) , full_ies_match=>$full_ies_match};
                
                my $ambiguous_match = (!$self->{REPEATED_BOUNDARIES}->{length($part->{seq})} or !$self->{REPEATED_BOUNDARIES}->{length($part->{seq})}->{$part->{seq}}) ? 'FALSE' : 'TRUE';
                $part->{ambiguous_match} = $ambiguous_match;
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{ambiguous_matches} = 'TRUE' if($ambiguous_match eq 'TRUE');
                
                push @{$self->{ALIGNMENTs}->{$rname}->{$file_idx}->{parts}}, $part;
                push @unmapped_parts , {id=>join("_",$rname,'part',$part_start ,$part_start+length(substr($read_seq,0,$1+2))), seq => substr($read_seq,0,$1+2) };
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{aln_length} += length($part->{seq});
            }

        } elsif ($cigar =~ /^(\d+)S(\d+)M(\d+)S$/ and $aln->start ==1 and $aln->end == length($self->{IES_SEQ}->{$ies_id})){
            
            if($2 > $self->{MIN_MATCH_LENGTH}) {
                my $part = { type=>'SMS',tname =>$ies_id, tstart=>$aln->start, tend=>$aln->end, strand=>$aln->strand , qstart=> $aln->query->start,qend=> $aln->query->end, seq => substr($read_seq,$1+2,$2-4), full_ies_match=>$full_ies_match };
                  
                my $ambiguous_match = (!$self->{REPEATED_BOUNDARIES}->{length($part->{seq})} or !$self->{REPEATED_BOUNDARIES}->{length($part->{seq})}->{$part->{seq}}) ? 'FALSE' : 'TRUE';
                $part->{ambiguous_match} = $ambiguous_match;
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{ambiguous_matches} = 'TRUE' if($ambiguous_match eq 'TRUE');
                
                push @{$self->{ALIGNMENTs}->{$rname}->{$file_idx}->{parts}}, $part;
                
                $self->{ALIGNMENTs}->{$rname}->{$file_idx}->{aln_length} += length($part->{seq});
                
                push @unmapped_parts , {id=>join("_",$rname,'part',1,length(substr($read_seq,0,$1+2))), seq => substr($read_seq,0,$1+2) };
                push @unmapped_parts , {id=>join("_",$rname,'part',$aln->query->end+1,length(substr($read_seq,$1+$2-2))+$aln->query->end), seq => substr($read_seq,$1+$2-2) };
            }
            
        }
            
    }
    $self->stderr("Done\n" ); 
    return \@unmapped_parts;
    
    
}

   
=head2 finish

 Usage   : $factory->finish
 Function: Finish all the procedure and/or comput all results
 Returns : Nothing
 Args    : Nothing

=cut

sub finish {
    my ($self,$results) = @_;  

    my @fastq_files=@{$self->{INPUT_FILEs}};

    foreach my $file_handler (@fastq_files) {
        my $file_basename = $file_handler->{basename};
        my $input_file = $file_handler->{input_file};

        my $unmapped_parts = $self->_mapping($file_handler,$input_file);


        while($unmapped_parts and scalar @{$unmapped_parts} != 0) {
            $file_handler->{round}++;
            my $round = $file_handler->{round};
            my $new_input_file= $self->{PATH}."/$file_basename.not_well_mapped.round$round.fa";

            # write new input file
            open(OUT,">$new_input_file") or die $new_input_file;
            foreach my $r (@{$unmapped_parts}) {
                print OUT ">",$r->{id},"\n",$r->{seq},"\n";
            }
            close OUT;

            $unmapped_parts =$self->_mapping($file_handler,$new_input_file);
        }
       # die;
    }
   
    my $mode = $self->get_mode;
    my $outfile = $self->{PATH}.'/'.$mode.".tab";
    $self->stderr("Write $outfile ... \n");
    open(OUT,">$outfile") or die $outfile;
    print OUT join("\t",qw(READ_NAME FILE_IDX READ_SEQ MONOMER AMBIGUOUS_MATCHES NB_IES IES_IDs IES_LENGTH DETAILS)),"\n";
    foreach my $rname (sort keys %{$self->{ALIGNMENTs}}) {
        foreach my $fh (sort keys %{$self->{ALIGNMENTs}->{$rname}}) {

            next if(!$self->{ALIGNMENTs}->{$rname}->{$fh}->{parts} or scalar @{$self->{ALIGNMENTs}->{$rname}->{$fh}->{parts}} == 1);
            my ($read_seq,$has_a_partial_full_match,$aln_length) = ($self->{ALIGNMENTs}->{$rname}->{$fh}->{seq} , $self->{ALIGNMENTs}->{$rname}->{$fh}->{has_a_partial_full_match},$self->{ALIGNMENTs}->{$rname}->{$fh}->{aln_length});
            next if(!$has_a_partial_full_match);
            my $read_length = length($read_seq);
            next if($read_length != $aln_length);
            #print STDERR join(" ",$rname,"R$fh",$read_seq,$read_length,$aln_length),"\n";
            my @ies_ids;
            my @details_matches;
            my @ies_lengths;
            my @parts = sort { $a->{qstart} <=> $b->{qstart} } @{$self->{ALIGNMENTs}->{$rname}->{$fh}->{parts}};
            
            my $uncoherent_concatemer= 0;
            if(scalar @parts > 2) {
                for(my $i=1;$i<(scalar(@parts) -1);$i++) {
                    $uncoherent_concatemer = 1 if($parts[$i]->{full_ies_match} eq 'FALSE');
                }
            }
            next if($uncoherent_concatemer);
                
            foreach my $part (@parts) {

                my ($aln_type,$ies_id,$tstart,$tend,$qstart,$qend,$match_seq,$ambiguous_match,$full_ies_match) = ($part->{type},$part->{tname},$part->{tstart},$part->{tend},$part->{qstart},$part->{qend},$part->{seq},$part->{ambiguous_match},$part->{full_ies_match});
                my $strand = ($part->{strand} > 0) ? '+' : '-';
                push @ies_ids, $ies_id;
                push @ies_lengths, length($self->{IES_SEQ}->{$ies_id});
                push @details_matches, join(":", $tstart,$tend,$strand,$ambiguous_match);
                
                #$unique_ies_ids{$ies_id}=1;
            }
            my $is_monomer = (scalar(@{$self->{ALIGNMENTs}->{$rname}->{$fh}->{parts}}) == 2 and scalar uniq(@ies_ids) ==1) ? 'TRUE' : 'FALSE';

            print OUT join("\t",$rname,"R$fh",$read_seq,$is_monomer,$self->{ALIGNMENTs}->{$rname}->{$fh}->{ambiguous_matches},scalar(@{$self->{ALIGNMENTs}->{$rname}->{$fh}->{parts}}),join(" ",@ies_ids),join(" ",@ies_lengths),join(" ",@details_matches)),"\n";
            #die ;
        }
    }
    close OUT;
        
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
  # nothing to do in paralel
  return;

}


1;
