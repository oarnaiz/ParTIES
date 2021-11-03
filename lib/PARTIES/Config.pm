package PARTIES::Config;

use strict;


use File::Which qw(which);
use Bio::GFF3::LowLevel qw/ gff3_parse_feature /;
use Statistics::R;
use Bio::SeqIO;
use Bio::DB::Sam;


use PARTIES::Map;
use PARTIES::MIRAA;
use PARTIES::MICA;
use PARTIES::Insert;
use PARTIES::Assembly;
use PARTIES::Run;
use PARTIES::MIRET;
use PARTIES::MILORD;
use PARTIES::Compare;
use PARTIES::Concatemer;
use PARTIES::Utils;



=head2 get_parties_modes

 Title   : Get all availables modes
 Usage   : get_parties_modes;
 Function: 
 Returns : Hash with all modes
 Args    : Nothing

=cut

sub get_parties_modes {
   my %modes = ( 
		Run 		=> { N=> 1, DESC=>'Run ParTIES using the configuration file'},
		Map 		=> { N=> 2, DESC=>'Map reads on a reference using bowtie2'},
		MIRAA		=> { N=> 3, DESC=>'Method of Identification by Read Alignment Anomalies'},
		Assembly 	=> { N=> 6, DESC=>'Filter reads and assemble them'},
		MICA 		=> { N=> 4, DESC=>'Method of Identification by Comparison of Assemblies'},
		Insert 		=> { N=> 5, DESC=>'Insert IES within a genome to create an IES containing reference'},		
		MIRET 		=> { N=> 6, DESC=>'Method of Ies RETention'},
		MILORD 		=> { N=> 7, DESC=>'Method of Identification and Localization of Rare Deletion'},
		Compare 	=> { N=> 8, DESC=>'Compare IES/InDel datasets'},
		Concatemer 	    => { N=> 9, DESC=>'IES Concatemer detection'},
               );
   return %modes;
} 


=head2 get_program_path

 Title   : Get program paths
 Usage   : get_program_path;
 Function: 
 Returns : path
 Args    : program name

=cut


sub get_program_path {
  my ($self,$program) = @_;  
  die "Need program" if(!$program);
    
  # change here 
  my %program_paths = (
			velveth => '',
			velvetg => '',
			gzip => '',
			bowtie2 => '',
			samtools => '',
			RepeatMasker => '',
			blat => '',
			muscle => '',
			bwa => ''
			);
  
  
  if($program_paths{$program} and -e $program_paths{$program}) {
     return $program_paths{$program};
  } else {
    my ($path) = which($program);
    if($path) {
       return $path;
    }   
  } 
  die "$program seems to be not installed";
}


1;
