package Bio::Graphics::Browser2::Plugin::BlastSearch;

use strict;
use warnings;
use Bio::Graphics::Browser2::Plugin;
use Bio::Graphics::Feature;
use File::Spec::Unix;
use Config::Tiny;
use CGI qw(:standard *table);
use vars '$VERSION','@ISA','$blast_executable';
use Config::Tiny;
use File::Spec::Unix;

=head1 NAME

BIO::Graphics::Browser2::Plugin::BlastSearch -- a plugin that uses BLAST to search a given sequence against the scaffold

=head1 SYNOPSIS


=head1 DESCRIPTION

This Gbrowse plugin will take a sequence (entered in teh configuration screen)
and BLAST it against the full scaffold, creating hits as a new sequence.

=cut

# GBrowse plugins have several subroutines that are called automatically by the browser. Which subroutines are called depend on the type of plugin.
# Every plugin type calls config_defaults when the plugin first runs, configure_form when the user configures the plugin, and config_defaults when the user saves the configuration
# Finders call the find subroutine and accept a list of objects implementing the featurel interface.
# More information about plugins is available at: http://search.cpan.org/~lds/GBrowse-2.31/lib/Bio/Graphics/Browser2/Plugin.pm

$blast_executable = "";

$VERSION = '1.00';

@ISA = qw(Bio::Graphics::Browser2::Plugin);

my @COLORS = qw(red green blue orange cyan black
	turquoise brown indigo wheat yellow emerald);

sub name { "Alignments With BLAST" }

sub description {
	p("This plugin will take an input sequence - entered in the 'Configure' screen - and run blast to find hits,",
	"hits will be drawn as features").
	p("This plugin was written by Brad Covey.");
}

sub type { 'finder' }
sub init {
	my $self = shift;
	my $conf = $self->browser_config;

}


sub config_defaults {
	my $cfg = Config::Tiny->new;
	$cfg = Config::Tiny->read('/etc/gbrowse_addons.ini');
	my $blastdb = $cfg->{main}{blastdb_location};
	$blastdb = File::Spec->catfile($blastdb, "Toceanica2");
	my $self = shift;
	return {sequence_to_blast => '',
			"job_name" => 'BlastSearch',
			"job_description" => 'BlastSearch',
			"file_upload" => 'FASTA file',
			"query_type" => 'blastn',
			task => 'blastn',
			db => "$blastdb",
			evalue => '',
			"word_size" => '',
			"culling_limit" => '',
			"max_target_seqs" => ''}
}

sub reconfigure {
	my $self = shift;
	my $current = $self->configuration;

	my $filename = $self->config_param('file_upload');
	if( $filename eq 'FASTA file' || $filename eq '') {
		$current->{'sequence_to_blast'}    = $self->config_param('sequence_to_blast');
	} else {
		use File::Temp;
		my ($fh_in, $filename_in) = tempfile();
		my ($bytesread, $buffer);
		my $numbytes = 1024;
	
		while($bytesread = read($filename, $buffer, $numbytes)) {
			print $fh_in $buffer;
		}
		$current->{'file_upload_filename'} = $filename_in;
	}
	
	$current->{'query_type'}      = $self->config_param('query_type');
	$current->{'task'}            = $self->config_param('task');
	$current->{'evalue'}          = $self->config_param('evalue');
	$current->{'word_size'}       = $self->config_param('word_size');
	$current->{'culling_limit'}   = $self->config_param('culling_limit');
	$current->{'max_target_seqs'} = $self->config_param('max_target_seqs');

	my $track_name = $self->config_param('job_name');
	$track_name =~ s/ /-/g;
	unless( $track_name eq '') {
		$current->{'job_name'} = $track_name;
	}
	$current->{'job_description'} = $self->config_param('job_description');
}

sub configure_form {
	my $self = shift;
	my $current_config = $self->configuration;
	my $form;

	my %labels = (
		'blastn' => 'Nucleotide Sequence',
		'tblastn' => 'Protein Sequence'
	);

	$form = h3("Leave fields blank for default BLAST values").
		table({-border => 0},TR([
			td([b("Sequence To Align")]),
			td(["Enter a Sequence: ",	textfield(-name => $self->config_name('sequence_to_blast'), -size=> 100, -value=>$current_config->{'sequence_to_blast'})]),
			td(["<div style=\"text-align:center;\">or</div>"]),
			td(["FASTA File Upload: ", filefield(-name=> $self->config_name('file_upload'), -default=>'FASTA file',-size=>50,-maxlength=>150)]),
			td("Query Type: ").td(radio_group( -name=>$self->config_name('query_type'), -values=>['blastn','tblastn'], -default=>$current_config->{'query_type'},-labels=>\%labels)),
			td(["<br><br>"]),
			td([b("Job Name")]),
			td(["Name(optional): ", textfield(-name => $self->config_name('job_name'), -default=>'', -size=>20, -maxlength=>150, -value=>$current_config->{'job_name'})]),
			td(["Description(optional, used as searchable tag): ", textfield(-name => $self->config_name('job_description'), -defaults=>'', -size=>20, -maxlength=>150, -value=>$current_config->{'job_description'})]),
			td(["<br>"]),
			td([b("Advanced Settings")]),
			td("Leave blank for BLAST defaults"),
			td(["Task(Nucleotide Sequences Only): ",	popup_menu($self->config_name('task'),['blastn','blastn-short','dc-megablast','megablast','rmblastn'], $current_config->{'task'})]),
			td(["evalue: ",	textfield(-name=>$self->config_name('evalue'), -default=>$current_config->{'evalue'},-size=>10,-maxlength=>10)]),
			td(["Word Size: ", textfield(-name=>$self->config_name('word_size'), -default=>$current_config->{'word_size'},-size=>5,-maxlength=>5,-title=>'Word size for wordfinder algorithm (length of best perfect match)')]),
			td(["Culling Limit: ", textfield(-name=>$self->config_name('culling_limit'), -default=>$current_config->{'culling_limit'},-size=>5,-maxlength=>5,-title=>'If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit')]),
			td(["Max Target Sequences: ", textfield(-name=>$self->config_name('max_target_seqs'), -default=>$current_config->{'max_target_seqs'},-size=>7,-maxlength=>7,-title=>'Maximum number of aligned sequences to keep. Not applicable for outfmt <= 4')])
			]));

	return $form;
}

sub find {
	my $self = shift;
	my $segment  = shift;  # we ignore this!
	my $config   = $self->configuration;
	my @results;

	my $bl_result_file = $self->do_blast() || die "blast call failed";
	use Bio::SearchIO;
	my $searchio = new Bio::SearchIO(-format => 'blasttable',
										-file => $bl_result_file);
	# while( my $result = $searchio->next_result) {
	# 	while( my $hit = $result->next_hit ) {
	# 		while( my $hsp = $hit->next_hsp ) {
	# 			# collect feature parameters
	# 			my $seq_id = $hit->accession;
	# 			$seq_id = substr $seq_id,4;
	# 			my $start = $hsp->start('hit');
	# 			my $end = $hsp->end('hit');
	# 			my $name = $self->configuration->{'job_name'};
	# 			my $type = 'BlastSearch';
	# 			my $score = $hsp->score;
	# 			my $desc = 'Evalue: '.$hsp->evalue.' Hit Strand: '.$hsp->strand('hit').' Query Strand: '.$hsp->strand('query');
	# 			my $strand = $hsp->strand('hit');
	# 			# create feature and push onto list to be returned
 #                print STDERR "SEQID: $seq_id\n";
	# 			push @results, Bio::Graphics::Feature->new(
	# 				-seq_id => $seq_id,
	# 				-start => $start,
	# 				-end => $end,
	# 				-name => $name,
	# 				-type => $type,
	# 				-score => $score,
	# 				-desc => $desc,
	# 				-strand => "$strand");
	# 		}
	# 	}
	while( my $result = $searchio->next_result) {
		while( my $hit = $result->next_hit ) {
			while( my $hsp = $hit->next_hsp ) {
				# collect feature parameters
				my $seq_id = $hit->name;
				my @split = split(/\|/,$seq_id);
				$seq_id = $split[3];
				print STDERR "SEQ ID: $seq_id\n";
				# $seq_id = substr $seq_id,4;
				my $start = $hsp->start('hit');
				my $end = $hsp->end('hit');
				my $name = $self->configuration->{'job_name'};
				my $type = 'BlastSearch';
				my $score = $hsp->score;
				my $desc = $self->configuration->{'job_description'}.' Evalue: '.$hsp->evalue.' Hit Strand: '.$hsp->strand('hit').' Query Strand: '.$hsp->strand('query');
				my $strand = $hsp->strand('hit');
				# create feature and push onto list to be returned
				push @results, Bio::Graphics::Feature->new(
					-seq_id => $seq_id,
					-start => $start,
					-end => $end,
					-name => $name,
					-type => $type,
					-score => $score,
					-desc => $desc,
					-strand => "$strand");
			}
		}
        # unlink $bl_result_file;
		# return the array of features
		use Data::Dumper;
		print STDERR "RETURNING: ".Dumper(\@results);
		return (\@results, "BLAST Search");
	}	
}

# subroutine that performs a blast based on the configuration parameters
# returns the results of the blast in a tempfile
sub do_blast {
	my ($self) = @_;
	use File::Temp;
	use File::Temp qw/ tempfile tempdir /;
	# create tempfile and store the query in it
	my ($query, $fh_query, $filename_query);
	if( $self->config_param('file_upload') eq 'FASTA file' || $self->config_param('file_upload') eq '') {
		($fh_query, $filename_query) = tempfile();
		my $query_seq = $self->configuration->{'sequence_to_blast'};
		print $fh_query ">Query Sequence\n$query_seq\n";
		$query = "-query $filename_query";
	} else {
		$filename_query = $self->configuration->{'file_upload_filename'};
		$query = "-query $filename_query";
	}


	# fetch configuration parameters
	my $task            = $self->configuration->{'task'};
	my $evalue          = $self->configuration->{'evalue'};
	my $db              = $self->configuration->{'db'};
	my $conf_file = "/etc/gbrowse_addons.ini";
	my $cfg = Config::Tiny->new;
	$cfg = Config::Tiny->read($conf_file);
    $db = File::Spec->catfile($cfg->{main}{blastdb_location}, $self->page_settings->{'source'});
	my $word_size       = $self->configuration->{'word_size'};
	my $culling_limit   = $self->configuration->{'culling_limit'};
	my $max_target_seqs = $self->configuration->{'max_target_seqs'};
	my $query_type      = $self->configuration->{'query_type'};

	# fetch a flag for query_type
	my $is_nucl = $query_type eq 'blastn' ? 1 : 0;

	# format parameters to be passed to BLAST leaving unused parameters as empty strings
	if($task && $is_nucl){
		$task = "-task $task";
	} else {
		$task = "";
	}

	if($evalue) {
		$evalue = "-evalue $evalue";
	} else {
		$evalue = "";
	}

	if($word_size) {
		$word_size = "-word_size $word_size";
	} else {
		$word_size = "";
	}

	if($culling_limit) {
		$culling_limit = "-culling_limit $culling_limit";
	} else {
		$culling_limit = "";
	}

	if($max_target_seqs) {
		$max_target_seqs = "-max_target_seqs $max_target_seqs";
	} else {
		$max_target_seqs = "";
	}

	# call blast with the query tempfile and the configuration parameters
	# my $blast_result = `blastn -query $filename_query -task $task -db $db -evalue $evalue`;
	my ($fh_results, $filename_results) = tempfile();
	my $blast_result = `$query_type $query $task -outfmt 6 -db $db $evalue $word_size $culling_limit $max_target_seqs -out $filename_results`;
	print STDERR "BLAST OUTPUT: $filename_results\n";
	
	# store blast results in a tempfile, unlink query file and return results file
	# print $fh_results $blast_result;
	unlink $filename_query;
	return $filename_results;
}

1;
