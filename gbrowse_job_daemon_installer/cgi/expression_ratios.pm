#!/usr/bin/perl
package expression_ratios;

# FIX RECIPROCAL BLAST
# This module does the heavy lifting in the expression ratio calculation
# it has different modes for the different types of expression calculation

# TODO:
# 	implement control vs experimental with multiple input files
# 	implement DNA vs RNA with multiple input files

use warnings;
use strict;
use File::Temp qw/ tempfile /;
use File::Basename;
use POSIX;
use Data::Dumper;
use Bio::SeqIO;
use Exporter;
use Config::Tiny;

# boilerplate for exporting functions
our @ISA = qw( Exporter );

our @EXPORTER_OK = qw(calculate_expression);

our @EXPORT = qw(calculate_expression);

my $working_dir = "/usr/share/gbrowse_job_daemon/";
chdir $working_dir;

my $conf_file = "/etc/gbrowse_addons.ini";
my $cfg = Config::Tiny->new;
$cfg = Config::Tiny->read($conf_file);

my $store_location = File::Spec->catdir($cfg->{main}{gbrowse_job_daemon_location}, "expression");

# how to call this function:
# calculate_expression($reference, 'cve', $experimental, $control)
# calculate_expression($reference, 'dvr', $dna, $rna)
# calculate_expression($reference, 'dor', $file )
# calculate_expression($reference, 'dorm', $file1, $file2, $file3, $file4...)
# cve = control vs. experimental
# dvr = DNA vs. RNA
# dor = DNA or RNA, single sample
# dorm = DNA or RNA, multiple sample
sub calculate_expression {
	my $title     = shift;
	my $reference = shift;
	my $confid	  = shift;
	my $task      = shift;
	my $query1    = shift;
	my $query2    = shift;
	my @inputs;
	while( my $tmp = shift ) {
		push(@inputs, $tmp);
	}


	if( $task eq "cve" ) {
		# remove sequences that don't hit the reference
		# potentially perform alignment 
		# my $experimental = trim_input( $query1, $reference, 0 );
		# my $control      = trim_input( $query2, $reference, 0 );

		# get hits
		my $experimental_blast = blast($query1, $reference);
		my $control_blast      = blast($query2, $reference);

		# average hte hits from the previous step over 300bp windows on the reference genome
		my %experimental_hash = %{average_over_windows($experimental_blast)};
		my %control_hash      = %{average_over_windows($control_blast)};

		# find the total number of sequences(after trimming) found in the experimental and the control
		my $experimental_seq  = Bio::SeqIO->new(-file => $query1,
		                           		-format => 'Fasta'
		);
		my $experimental_length = 0;
		while( my $seq = $experimental_seq->next_seq ) {
			$experimental_length += 1;
		}

		my $control_seq  = Bio::SeqIO->new(-file => $query2,
		                           		-format => 'Fasta'
		);
		my $control_length = 0;
		while( my $seq = $control_seq->next_seq ) {
			$control_length += 1;
		}

		# Step through the hashes and merge the averaged results, normalize using the sequence counts
		# Also builds a hash containing the total number of hits in each reference
		my %hit_count;
		my %variance;
		my %merged_hash;
		foreach my $ref_key ( keys %experimental_hash ) {
			my $sum = 0;
			my $var_max = 0;
			foreach my $key ( keys %{$experimental_hash{$ref_key}} ) {
				if( !(exists $control_hash{$ref_key}{$key}) ) {
					$merged_hash{$ref_key}{$key} = log(1)/log(2);
					delete $experimental_hash{$ref_key}{$key};
				} else {
					# dna_length and rna_length are equal to the total number of reads that mapped to any scaffold in the initial BLAST search
					$merged_hash{$ref_key}{$key} = log((($experimental_hash{$ref_key}{$key}/300)/($control_hash{$ref_key}{$key}/300))*($control_length/$experimental_length))/log(2);
					$sum += ($experimental_hash{$ref_key}{$key} + $control_hash{$ref_key}{$key});
					if( $var_max < abs($merged_hash{$ref_key}{$key}) ) {
						$var_max = abs($merged_hash{$ref_key}{$key});
					}
					delete $experimental_hash{$ref_key}{$key};
					delete $control_hash{$ref_key}{$key};
				}
			}
			if( $sum > 0 ) {
				$hit_count{$ref_key} += $sum;
			}
			if( $var_max > 0 ) {
				$variance{$ref_key} += $var_max;
			}
		}
		foreach my $ref_key ( keys %control_hash ) {
			my $sum = 0;
			my $var_max = 0;
			foreach my $key ( keys %{$control_hash{$ref_key}} ) {
				if( !(exists $experimental_hash{$ref_key}{$key}) ) {
					$merged_hash{$ref_key}{$key} = log(1)/log(2);
					delete $control_hash{$ref_key}{$key};
				} else {
					# dna_length and rna_length are equal to the total number of reads that mapped to any scaffold in the initial BLAST search
					$merged_hash{$ref_key}{$key} = log((($experimental_hash{$ref_key}{$key}/300)/($control_hash{$ref_key}{$key}/300))*($control_length/$experimental_length))/log(2);
					$sum += ($experimental_hash{$ref_key}{$key} + $control_hash{$ref_key}{$key});
					if( $var_max < abs($merged_hash{$ref_key}{$key}) ) {
						$var_max = abs($merged_hash{$ref_key}{$key});
					}
					delete $experimental_hash{$ref_key}{$key};
					delete $control_hash{$ref_key}{$key};
				}
			}
			if( $sum > 0 ) {
				$hit_count{$ref_key} = $sum;
			}
			if( $var_max > 0 ) {
				$variance{$ref_key} = $var_max;
			}
		}
		open my $hit_count_fh, ">", File::Spec->catfile($store_location, $confid."_hits") or die "Unable to open ".File::Spec->catfile($store_location, $confid."_hits")." for writing: $!";
		foreach my $val (reverse sort {$hit_count{$a} <=> $hit_count{$b} } keys %hit_count) {
			print $hit_count_fh $val."\t".$hit_count{$val}."\n";
		}
		close $hit_count_fh;
		open my $var_fh, ">", File::Spec->catfile($store_location, $confid."_variance") or die "Unable to open file for writing: $!";
		foreach my $val (reverse sort {$variance{$a} <=> $variance{$b} } keys %variance) {
			print $var_fh $val."\t".$variance{$val}."\n";
		}
		close $var_fh;

		# build the hashes into a feature file to be uploaded to gbrowse. The file returned includes a simple configuration for drawing as an XYPlot
		my $feature_file = generate_featurefile(\%merged_hash, $reference, $title, $query1, $reference, $confid);
		return $feature_file, File::Spec->catfile($store_location, $confid."_hits"), File::Spec->catfile($store_location, $confid."_variance");
	} elsif( $task eq "dvr" ) {
		# remove sequences that don't hit the reference
		# potentially perform alignment 
		# my $dna = trim_input( $query1, $reference, 0 );
		# my $rna = trim_input( $query2, $reference, 0 );

		# get hits
		my $dna_blast = blast($query1, $reference);
		my $rna_blast = blast($query2, $reference);

		# average hte hits from the previous step over 300bp windows on the reference genome
		my %dna_hash = %{average_over_windows($dna_blast)};
		my %rna_hash = %{average_over_windows($rna_blast)};

		# find the total number of sequences(after trimming) found in the experimental and the control
		my $dna_seq  = Bio::SeqIO->new(-file => $query1,
		                           		-format => 'Fasta'
		);
		my $dna_length = 0;
		while( my $seq = $dna_seq->next_seq ) {
			$dna_length += 1;
		}

		my $rna_seq  = Bio::SeqIO->new(-file => $query2,
		                           		-format => 'Fasta'
		);
		my $rna_length = 0;
		while( my $seq = $rna_seq->next_seq ) {
			$rna_length += 1;
		}

		# step through the hashes and merge the averaged results, normalizing with the sequence count
		# Also builds a hash containing the total number of hits in each reference
		my %hit_count;
		my %merged_hash;
		my %variance;
		foreach my $ref_key ( keys %dna_hash ) {
			my $sum = 0;
			my $var_max = 0;
			foreach my $key ( keys %{$dna_hash{$ref_key}} ) {
				if( !(exists $rna_hash{$ref_key}{$key}) ) {
					$merged_hash{$ref_key}{$key} = log(1)/log(2);
					delete $dna_hash{$ref_key}{$key};
				} else {
					# dna_length and rna_length are equal to the total number of reads that mapped to any scaffold in the initial BLAST search
					$merged_hash{$ref_key}{$key} = log((($rna_hash{$ref_key}{$key}/300)/($dna_hash{$ref_key}{$key}/300))*($dna_length/$rna_length))/log(2);
					$sum += $dna_hash{$ref_key}{$key} + $rna_hash{$ref_key}{$key};
					if( $var_max < abs($merged_hash{$ref_key}{$key}) ) {
						$var_max = abs($merged_hash{$ref_key}{$key});
					}
					delete $dna_hash{$ref_key}{$key};
					delete $rna_hash{$ref_key}{$key};
				}
			}
			$hit_count{$ref_key} = $sum;
			$variance{$ref_key} = $var_max;
		}
		foreach my $ref_key ( keys %rna_hash ) {
			my $sum = 0;
			my $var_max = 0;
			foreach my $key ( keys %{$rna_hash{$ref_key}} ) {
				if( !(exists $dna_hash{$ref_key}{$key}) ) {
					$merged_hash{$ref_key}{$key} = log(1)/log(2);
					delete $rna_hash{$ref_key}{$key};
				} else {
					# dna_length and rna_length are equal to the total number of reads that mapped to any scaffold in the initial BLAST search
					$merged_hash{$ref_key}{$key} = log((($rna_hash{$ref_key}{$key}/300)/($dna_hash{$ref_key}{$key}/300))*($dna_length/$rna_length))/log(2);
					$sum += $dna_hash{$ref_key}{$key} + $rna_hash{$ref_key}{$key};
					if( $var_max < abs($merged_hash{$ref_key}{$key}) ) {
						$var_max = abs($merged_hash{$ref_key}{$key});
					}
					delete $dna_hash{$ref_key}{$key};
					delete $rna_hash{$ref_key}{$key};
				}
			}
			$hit_count{$ref_key} = $sum;
			$variance{$ref_key} = $var_max;
		}
		open my $hit_count_fh, ">", File::Spec->catfile($store_location, $confid."_hits") or die "Unable to open file for writing: $!";
		foreach my $val (reverse sort {$hit_count{$a} <=> $hit_count{$b} } keys %hit_count) {
			print $hit_count_fh $val."\t".$hit_count{$val}."\n";
		}
		open my $var_fh, ">", File::Spec->catfile($store_location, $confid."_variance") or die "Unable to open file for writing: $!";
		foreach my $val (reverse sort {$variance{$a} <=> $variance{$b} } keys %variance) {
			print $var_fh $val."\t".$variance{$val}."\n";
		}
		close $var_fh;

		# build the hashes into a feature file to be uploaded to gbrowse. The file returned includes a simple configuration for drawing as an XYPlot
		my $feature_file = generate_featurefile(\%merged_hash, $reference, $title, $query1, $reference, $confid);
		return $feature_file, File::Spec->catfile($store_location, $confid."_hits"), File::Spec->catfile($store_location, $confid."_variance");
	} elsif( $task eq "dor" ) {
		# remove sequences that don't hit the reference
		# potentially perform alignment 
		# my $input = trim_input( $query1, $reference, 0 );
		# get hits
		my $input_blast = blast($query1, $reference);
		# average hte hits from the previous step over 300bp windows on the reference genome
		# Also builds a hash containing the total number of hits in each reference
		my %input_hash = %{average_over_windows($input_blast)};
		my %hit_count;
		my %variance;
		foreach my $ref_key ( keys %input_hash ) {
			my $sum = 0;
			my $var_max = 0;
			foreach my $key ( keys %{$input_hash{$ref_key}} ) {
				$input_hash{$ref_key}{$key} = $input_hash{$ref_key}{$key}/300;
				$sum += $input_hash{$ref_key}{$key};
				if( $var_max < abs($input_hash{$ref_key}{$key}) ) {
					$var_max = abs($input_hash{$ref_key}{$key});
				}
			}
			$hit_count{$ref_key} = $sum;
		}
		open my $hit_count_fh, ">", File::Spec->catfile($store_location, $confid."_hits") or die "Unable to open file for writing: $!";
		foreach my $val (reverse sort {$hit_count{$a} <=> $hit_count{$b} } keys %hit_count) {
			print $hit_count_fh $val."\t".$hit_count{$val}."\n";
		}
		open my $var_fh, ">", File::Spec->catfile($store_location, $confid."_variance") or die "Unable to open file for writing: $!";
		foreach my $val (reverse sort {$variance{$a} <=> $variance{$b} } keys %variance) {
			print $var_fh $val."\t".$variance{$val}."\n";
		}
		close $var_fh;

		# build the hashes into a feature file to be uploaded to gbrowse. The file returned includes a simple configuration for drawing as an XYPlot
		my $feature_file = generate_featurefile(\%input_hash, $reference, $title, $query1, $reference, $confid);
		return $feature_file, File::Spec->catfile($store_location, $confid."_hits"), File::Spec->catfile($store_location, $confid."_variance");
	} elsif( $task eq "dorm" ) {
		push(@inputs, $query1);
		push(@inputs, $query2);

		my ($output_features_fh, $output_features_filename) = tempfile( DIR => $working_dir."output_files/");

		my $count = 0;
		my %hit_count;
		my %variance;
		foreach my $file (@inputs) {
			# remove sequences that don't hit the reference
			# potentially perform alignment 
			# my $input = trim_input( $file, $reference, 0 );
			# get hits
			my $input_blast = blast($file, $reference);
			# average hte hits from the previous step over 300bp windows on the reference genome
			my %input_hash = %{average_over_windows($input_blast)};
			foreach my $ref_key ( keys %input_hash ) {
				my $sum = 0;
				my $var_max = 0;
				foreach my $key ( keys %{$input_hash{$ref_key}} ) {
					$sum += $input_hash{$ref_key}{$key};
					$input_hash{$ref_key}{$key} = $input_hash{$ref_key}{$key}/300;
					if( $var_max < abs($input_hash{$ref_key}{$key}) ) {
						$var_max = abs($input_hash{$ref_key}{$key});
					}
				}
				$hit_count{$ref_key} += $sum;
				if( $hit_count{$ref_key} < $var_max ) {
					$hit_count{$ref_key} = $var_max;
				}
			}
			# build the hashes into a feature file to be uploaded to gbrowse. The file returned includes a simple configuration for drawing as an XYPlot
			my $feature_file = generate_featurefile(\%input_hash, $reference, "file-$count", $query1, $reference, $confid);

			open( FILE, "<", $feature_file );
			while (<FILE>) {
				print $output_features_fh $_
			}
			print $output_features_fh "\n";
			$count++;
		}
		open my $hit_count_fh, ">", File::Spec->catfile($store_location, $confid."_hits") or die "Unable to open file for writing: $!";
		foreach my $val (reverse sort {$hit_count{$a} <=> $hit_count{$b} } keys %hit_count) {
			print $hit_count_fh $val."\t".$hit_count{$val}."\n";
		}
		close $hit_count_fh;
		open my $var_fh, ">", File::Spec->catfile($store_location, $confid."_variance") or die "Unable to open file for writing: $!";
		foreach my $val (reverse sort {$variance{$a} <=> $variance{$b} } keys %variance) {
			print $var_fh $val."\t".$variance{$val}."\n";
		}
		close $var_fh;

		return $output_features_filename, File::Spec->catfile($store_location, $confid."_hits"), File::Spec->catfile($store_location, $confid."_variance");
	} else {
		die "Invalid task";
	}

	# perform a BLAST, store the results in a tempfile, return the tempfile
	sub blast {
		my $query     = shift;
		my $reference = shift;

		my ($blast_fh, $blast_filename) = tempfile ( DIR => $working_dir."temp_results/" );
		my $result = `blastn -query $query -db $reference -max_target_seqs 1 -outfmt 6 -out $blast_filename`;
		close $blast_filename;

		return $blast_filename;
	}

	# trim_input( $query_filename, $db_filename, $reciprocal_blast ); returns trimmed query filename
	sub trim_input {
		use File::Temp;
		use File::Basename;
		use Bio::SearchIO;

		my $query_filename   = shift;
		my $db_filename      = shift;
		my $reciprocal_blast = shift;

		my $query_title = basename($query_filename);

		open my $in, $query_filename or die "Could not open $query_filename: $!";
		my( $query_db_fh, $query_db_filename ) = tempfile(DIR => $working_dir."temp_blastdb/");

		while( my $line = <$in>)  {
		    print $query_db_fh $line;
		}
		close $in;
		close $query_db_filename;

		my $result = `makeblastdb -dbtype nucl -parse_seqids -title $query_title  -in $query_db_filename -out $query_db_filename`;
		# Can't easily check for errors from makeblastdb because it returns a lot of output on success as well
		# Figure something out if this becomes an issue

		my ($blast1_fh, $blast1_filename) = tempfile( DIR => $working_dir."temp_results/");
		my ($blast2_fh, $blast2_filename) = tempfile( DIR => $working_dir."temp_results/");
		close $blast1_filename;
		close $blast2_filename;

		# blast the query against the scaffold
		$result = `blastn -query $query_db_filename -db $db_filename -outfmt 6 -out $blast1_filename`;
		if($result ne '') {
			die "Error from blast: $result";
		}

		my %unique_hits;
		# FIX ME
		if($reciprocal_blast) {
			# blast the scaffold against the query
			$result = `blastn -query $db_filename -db $query_db_filename -outfmt 6 -out $blast2_filename`;
			if($result ne '') {
				die "Error from blast: $result";
			}

			my $parser2 = new Bio::SearchIO(
				-format => 'blasttable',
				-file   => $blast2_filename
			);

			while ( my $result = $parser2->next_result() ) {
			    while( my $hit = $result->next_hit ) {
					while( my $hsp = $hit->next_hsp ) {
						$unique_hits{$hit->name()} = 1;
					}
			    }
			}
		}

		# Create a hash that will be used to chop sequences without hits out of the query fasta
		my $parser1 = new Bio::SearchIO(
			-format => 'blasttable',
			-file   => $blast1_filename
		);
		print STDERR "BLAST RESULT: $blast1_filename\n";
		while ( my $result = $parser1->next_result() ) {
		    while( my $hit = $result->next_hit ) {
				while( my $hsp = $hit->next_hsp ) {
					$unique_hits{$result->query_name} = 1;
				}
		    }
		}

		# convert the hash into a string
		my $query_string = '';
		foreach my $key ( keys %unique_hits )
		{
			$query_string .= "$key,";
		}
		$query_string = substr $query_string, 0, -1;

		# call blastdbcmd to chop the sequences out of the fasta
		chdir $working_dir."temp_blastdb/";
		my $db_name = basename($query_db_filename);
		my ($unique_fasta_fh, $unique_fasta_filename) = tempfile( DIR => $working_dir."temp_blastdb/" );
		$result = `blastdbcmd -db $db_name -entry $query_string -out $unique_fasta_filename`;
		if($result ne '') {
			die "Error from blastdbcmd: $result";
		}

		return $unique_fasta_filename;
	}

	sub average_over_windows {
		my $blast_filename = shift;

		my $parser = new Bio::SearchIO(
			-format => 'blasttable',
			-file   => $blast_filename
		);

		my %window_hash;
		while ( my $result = $parser->next_result() ) {
		    while( my $hit = $result->next_hit ) {
				while( my $hsp = $hit->next_hsp ) {
					# initializing new parts of the hash
					if( !(exists $window_hash{$hit->name}) ) {
						# figure out where the window starts and ends
						my $start_id = floor($hsp->start('hit')/300);
						my $end_id = ceil($hsp->end('hit')/300);
						my $i;
						# count the total # of bp within the current window covered by the current hit, add to the total
						for( $i = $start_id; $i <= $end_id; $i++ ) {
							my $bp_in_window;
							if( $hsp->start('hit') >= ($i*300+1) && $hsp->end('hit') >= (($i+1)*300) ) {
								$bp_in_window = ($i+1)*300 - $hsp->start('hit');
							} elsif( $hsp->end('hit') <= (($i+1)*300 ) && $hsp->start('hit') <= ($i*300+1) && $hsp->end('hit') >= ($i*300+1)) {
								$bp_in_window = $hsp->end('hit') - ($i*300+1);
							} elsif($hsp->start('hit') >= ($i*300+1) && $hsp->end('hit') <= (($i+1)*300)) {
								$bp_in_window = $hsp->end('hit') - $hsp->start('hit');
							}else {
								$bp_in_window = 300;
							}
							$window_hash{$hit->name} = {
														$i => $bp_in_window
													};
						}
					# adding to parts that have already been initialized
					} else {
						# figure out where hte window starts and ends
						my $start_id = floor($hsp->start('hit')/300);
						my $end_id = ceil($hsp->end('hit')/300);
						my $i;
						# count the total # of bp within the current window covered by the current hit, add to the total
						for( $i = $start_id; $i <= $end_id; $i++ ) {
							# initialize new parts of the hash
							if( !(exists $window_hash{$hit->name}{$i}) ) {
								my $bp_in_window;
								if( $hsp->start('hit') >= ($i*300+1) && $hsp->end('hit') >= (($i+1)*300) ) {
									$bp_in_window = ($i+1)*300 - $hsp->start('hit');
								} elsif( $hsp->end('hit') <= (($i+1)*300 && $hsp->start('hit') <= ($i*300+1) && $hsp->end('hit') >= ($i*300+1)) ) {
									$bp_in_window = $hsp->end('hit') - ($i*300+1);
								} elsif($hsp->start('hit') >= ($i*300+1) && $hsp->end('hit') <= (($i+1)*300)) {
									$bp_in_window = $hsp->end('hit') - $hsp->start('hit');
								}else {
									$bp_in_window = 300;
								}
								$window_hash{$hit->name}{$i} = $bp_in_window;
							# adding to parts that have already been initalized
							} else {
								my $bp_in_window;
								if( $hsp->start('hit') >= ($i*300+1) && $hsp->end('hit') >= (($i+1)*300) ) {
									$bp_in_window = ($i+1)*300 - $hsp->start('hit');
								} elsif( $hsp->end('hit') <= (($i+1)*300) && $hsp->start('hit') <= ($i*300+1) && $hsp->end('hit') >= ($i*300+1)) {
									$bp_in_window = $hsp->end('hit') - ($i*300+1);
								} elsif($hsp->start('hit') >= ($i*300+1) && $hsp->end('hit') <= (($i+1)*300)) {
									$bp_in_window = $hsp->end('hit') - $hsp->start('hit');
								}else {
									$bp_in_window = 300;
								}
								$window_hash{$hit->name}{$i} += $bp_in_window;
							}
						}
					}
				}
		    }
		}

		return \%window_hash;
	}

	# Accepts a hash of averaged windows
	# converts the hash into a featurefile ready for gbrowse(basic xyplot config included at head of file)
	# returns the naem of the featurefile
	sub generate_featurefile {
		my $temp        = shift;
		my %window_hash = %{$temp};
		my $ref_fasta   = shift;
		my $title       = shift;
		my $query1      = shift;
		my $db          = shift;
		my $confid      = shift;

		my $ref_seq  = Bio::SeqIO->new(-file => $ref_fasta,
		                           		-format => 'Fasta'
		);

		my $max_score = 0;
		my $min_score = 0;
		my $count     = 1;
		my ($feature_fh, $feature_filename) = tempfile( DIR => $working_dir."output_files/" );
		select((select($feature_fh), $|=1)[0]);
		# print $feature_fh "##gff-version 3\n";
		my $seq;
		my $type = "region";
		my $last_seq;
		# step through the reference sequences to ensure that the xyplot covers all windows on the reference, not just ones that had hits(this is what would happen if we stepped through the hash instead)
		while( $seq = $ref_seq->next_seq ) {
			# split the current reference sequence into windows
			my $windows = floor($seq->length/300);
			my $ref_length = $seq->length;

			# printint out the results formatted into a featurefile
			my $reference = $seq->display_id;
			print $feature_fh "reference=$reference\n";
			# print $feature_fh "$reference\texpression\t$type\t1\t$ref_length\t0\t0\t0\tName=$title;ID=$count\n";
			my $i;
			for( $i = 0; $i<=$windows; $i++ ) {
				my $start = $i*300+1;
				my $end   = ($i+1)*300;
				if( $i == $windows) {
					$end = $seq->length;
				}
				my $score; 
				unless( $score = $window_hash{$reference}{$i} ) {
					$score = 0;
				}
				# calculate the minimum and maximum score shown in the results(needed for the configuration)
				$max_score = ($score > $max_score) ? $score : $max_score;
				$min_score = ($score < $min_score) ? $score : $min_score;

				print $feature_fh "expression\t$title\t$start..$end\tscore=$score\n";
				# print $feature_fh "$reference\texpression\t$type\t$start\t$end\t$score\t0\t0\tName=$title;Parent=$count\n";
			}
			print $feature_fh "\n";
			$count++;

			$last_seq = $seq->display_id;
		}

		print STDERR "LAST SEQUENCE ENCOUNTERED: $last_seq\n";

		close $feature_filename;

		my $config_string = "";	
		# add the basic configuration to the head of the featurefile
		$config_string .= "[$confid]\n";
		$config_string .= "feature = expression\n";
		$config_string .= "glyph = xyplot\n";
		$config_string .= "graph_type = boxes\n";
		$config_string .= "fgcolor = black\n";
		$config_string .= "bgcolor = red\n";
		$config_string .= "height = 100\n";
		$config_string .= "min_score = $min_score\n";
		$config_string .= "max_score = $max_score\n";
		$config_string .= "label = 1\n";
		$config_string .= "key = Expression Level\n\n";
		$config_string .= "\n";

		print STDERR "Config String: $config_string";

		my $file = $feature_filename;

		open my $in,  '<',  $file      or die "Can't read old file: $!";
		open my $out, '>', "$file.new" or die "Can't write new file: $!";

		print $out $config_string;

		while( <$in> ) {
		    print $out $_;
		}
		close $out;
		close $in;

		unlink($file);
		rename("$file.new", $file);


		return $file;
	}
}

1;