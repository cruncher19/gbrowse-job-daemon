#!/usr/bin/perl

# This script should be run as root. It will install the gbrowse_job_daemon, configure it as an upstart service, and create any folders necessary for the daemon to operate
# It will also optionall install the BlastSearch plugin for gbrowse_job_daemon

use Carp qw/croak/;
use File::Copy "cp";

umask 000;

if ($<) {
    croak "Error: exiting program, must be run as root\n";
}

my $temp;
my $daemon_dir = "/usr/share/gbrowse_job_daemon/";
my $daemon_exec = "gbrowse_job_daemon";
my $upstart_config = "gbrowse_job_daemon.conf";
my $upstart_location = "/etc/init/";
# my $daemon_dir = "/home/guest/installer_test/";

my $app_dir = "/usr/bin/";
print "Application installation directory[$app_dir]: ";
if( ($temp = <>) ne "\n" ) {
	$app_dir = add_slash($temp);
}
unless( -d $app_dir ) {
	croak "Error: directory: $app_dir does not exist!\n";
}

my $gbrowse_dir = "/etc/gbrowse2/";
print "GBrowse installation directory[$gbrowse_dir]: ";
if( ($temp = <>) ne "\n" ) {
	$gbrowse_dir = add_slash($temp);
}
unless( -d $gbrowse_dir ) {
	croak "Error: directory: $gbrowse_dir does not exist!\n";
}
unless( -d $gbrowse_dir."/plugins/" ) {
	croak "Error: there is no plugins directory in the gbrowse directory: $gbrowse_dir\n";
}


my $cgi_dir = "/usr/lib/cgi-bin/";
print "GBrowse cgi-bin directory[$cgi_dir]: ";
if( ($temp = <>) ne "\n" ) {
	$cgi_dir = add_slash($temp);	
}
unless( -d $cgi_dir ){
	croak "Error: directory: $cgi_dir does not exist!\n";
}
unless( -d $cgi_dir."/gb2/" ) {
	croak "Error: there is no gb2 folder in the cgi-bin directory: $cgi_dir\n";
}


my $www_dir = "/var/www/";
print "Where is your apache document root?[$www_dir]: ";
if( ($temp = <>) ne "\n" ) {
	$www_dir = add_slash($temp);
}
unless( -d $www_dir ) {
	croak "Error: directory: $www_dir does not exist!\n";
}

print "Directories:\n\tGBrowse install dir:\t$gbrowse_dir\n\tcgi_bin dir:\t$cgi_dir\n\tapache root:\t$www_dir\n";

mkdir $daemon_dir, 0777 or croak "Error setting up the daemon directory $daemon_dir/";
mkdir $daemon_dir."blast_targets/", 0777 or croak "Error setting up the daemon directory $daemon_dir/";
mkdir $daemon_dir."configs/", 0777 or croak "Error setting up the daemon directory $daemon_dir/";
mkdir $daemon_dir."cookies/", 0777 or croak "Error setting up the daemon directory $daemon_dir/";
mkdir $daemon_dir."output_files/", 0777 or croak "Error setting up the daemon directory $daemon_dir/";
mkdir $daemon_dir."sessions/", 0777 or croak "Error setting up the daemon directory $daemon_dir/";
mkdir $daemon_dir."temp_blastdb/", 0777 or croak "Error setting up the daemon directory $daemon_dir/";
mkdir $daemon_dir."temp_results/", 0777 or croak "Error setting up the daemon directory $daemon_dir/";

my @cgi_scripts = @{get_file_list("cgi/")};
map {
	print "Copying $_ to $cgi_dir"."gb2/\n"; 
	cp "cgi/".$_, $cgi_dir."gb2/" or croak "Unable to copy file $_: $!";
	unless( $_ =~ /.*\..*/){
		chmod 0755, $cgi_dir."gb2/$_";
	}
} @cgi_scripts;

my @www_files = @{get_file_list("www/")};
map {
	print "Copying $_ to $www_dir\n";
	cp "www/".$_, $www_dir or croak "Unable to copy file $_: $!";
} @www_files;

my $plugin_flag = 'y';
print "Install the BlastSearch GBrowse plugin(y/n)?[Y]: ";
if( ($temp = <>) ne "\n" ) {
	$plugin_flag = lc $temp;
}
unless( $plugin_flag eq 'n' || $plugin_flag eq 'y' ) {
	croak "Error: Invalid answer given\n";
}
if($plugin_flag eq 'y') {
	my @plugin_files = @{get_file_list("plugins/")};
	map {
		print "Copying $_ to $gbrowse_dir"."plugins/\n";
		cp "plugins/".$_, $gbrowse_dir."plugins/" or croak "Unable to copy file $_: $!";
	} @plugin_files;
}


print "Copying $daemon_exec to $app_dir\n";
cp $daemon_exec, $app_dir;
print "Copying $upstart_config to $upstart_location\n";
cp $upstart_config, $upstart_location;
print "registering daemon as an upstart service\n";
$temp = `initctl reload-configuration`;
print "Service Regestration complete, output: $temp\n";

print "\n\nInstallation completed successfully!\n";



sub add_slash {
	my $path = shift;
	chomp($path);

	if( $path =~ /\/$/ ) {
		return $path;
	} else {
		return $path."/";
	}
}

sub get_file_list {
	my $dir = shift;

	opendir DIR, $dir or die "cannot open dir $dir: $!";
	my @file= grep { !/^\./ } readdir DIR;
	closedir DIR;

	return \@file;
}

1;