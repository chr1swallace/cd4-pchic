#!/usr/bin/perl

## Macd location (get from my github)
use lib '/home/oliver/GIT_REPOS/macd/lib';
use Macd;
use FindBin qw($Bin);
use File::Path qw(make_path);
use File::Find;
use strict;
use Data::Dumper;

## get these from environment variables

my $HOME = $ENV{HOME} || die "NO \$HOME variable set\n";
my $GRPATH = $ENV{GRPATH} || die "NO \$GRPATH variable set\n";
my $R_LIBS = $ENV{R_LIBS} || die "NO \$R_LIBS variable\n";
my $CHIGP_PATH = "$ENV{GRPATH}/cd4chic/gwas_paper/";

## site specific conf file for macd - see http://github.com/ollyburren/macd
## for more details
my $grid_cnf = "$GRPATH/macd/example/ini/example.cnf";

my $r_lib_dir = $R_LIBS;
	
my $DRIVER = Macd::GRIDDriverFactory->instantiate('SGE',inifile => $grid_cnf);


###############
#CONFIGURATION#
###############
## where to mail to once finished.
my $BASE_DIR = "$GRPATH/cd4chic/gwas_paper/DATA/";
my $DOTEST = 0; #IF set to true allows us to test the script by running only a few jobs
my $MANIFEST_FILE = "$BASE_DIR/support/gwas_manifest.csv";
my $ODIR = "$BASE_DIR/out/highscore_fragments/";
my $RSCRIPT="/home/oliver/bin/Rscript  $CHIGP_PATH/R/frag_level_GWAS.R";
## prior
my $PI_I=1e-4;
my $LOG_DIR = "$GRPATH/cd4chic/gwas_paper/DATA/log/highscore_frags/";

my $step = Macd::Step->new(
	logdir=> $LOG_DIR,
	driver=>$DRIVER,
	env_settings=>"R_LIBS=$r_lib_dir,GRPATH=$GRPATH"
);







my @GWAS;
## open and parse the manifest file
open(GWAS,$MANIFEST_FILE) || die "Cannot open $MANIFEST_FILE\n";
my @header=();
my @TRAIT;
while(<GWAS>){
	chomp;
	if(!@header){
		@header=split(",",$_);
		next;
	}else{
		my %tmp;
		my @vals=split(",",$_);
		for(my $x=0;$x<@header;$x++){
			$tmp{$header[$x]}=$vals[$x];
		}
		$tmp{n_samples}=$tmp{cases} + $tmp{controls};
		if($tmp{type} eq 'CC'){
			$tmp{prop_cases}=sprintf("%.3f",($tmp{cases}/$tmp{n_samples}));
		}else{
			$tmp{prop_cases}=1;
		}
		next if $tmp{label} eq 'RA';
		push @TRAIT,\%tmp;
	}
}

#print Dumper(\@TRAIT);
my $TEST==0;
my $addcount=0;
foreach my $t(@TRAIT){
	my $disease=$t->{label};
	next if -e "$ODIR/${disease}_17_05_16_gs_0.5.csv";
	next if $TEST > 1;
	my @param = "$RSCRIPT";
	push @param, "--disease=$disease";
	push @param, "--out_dir=$ODIR";
	my $cmd = join(" ",@param);
	#die "$cmd\n";
	my $job =  Macd::Step::Job->new(command=>$cmd);
	$step->add_job($job);
	$addcount++;
	$TEST++ if $DOTEST;
	##last;
	
}
print "Added $addcount jobs\n";


if($step->execute()){
	print "Step submitted successfully\n";
	## we can hold up prog execution as follows
	## in case we have a step below that requires the output
	$step->wait_on_complete();
	print "Step completed successfully\n";
}else{
	print "Step did not complete successfully\n";
}

