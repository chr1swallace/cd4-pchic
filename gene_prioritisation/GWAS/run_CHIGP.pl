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

my $MODE=shift;


my $HOME = $ENV{HOME} || die "NO \$HOME variable set\n";
my $GRPATH = $ENV{GRPATH} || die "NO \$GRPATH variable set\n";
my $R_LIBS = $ENV{R_LIBS} || die "NO \$R_LIBS variable\n";
my $CHIGP_PATH = "$ENV{GRPATH}/CHIGP/";


## site specific conf file for macd - see http://github.com/ollyburren/macd
## for more details
my $grid_cnf = "$GRPATH/macd/example/ini/example.cnf";

my $r_lib_dir = $R_LIBS;
	
my $DRIVER = Macd::GRIDDriverFactory->instantiate('SGE',inifile => $grid_cnf);



###############
#CONFIGURATION#
###############
## where to mail to once finished.
my $PREFIX = 'javierre_tnact_';
my $BASEDIR = "$GRPATH/cd4chic/gwas_paper/DATA/";
my $DOTEST = 0; #IF set to true allows us to test the script by running only a few jobs
my $CHIC_THRESH=5;
my $PMI_DIR = "$BASEDIR/out/pmi/";
my $INT_FILE="$BASEDIR/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75_TNACT.tab";
my $CSNP_FILE="$BASEDIR/support/cSNPs_w_ENSG.e75.bed";
my $FRAG_FILE="$BASEDIR/support/Digest_Human_HindIII.bed";
my $REG_FILE="$BASEDIR/support/0.1cM_regions.b37.bed";
my $RDATA_DIR="$BASEDIR/RDATA/";
my $CSNPS="$RDATA_DIR/${PREFIX}csnps.by.ld.RData";
my $INTS="$RDATA_DIR/${PREFIX}interactions.RData";
my $FRAGS="$RDATA_DIR/${PREFIX}frags.by.ld.RData";
my $TARGET_GENE_CSNPS_ONLY=1;
my $INCLUDE_INTERACTIONS_ONLY=1;
my $DECOMPOSE=1;
my $PPI_THRESHOLD=0.5;
my $BF_THRESHOLD=3;
my $SETS_YAML_FILE="$BASEDIR/support/tnact_tree.yaml";
my $SUPPORT_RSCRIPT="/home/oliver/bin/Rscript $CHIGP_PATH/R/generateResourceFiles.R";
my $RSCRIPT="/home/oliver/bin/Rscript $CHIGP_PATH/R/computeGeneScoreH.R";
my $OUTDIR = "$BASEDIR/out/tnact_hierarchical_geneScore/";
my $LOG_DIR = "$BASEDIR/log/tnact_hierachical_geneScore/";





my $step = Macd::Step->new(
	logdir=> $LOG_DIR,
	driver=>$DRIVER,
	env_settings=>"R_LIBS=$r_lib_dir,GRPATH=$GRPATH"
);

## if our support files don't exist create them first.

if(! -e $CSNPS){
	my $cmd = "$SUPPORT_RSCRIPT --prefix=$PREFIX --cSNP_file=$CSNP_FILE --interaction_file=$INT_FILE --pchic.thresh=$CHIC_THRESH --res_frag_file=$FRAG_FILE --region_bed_file=$REG_FILE out_dir=$RDATA_DIR";
	#die $cmd;
	`$cmd`
}


if(! -e $OUTDIR){
	print "Making path $OUTDIR\n";
	make_path($OUTDIR);
}

my $CONTROL;
my @GWAS;

## get all PMI 
find(sub{
	if(/\.pmi$/){                                 
		(my $stub = $_)=~s/\.pmi//;
		push @GWAS,$stub;
	}
},($PMI_DIR));



my $addcount;
foreach my $g(@GWAS){
	my $ofile = "$OUTDIR/${g}.pmi_prioritised.tab";
	print "$ofile\n";
	next if -e $ofile;
	my @param = "$RSCRIPT";
	push @param, "--pmi_file=$PMI_DIR$g.pmi";
	push @param, "--out_dir=$OUTDIR";
	push @param, "--int=$INTS";
	push @param, "--frags=$FRAGS";
	push @param, "--csnps=$CSNPS";
	push @param, "--target.gene.cSNPs.only=1";
	push @param, "--sets=$SETS_YAML_FILE";
	push @param, "--include.interactions.only=$INCLUDE_INTERACTIONS_ONLY";
	push @param, "--decompose=1";
	push @param, "--ppi.thresh=$PPI_THRESHOLD";
	push @param, "--BF.thresh=$BF_THRESHOLD";
	my $cmd = join(" ",@param);
	#die "$cmd";
	my $job =  Macd::Step::Job->new(command=>$cmd);
	$step->add_job($job);
	$addcount++;
	last if $DOTEST;
}
print "Added $addcount jobs\n";




if($step->execute()){
	print "Step submitted successfully\n";
	## we can hold up prog execution as follows
	## in case we have a step below that requires the output
	$step->wait_on_complete();
	print "Step completed successfully\n";
}else{
	print "Error submitting step\n";
}
