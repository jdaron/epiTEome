#!/usr/bin/env perl

use strict;
use Getopt::Long ;
use warnings ;
use File::Basename ;
use Bio::SeqIO ;
use Bio::Seq::Quality;
use Bio::Tools::GFF;
use Bio::DB::Sam;
use Bio::DB::Bam::Alignment;
use Cwd;
use File::Which;
use Set::IntervalTree ;
use Statistics::Descriptive;
# use Devel::Size qw(total_size); # Need to be remove for public version
use Data::Dumper ;

my $prog = basename($0) ;
my $VERSION = '1.0';
my $lastmodif = '2017-02-01';
my $help;
my ($gffFile, $unmappedFastqFile, $targetList, $refFasta) ;
my $threads = 1 ;
my $bamFile ;
my $batchNb = 5000 ;
my $samOutputFile ;
my $lengthPairedEnd = '25,30,40' ;
my $windowSize = 10 ;
my $displayMet = "metaplot" ;
my $debug ; 
my $noclean ;
my $segemehl_Eval = 5 ;
my $segemehl_D = 1 ;

&GetOptions(
			"gff:s"                 => \$gffFile,
			"un:s"                  => \$unmappedFastqFile,
			"t:s"                   => \$targetList,
			"ref:s"                 => \$refFasta,
			"noclean"               => \$noclean,
			"p:i"                   => \$threads,
			"b:i"                   => \$batchNb,
			"chop:s"                => \$lengthPairedEnd,
			"w:i"                   => \$windowSize,
			"metO:s"                => \$displayMet,
			"debug"                 => \$debug,
			"E:i"                   => \$segemehl_Eval,
			"D:i"                   => \$segemehl_D,
			"h|help"                => \$help);


$help and &help;

&main();
sub main {
	my $self = {};
	bless $self;
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Start program!\n" ;

	$self->setOptions();
	$self->checkPrograms();
	$self->readTargetFile();
	$self->sF_TSDsize();
	$self->readGffFile();

	#---- Uncompress INPUT FASTQ if necessary
	my @unCompress ;
	if($self->{option}->{unmappedFastqFile}->{file}=~/(.*).bz2$/){
		@unCompress = ($self->{program}->{bzip2}, '-d', $self->{option}->{unmappedFastqFile}->{file}) ;
		$self->{option}->{unmappedFastqFile}->{file} = $1 ;
	}
	elsif($self->{option}->{unmappedFastqFile}->{file}=~/(.*).gz$/){
		@unCompress = ($self->{program}->{gzip}, '-d', $self->{option}->{unmappedFastqFile}->{file}) ;
		$self->{option}->{unmappedFastqFile}->{file} = $1 ;
	}
	# print STDERR join ( " ", @unCompress),"\n" ;
	system(join" ", @unCompress) ;
	#---- END

	#---- STEP 1: mapping of 'paired-end'
	$datestring = localtime();
	print STDERR "\nINFO $prog $datestring STEP 1: read ends mapping.\n" ;
	$self->{fileName}->{mappingOutput}->{bam} = $self->{option}->{unmappedFastqFile}->{prefix}.".step1.sort.bam" ;
	if(-e $self->{fileName}->{mappingOutput}->{bam}){
		print STDERR "INFO $prog: --> Skiping Step 1: $self->{fileName}->{mappingOutput}->{bam} allready exist.\n" ;
		$self->{fileName}->{splitFastq}->{inFastq} = $self->{option}->{unmappedFastqFile}->{file} ;
		$self->splitFastq("PaireEnd", "subset"); # splitFastq sub to find reads length of fastq file.
	}
	else{
		$self->{fileName}->{splitFastq}->{inFastq} = $self->{option}->{unmappedFastqFile}->{file} ;
		$self->splitFastq("PaireEnd");
		$self->{fileName}->{segemehl}->{-i}       = $self->{fileName}->{refFasta}->{file}.'.ctidx' ;
		$self->{fileName}->{segemehl}->{-j}       = $self->{fileName}->{refFasta}->{file}.'.gaidx' ;
		$self->{fileName}->{segemehl}->{-d}       = $self->{fileName}->{refFasta}->{file} ;
		$self->{fileName}->{segemehl}->{-q}       = $self->{fileName}->{splitFastq}->{outFastq} ;
		$self->{fileName}->{segemehl}->{-o}       = $self->{option}->{unmappedFastqFile}->{prefix}.".step1.sam" ;
		$self->{fileName}->{segemehl}->{-Eval}    = $segemehl_Eval ;
		$self->{fileName}->{mappingOutput}->{bam} = $self->{option}->{unmappedFastqFile}->{prefix}.".step1.sort.bam" ;
		$self->mappingMethylCSeq();

		if(-z $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.sort.bam'){
			print STDERR "DIE $prog: File $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.sort.bam' is empty... die!\n" and die ;
		}
		my @rm = ("rm", $self->{fileName}->{splitFastq}->{outFastq});
		system(join" ", @rm) ;
	}
	#---- END STEP 1

	#---- STEP 2: mapping split-reads
	$datestring = localtime();
	print STDERR "\nINFO $prog $datestring STEP 2: split-reads mapping.\n" ;
	for(1..$self->{batch}){
		$datestring = localtime();
		print STDERR "\nINFO $prog $datestring BATCH #$_.\n" ;
		$self->{fileName}->{mappingOutput}->{bam} = $self->{option}->{unmappedFastqFile}->{prefix}.".step1.sort.bam" ;

		foreach my $feat(sort{$self->{targetFeatTmp}->{$a}->{batch} <=> $self->{targetFeatTmp}->{$a}->{batch}} keys %{$self->{targetFeatTmp}}){
			if($self->{targetFeatTmp}->{$feat}->{batch} == $_){
				$self->{targetFeat}->{$feat} = $self->{targetFeatTmp}->{$feat} ;
			}
		}
		
		$self->step2() ;
		$self->{cleanNum} = 1 ;
		if(not defined $noclean){
			$self->clean($self->{cleanNum}) ;
		}
		delete($self->{targetFeat}) ;
	}
	#---- END STEP 2

	$self->makeWindows();
	$self->print() ;
	$self->printBam();
	$self->printRowDataMet();
	$self->calcMetPercentage();
    $self->calcMeanAndCi();
    $self->printMetaplot();
	$self->{cleanNum} = 2 ;
	if(not defined $noclean){
		$self->clean($self->{cleanNum}) ;
	}

	# compress starting fastq input
	my @bzip2 = ($self->{program}->{bzip2}, $self->{option}->{unmappedFastqFile}->{file} ) ;
	print STDERR join ( " ", @bzip2),"\n" ;
	system(join" ", @bzip2) ;
	$datestring = localtime();
	print STDERR "INFO $prog $datestring Finishied succesfully !\n" ;
	exit(0);
} # END OF main


###########
### sub ###
###########

sub setOptions {
	my $self = shift;
	$self->{dirName}->{currentDir} = cwd ;
	$self->{option}->{maxReadsLength} = 0 ;
	if (defined ($gffFile)){
		unless (-e $gffFile){ print STDERR "DIE $prog: Could not find file: " . $gffFile . "\n" and die ; }
		$self->{option}->{gff} = $gffFile ;
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -gff.\n" and &help ;
	}
	if (defined ($unmappedFastqFile)){
		unless (-e $unmappedFastqFile){ print STDERR "DIE $prog: Could not find file: " . $unmappedFastqFile . "\n" and die ; }
		$self->{dirName}->{unDir} = dirname($unmappedFastqFile) ;
		$self->{option}->{unmappedFastqFile}->{file} = $unmappedFastqFile ;
		 if($unmappedFastqFile =~/(.*).fastq$/){
		 	$self->{option}->{unmappedFastqFile}->{prefix} = $1 ;
		 }
		 elsif($unmappedFastqFile =~/(.*).fq$/){
		 	$self->{option}->{unmappedFastqFile}->{prefix} = $1 ;
		 }
		 elsif($unmappedFastqFile =~/(.*).fastq.bz2$/){
		 	$self->{option}->{unmappedFastqFile}->{prefix} = $1 ;
		 }
		 elsif($unmappedFastqFile =~/(.*).fq.bz2$/){
		 	$self->{option}->{unmappedFastqFile}->{prefix} = $1 ;
		 }
		 elsif($unmappedFastqFile =~/(.*).fastq.gz$/){
		 	$self->{option}->{unmappedFastqFile}->{prefix} = $1 ;
		 }
		 elsif($unmappedFastqFile =~/(.*).fq.gz$/){
		 	$self->{option}->{unmappedFastqFile}->{prefix} = $1 ;
		 }
		 else{
			print STDERR "ERROR $prog: File give in -un option should have extention fastq or fq.\n" and &help ;
		 }
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -un.\n" and &help ;
	}
	if (defined ($targetList)){
		unless (-e $targetList){ print STDERR "DIE $prog: Could not find file: " . $targetList . "\n" and die ; }
		$self->{option}->{targetList} = $targetList ;
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -t.\n" and &help ;
	}
	if (defined ($targetList)){
		unless (-e $targetList){ print STDERR "DIE $prog: Could not find file: " . $targetList . "\n" and die ; }
		$self->{option}->{targetList} = $targetList ;
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -t.\n" and &help ;
	}
	if (defined ($refFasta)){
		unless (-e $refFasta){ print STDERR "DIE $prog: Could not find file: " . $refFasta . "\n" and die ; }
		 if($refFasta =~/(.*).fasta/ or $refFasta =~/(.*).fna/){
			unless (-e $refFasta.".ctidx"){ print STDERR "DIE $prog: Could not find segemehl index file .ctidx \n" and &help ; }
			unless (-e $refFasta.".gaidx"){ print STDERR "DIE $prog: Could not find segemehl index file .gaidx \n" and &help ; }
			$self->{fileName}->{refFasta}->{file} = $refFasta ;
	 }
		 else{
			print STDERR "ERROR $prog: File give in -fasta option should have extention fasta or fna.\n" and &help ;
		 }
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -fasta.\n" and &help ;
	}
	if(defined $lengthPairedEnd){
		if($lengthPairedEnd =~/,/){
			my @array = split(",", $lengthPairedEnd) ;
			foreach my $nb (sort  @array){
				push(@{$self->{option}->{lengthPairedEnd}}, $nb) ;
			}
		}
		else{
			push(@{$self->{option}->{lengthPairedEnd}}, $lengthPairedEnd) ;
		}
	}
}
sub checkPrograms{
	my $self = shift;
	$self->{program}->{bzip2} = which 'bzip2';
	if(not defined $self->{program}->{bzip2}){ print STDERR "DIE $prog: Missing program: bzip2.\n" and die ; }
	$self->{program}->{gzip} = which 'gzip';
	if(not defined $self->{program}->{gzip}){ print STDERR "DIE $prog: Missing program: gzip.\n" and die ; }
	$self->{program}->{samtools} = which 'samtools';
	if(not defined $self->{program}->{samtools}){ print STDERR "DIE $prog: Missing program: samtools.\n" and die ; }
	$self->{program}->{maskFastaFromBed} = which 'maskFastaFromBed';
	if(not defined $self->{program}->{maskFastaFromBed}){ print STDERR "DIE $prog: Missing program: maskFastaFromBed.\n" and die ; }
	$self->{program}->{fastqutils} = which 'fastqutils';
	if(not defined $self->{program}->{fastqutils}){ print STDERR "DIE $prog: Missing program: fastqutils.\n" and die ; }
	$self->{program}->{bamutils} = which 'bamutils';
	if(not defined $self->{program}->{bamutils}){ print STDERR "DIE $prog: Missing program: bamutils.\n" and die ; }
	$self->{program}->{segemehl} = which 'segemehl.x';
	if(not defined $self->{program}->{segemehl}){ print STDERR "DIE $prog: Missing program: segemehl.x.\n" and die ; }
	$self->{program}->{rm} = which 'rm';
	if(not defined $self->{program}->{rm}){ print STDERR "DIE $prog: Missing program: rm.\n" and die ; }
	$self->{program}->{awk} = which 'awk';
	if(not defined $self->{program}->{awk}){ print STDERR "DIE $prog: Missing program: awk.\n" and die ; }
	# print STDERR "INFO $prog: All programs have been found succesfully.\n" ;
}

sub step2{
	my $self = shift;
	
	######################
	#### If only for development 

	if(-e $self->{option}->{unmappedFastqFile}->{prefix}.".step2.sort.bam"){ # Only for developpement
		print STDERR "INFO $prog: --> Skiping Step 2: $self->{option}->{unmappedFastqFile}->{prefix}.step2.sort.bam allready exist.\n" ;
		$self->{fileName}->{mappingOutput}->{bam} = $self->{option}->{unmappedFastqFile}->{prefix}.".step2.sort.bam" ;

		#check if $self->{fileName}->{mappingOutput}->{bam} isn't either empty or troncated
		if(-z $self->{fileName}->{mappingOutput}->{bam}){
			print STDERR "INFO $prog: $self->{fileName}->{mappingOutput}->{bam} is empty.\n" ;
			return (0) ;
		}
		my $out=`samtools view -H $self->{fileName}->{mappingOutput}->{bam} | egrep -c "\@SQ"`;
	 	chomp $out ;
		if($out == 0){
			print STDERR "INFO $prog: $self->{fileName}->{mappingOutput}->{bam} header troncated.\n" ;
			return (0) ;
		}

		$self->selectSplitReads_v2(); # select split-reads mapped on the targeted TEs.
		$self->splitReadsGrab_v2();
		$self->groupOverlapingmNeoReads();
		$self->neoInsertionFinder();
	}
	else{
		$self->selectPairedEndReads();
		$self->pairedEndGrab();
		$self->checkPairedEndPosition_v3(); # regarding other copies in the same family

		my @fastqutils = ($self->{program}->{fastqutils} ,'filter -whitelist', $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.list', $self->{option}->{unmappedFastqFile}->{file}, ' > ',  $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.fastq') ;
		print STDERR join ( " ", @fastqutils),"\n" ;
		system(join" ", @fastqutils) ;
		if(-z $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.fastq'){
			print STDERR "INFO $prog: $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.fastq' is empty.\n" ;
			return (0) ;
		}

		# step 2 mapping 
		$self->{fileName}->{mappingOutput}->{bam} = $self->{option}->{unmappedFastqFile}->{prefix}.".step2.sort.bam" ;
		$self->{fileName}->{splitFastq}->{inFastq} = $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.fastq' ;
		$self->splitFastq("eachNucl");
		$self->{fileName}->{segemehl}->{-i}       = $self->{fileName}->{refFasta}->{file}.'.ctidx' ;
		$self->{fileName}->{segemehl}->{-j}       = $self->{fileName}->{refFasta}->{file}.'.gaidx' ;
		$self->{fileName}->{segemehl}->{-d}       = $self->{fileName}->{refFasta}->{file} ;
		$self->{fileName}->{segemehl}->{-q}       = $self->{fileName}->{splitFastq}->{outFastq} ;
		$self->{fileName}->{segemehl}->{-o}       = $self->{option}->{unmappedFastqFile}->{prefix}.".step2.sam" ;
		$self->{fileName}->{mappingOutput}->{bam} = $self->{option}->{unmappedFastqFile}->{prefix}.".step2.sort.bam" ;
		$self->{fileName}->{segemehl}->{-Eval}    = $segemehl_Eval ;
		$self->mappingMethylCSeq();

		#check if $self->{fileName}->{mappingOutput}->{bam} isn't either empty or troncated
		if(-z $self->{option}->{unmappedFastqFile}->{prefix}.".step2.sort.bam"){
			print STDERR "INFO $prog: $self->{fileName}->{mappingOutput}->{bam} is empty.\n" ;
			return (0) ;
		}
		my $out=`samtools view -H $self->{fileName}->{mappingOutput}->{bam} | egrep -c "\@SQ"`;
	 	chomp $out ;
		if($out == 0){
			print STDERR "INFO $prog: $self->{fileName}->{mappingOutput}->{bam} header troncated.\n" ;
			return (0) ;
		}

		$self->selectSplitReads_v2(); # select split-reads mapped on the targeted TEs.
		$self->splitReadsGrab_v2();
		$self->groupOverlapingmNeoReads();
		$self->neoInsertionFinder();
	}
}
sub clean {
	my $self = shift ;
	my $num = shift ;
	my $file ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".pe.fastq"} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".matePairedEnd.bam"} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".matePairedEnd.bam.bai"} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".fishReads.list"} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".step1.list"} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.'.step1.fastq'} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".step1.sr.fastq"} = 1 ;
 	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".step2.sort.bam"} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".step2.sort.bam.bai"} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".mateSplitReads.bam"} = 1 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".mateSplitReads.bam.bai"} = 1 ;

  	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".step1.sort.bam"} = 2 ;
	$file->{$self->{option}->{unmappedFastqFile}->{prefix}.".step1.sort.bam.bai"} = 2 ;
	$file->{"log"} = 2 ;

	# $file->{$self->{fileName}->{refSeqNewPrefix}."_refMasked*"} = 1 ;


	my @rm = $self->{program}->{rm} ;
	foreach my $fileName(sort keys %{$file}){
		if($file->{$fileName} <= $num){
			if(-e $fileName){
				push(@rm, $fileName) ;
			}
		}
	}
	if(scalar @rm>1){
		print STDERR "INFO $prog: Cleaning file: " ;
		print STDERR join ( " ", @rm),"\n" ;
		system(join" ", @rm) ;
	}
}
sub readTargetFile {
	my $self = shift;
	open("T",$self->{option}->{targetList}) or print STDERR "Die: Could not open file: " . $self->{option}->{targetList} . "\n" and exit ;
	while(<T>){
		chomp $_ ;
		$self->{targetList}->{$_} =  1 ;
	}
}
sub readGffFile {
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: readGffFile!\n" ;

	my $self = shift;
	my $gffio = Bio::Tools::GFF->new(-file => $self->{option}->{gff}, -gff_version => 3);
    my $feature;
    my $checkConcodanceGffTarget = "FALSE" ;
    my $nbTargetTEsPerFam ;

    # loop over the input stream
    while($feature = $gffio->next_feature()) {
    	# if ( $feature->strand eq "+" ){$feature->strand = 1 ; }
    	# if ( $feature->strand eq "-" ){$feature->strand = -1 ;}
    	if($feature->has_tag('ID') and $feature->has_tag('fam')){
    		my @id = $feature->get_tag_values('ID') ;
    		my @fam = $feature->get_tag_values('fam') ;
    		my $family = $fam[0] ;
    		if($fam[0]=~/ATHILA.*/){
    			$family = "ATHILA" ;
    			$feature->remove_tag('fam');
    			$feature->add_tag_value('fam',"ATHILA");
    		}
   			$self->{family}->{$family}->{$id[0]} = $feature ;

    	    # Index region into interval tree
	        my $obj = {
	            name   => $id[0],
	            start  => $feature->start-50,
	            stop   => $feature->end+50,
	            strand => $feature->strand
	        };
	        $self->{bedTree}->{$family}->{$feature->seq_id} = Set::IntervalTree->new unless ( defined $self->{bedTree}->{$family}->{$feature->seq_id} );
	        $self->{bedTree}->{$family}->{$feature->seq_id}->insert( $obj, $feature->start-50, $feature->end+50);

    		if(defined $self->{targetList}->{$id[0]}){
    			if($feature->end-$feature->start+1>$self->{option}->{maxReadsLength}){
	    			$checkConcodanceGffTarget = "TRUE" ;
	    			$self->{targetFeatTmp}->{$id[0]}->{$feature->primary_tag} = $feature ;
	    			$nbTargetTEsPerFam->{$family} ++ ;
    			}
    		}
    	}
    	elsif($feature->has_tag('Parent') and ($feature->primary_tag eq "LTR5" ) or ($feature->primary_tag eq "LTR3") ){
    		my @parent = $feature->get_tag_values('Parent') ;
    		if(defined $self->{targetList}->{$parent[0]}){
    			if($feature->end-$feature->start+1>$self->{option}->{maxReadsLength}){
	    			$self->{targetFeatTmp}->{$parent[0]}->{$feature->primary_tag} = $feature ;
				}
    		}
    	}
    }
    $gffio->close();

	# deal with batch 
	$self->{batch} = 1 ;
	my $nbCopiesWithinFam = 0 ;
	foreach my $family (sort keys %{$nbTargetTEsPerFam}){
		$nbCopiesWithinFam += $nbTargetTEsPerFam->{$family} ;
		if($nbCopiesWithinFam >= $batchNb){
			$nbCopiesWithinFam = 0 ;
			$self->{batch}++ ;
		}
		foreach my $teid (sort keys %{$self->{family}->{$family}}){
			if (defined $self->{targetFeatTmp}->{$teid}){
	    		$self->{targetFeatTmp}->{$teid}->{batch} = $self->{batch} ;
			}
		}	
	}
	if($checkConcodanceGffTarget eq "FALSE"){
		print STDERR "DIE $prog: Features ID in target file (-t) are not present in gff file. Die.\n" and die ;
	}
	delete($self->{targetList}) ;
	# $self->{targetList} = () ;
	# delete($nbTargetTEsPerFam) ;
	return(1) ;
}
sub splitFastq {
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: splitFastq!\n" ;

	my $self = shift ;
	my $method = shift ;
	my $subset = shift ;
	my $k = 1 ; 
	my $step = 1 ;
	my $min = $self->{option}->{lengthPairedEnd}[0] ;
	my $prefix ;

	unless (-e $self->{fileName}->{splitFastq}->{inFastq}){ print STDERR "DIE $prog: Could not find file: " . $self->{fileName}->{splitFastq}->{inFastq} . "\n" and die ; }
	if($self->{fileName}->{splitFastq}->{inFastq} =~/(.*).fastq$/){
		$prefix = $1 ;
	}
	elsif($self->{fileName}->{splitFastq}->{inFastq} =~/(.*).fq$/){
		$prefix = $1 ;
	}
	elsif($self->{fileName}->{splitFastq}->{inFastq} =~/(.*).fastq.bz2$/){
		$prefix = $1 ;
	}
	elsif($self->{fileName}->{splitFastq}->{inFastq} =~/(.*).fq.bz2$/){
		$prefix = $1 ;
	}
	else{
		die and print STDERR "DIE $prog: Only .fastq or .fq exention accepted in $self->{fileName}->{splitFastq}->{inFastq}\n" ;
	}
	my $r1 ;
	my $fileName ;
	if($method eq "eachNucl"){
		$fileName = $prefix.".sr.fastq" ; 
		$self->{fileName}->{splitFastq}->{outFastq} = $prefix.".sr.fastq" ;
		$self->{fileName}->{splitFastq}->{prefix} = $prefix ;
	}
	elsif($method eq "PaireEnd"){
		$fileName = $prefix.".pe.fastq" ; 
		$self->{fileName}->{splitFastq}->{outFastq} = $prefix.".pe.fastq" ;
		$self->{fileName}->{splitFastq}->{prefix} = $prefix ;
	}

	open($r1, '>', $fileName) or die("Could not open file.\n");
	my $in = Bio::SeqIO->new(-format  => 'fastq',
                             -file    => $self->{fileName}->{splitFastq}->{inFastq});

	while (my $data = $in->next_dataset) {
		if(defined $subset){
			$k++ ;
			if($k>1000){ last ;}
		}
  		my $seq =$data->{"-seq"} ;
  		my $length = length($seq) ;
  		if($length > $self->{option}->{maxReadsLength}){$self->{option}->{maxReadsLength} = $length ;}
		if($method eq "eachNucl"){

			if(defined $self->{storeChop}->{$data->{"-id"}}){ # if you delete checkPairedEndPosition_v3, you should delete also this part.
				$min = $self->{storeChop}->{$data->{"-id"}} ;
			}

			for (my $i = $min ; $i<$length-$min ; $i = $i+$step){
				my $id = $data->{'-id'}."_".$i."_".$length ;
	            print $r1 "@".$id."/1\n" ;
	            print $r1 (substr $data->{"-seq"}, 1, $i) ;
	            print $r1 "\n" ;
	            print $r1 "+\n" ;
	            print $r1 (substr $data->{"-raw_quality"}, 1, $i) ;
	            print $r1 "\n" ;
	            # h2
	            print $r1 "@".$id."/2\n" ;
	            print $r1 (substr $data->{"-seq"}, $i+1, $length);
	            print $r1 "\n" ;
	            print $r1 "+\n" ;
	            print $r1 (substr $data->{"-raw_quality"}, $i+1, $length) ;
	            print $r1 "\n" ;
			}
		}
		elsif($method eq "PaireEnd"){
			foreach my $lengthPairedEnd (@{$self->{option}->{lengthPairedEnd}}){
				my $id = $data->{'-id'} ;
				# h1
				print $r1 "@".$id."_".$lengthPairedEnd."/1\n" ;
				print $r1 (substr $data->{"-seq"}, 1, $lengthPairedEnd) ;
				print $r1 "\n" ;
				print $r1 "+\n" ;
				print $r1 (substr $data->{"-raw_quality"}, 1, $lengthPairedEnd) ;
				print $r1 "\n" ;
				# h2
				print $r1 "@".$id."_".$lengthPairedEnd."/2\n" ;
				print $r1 (substr $data->{"-seq"}, ($length-$lengthPairedEnd), $length);
				print $r1 "\n" ;
				print $r1 "+\n" ;
				print $r1 (substr $data->{"-raw_quality"}, ($length-$lengthPairedEnd), $length) ;
				print $r1 "\n" ;
			}
		}
	}
	close $r1 ;
	if($method eq "eachNucl"){ delete($self->{storeChop}) ; }
}
sub mappingMethylCSeq {
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: Segemehl mapping !\n" ;

	my $self = shift ;
	# check presence of file
	unless(-e $self->{fileName}->{segemehl}->{-i}){
		print STDERR "DIE $prog: Could not file Segemehl input file: $self->{fileName}->{segemehl}->{-i}.\n" and &help ;
	}
	unless(-e $self->{fileName}->{segemehl}->{-d}){
		print STDERR "DIE $prog: Could not file Segemehl input file: $self->{fileName}->{segemehl}->{-d}.\n" and &help ;
	}
	unless(-e $self->{fileName}->{segemehl}->{-q}){
		print STDERR "DIE $prog: Could not file Segemehl input file: $self->{fileName}->{segemehl}->{-q}.\n" and &help ;
	}

	my @map = ($self->{program}->{segemehl}, '--silent -t ',$threads, '-D', $segemehl_D, '-E ', $self->{fileName}->{segemehl}->{-Eval}, ' -i ', $self->{fileName}->{segemehl}->{-i}, '-j', $self->{fileName}->{segemehl}->{-j} , '-d', $self->{fileName}->{segemehl}->{-d}, '-q', $self->{fileName}->{segemehl}->{-q}, '-o', $self->{fileName}->{segemehl}->{-o}, ' -F 1 2>log') ;
	print STDERR join ( " ", @map),"\n" ;
	system(join" ", @map) ;
	# my @convertInBam = ($self->{program}->{awk}, '\'{if($6 ~/^[0-9]+M$/){print $0;} else if($1 ~/^@.*/){print $0;}}\'', $self->{fileName}->{segemehl}->{-o}, '|', $self->{program}->{samtools}, 'view -b - |', $self->{program}->{samtools},'sort - -o', $self->{fileName}->{segemehl}->{bam} ) ;
	my @convertInBam = ($self->{program}->{samtools}, 'view -b', $self->{fileName}->{segemehl}->{-o}, '|', $self->{program}->{samtools},'sort - -o', $self->{fileName}->{mappingOutput}->{bam} ) ;
	print STDERR join ( " ", @convertInBam),"\n" ;
	system(join" ", @convertInBam) ;
	my @rm = ($self->{program}->{rm}, $self->{fileName}->{segemehl}->{-o}) ;
	system(join" ", @rm) ;
	my @index = ($self->{program}->{samtools}, 'index ', $self->{fileName}->{mappingOutput}->{bam} ) ;
	print STDERR join ( " ", @index),"\n" ;
	system(join" ", @index) ;
}
sub selectPairedEndReads{
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: selectPairedEndReads !\n" ;
	my $self = shift ;
	$self->{sam} = Bio::DB::Sam->new(-bam => $self->{fileName}->{mappingOutput}->{bam} );
	my $list ;
	my $nbReadsStore = 0 ;
	open($list, '>', $self->{option}->{unmappedFastqFile}->{prefix}.".fishReads.list") or die("DIE: Could not open file.\n");
	foreach my $teid (sort keys %{$self->{targetFeat}}){
		my @fam ;
		if($self->{targetFeat}->{$teid}->{'te'}->has_tag('fam')){
			@fam = $self->{targetFeat}->{$teid}->{'te'}->get_tag_values('fam') ;
		}
		else{
			print STDERR "WARNING $prog: Not tag fam for feature $teid in Gff3.\n" ;
			next ;
		}
		my $teObj = $self->{targetFeat}->{$teid}->{te} ;
		my $inside = $teObj->start+$self->{option}->{maxReadsLength} ;
		if($inside>$teObj->end){ $inside = $teObj->end ; }

		##### for 5prime of the TEs.
		my $position ;
		if($teObj->start-20 < 1){
			$position = $teObj->seq_id.":1-".$inside ; # select reads mapped on each TEs of interested. Position = start-20, start+200
		}
		else {
			$position =  $teObj->seq_id.":".($teObj->start-20)."-".$inside ; # select reads mapped on each TEs of interested. Position = start-20, start+200

		}

		$self->{sam}->fetch($position,
		    sub {
				my $a = shift;
				if($a->strand eq "1"){
					my $readName = $a->display_name ;
					if($readName  =~/^(.*)_([0-9]+)\/2$/){
						my $mate = &reverseMate($readName) ;
						if(not defined $self->{targetFeat}->{$teid}->{reads}->{$mate}){
							print $list $mate."\n" ;
							$self->{targetFeat}->{$teid}->{reads}->{$mate} = 1 ;
						}
					}
				}
				elsif($a->strand eq "-1"){
					my $readName = $a->display_name ;
					if($readName  =~/^(.*)_([0-9]+)\/1$/){
						my $mate = &reverseMate($readName) ;
						if(not defined $self->{targetFeat}->{$teid}->{reads}->{$mate}){
							print $list $mate."\n" ;
							$self->{targetFeat}->{$teid}->{reads}->{$mate} = 1 ;
						}
					}
				}
			}
		);

		#### for 3prime of the TEs.
		$inside = $teObj->end-$self->{option}->{maxReadsLength} ;
		if($inside < $teObj->start){ $inside = $teObj->start ; }
		$position = $teObj->seq_id.":".$inside."-".($teObj->end+20) ; # select reads mapped on each TEs of interested. Position = end-200, end+50

		$self->{sam}->fetch($position,
		    sub {
				my $a = shift;
				if($a->strand eq "1"){
					my $readName = $a->display_name ;
					if($readName  =~/^(.*)_([0-9]+)\/1$/){
						my $mate = &reverseMate($readName) ;
						if(not defined $self->{targetFeat}->{$teid}->{reads}->{$mate}){
							print $list $mate."\n" ;
							$self->{targetFeat}->{$teid}->{reads}->{$mate} = 1 ;
						}
					}
				}
				elsif($a->strand eq "-1"){
					my $readName = $a->display_name ;
					if($readName  =~/^(.*)_([0-9]+)\/2$/){
						my $mate = &reverseMate($readName) ;
						if(not defined $self->{targetFeat}->{$teid}->{reads}->{$mate}){
							print $list $mate."\n" ;
							$self->{targetFeat}->{$teid}->{reads}->{$mate} = 1 ;
						}
					}
				}
			}
		);
		if(defined $self->{targetFeat}->{$teid}->{reads}){
			$nbReadsStore += scalar(keys %{$self->{targetFeat}->{$teid}->{reads}}) ;
		}
		print STDERR $teid."\t".$nbReadsStore."\r" ;
	}
	print STDERR "\n" ;
}
sub reverseMate{
	my $readName = shift ;
	my $mate ;
	if($readName  =~/^(.*)_([0-9]+)\/1$/){ # for h1
		$mate = $1."_".$2."/2" ;
		}
	elsif($readName =~/^(.*)_([0-9]+)\/2$/){ # for h2
		$mate = $1."_".$2."/1" ;
	}
	return($mate) ;
}
sub checkPairedEndPosition_v3{
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: checkPairedEndPosition !\n" ;

	my $self = shift ;
	my $list ;
	open($list, '>', $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.list') or die("DIE: Could not open file step1.list.\n");
	foreach my $teid (sort keys %{$self->{targetFeat}}){
		my @fam = $self->{targetFeat}->{$teid}->{'te'}->get_tag_values('fam') ;

		my $alreadyPassed ;
		READ:foreach my $readName (sort keys %{$self->{targetFeat}->{$teid}->{reads}}){
			my ($read, $chopSize) = $readName =~/^(.*)_([0-9]+)\/[1-2]$/ ;
			if(defined $alreadyPassed->{$read}){ next READ;}
			if(defined $self->{storeReads}->{allReads}->{$readName}){ # if mate is mapped.
				foreach my $num (sort keys %{$self->{storeReads}->{allReads}->{$readName}}){
					if(defined $alreadyPassed->{$read}){ next READ;}
					my $readChr = $self->{storeReads}->{allReads}->{$readName}->{$num}->seq_id ;
					my $readStart = $self->{storeReads}->{allReads}->{$readName}->{$num}->start ;
					my $readStop = $self->{storeReads}->{allReads}->{$readName}->{$num}->end ;

			        unless (defined $self->{bedTree}->{$fam[0]}->{$readChr}){ # read Chr different than copies Chr
						print $list $read."\n" ;
						if(defined $self->{storeChop}->{$read}){
							if($self->{storeChop}->{$read} < $chopSize){
								$self->{storeChop}->{$read} = $chopSize ;
							}
						}
						else{
							$self->{storeChop}->{$read} = $chopSize ;
						}
						if($chopSize == @{$self->{option}->{lengthPairedEnd}}[scalar(@{$self->{option}->{lengthPairedEnd}})-1]){
							$alreadyPassed->{$read} = 1 ;
						}
						next ; # because $bedTree->{$readChr} is not defined in this case
			        };
					
			        my $overlapping = $self->{bedTree}->{$fam[0]}->{$readChr}->fetch( $readStart, $readStop );
			        unless(defined $overlapping and @$overlapping){ # read Chr is not overlapping copies
						print $list $read."\n" ;
						if(defined $self->{storeChop}->{$read}){
							if($self->{storeChop}->{$read} < $chopSize){
								$self->{storeChop}->{$read} = $chopSize ;
							}
						}
						else{
							$self->{storeChop}->{$read} = $chopSize ;
						}
						if($chopSize == @{$self->{option}->{lengthPairedEnd}}[scalar(@{$self->{option}->{lengthPairedEnd}})-1]){
							$alreadyPassed->{$read} = 1 ;
						}
			        }
				}
			}
			else{ # the other mate isn't mapped => keep it.
				print $list $read."\n" ;
				$alreadyPassed->{$read} = 1 ;
			}
		}
		delete($self->{targetFeat}->{$teid}->{reads}) ;
	}
	$self->{storeReads} = () ;
}
sub checkPairedEndPosition_v2{
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: checkPairedEndPosition !\n" ;

	my $self = shift ;
	my $list ;
	open($list, '>', $self->{option}->{unmappedFastqFile}->{prefix}.'.step1.list') or die("DIE: Could not open file step1.list.\n");
	foreach my $teid (sort keys %{$self->{targetFeat}}){
		my @fam = $self->{targetFeat}->{$teid}->{'te'}->get_tag_values('fam') ;

		my $alreadyPassed ;
		READ:foreach my $readName (sort keys %{$self->{targetFeat}->{$teid}->{reads}}){
			my ($read) = $readName =~/^(.*)_[0-9]+\/[1-2]$/ ;
			if(defined $alreadyPassed->{$read}){ next READ;}
			if(defined $self->{storeReads}->{allReads}->{$readName}){ # if mate is mapped.
				foreach my $num (sort keys %{$self->{storeReads}->{allReads}->{$readName}}){
					if(defined $alreadyPassed->{$read}){ next READ;}
					my $readChr = $self->{storeReads}->{allReads}->{$readName}->{$num}->seq_id ;
					my $readStart = $self->{storeReads}->{allReads}->{$readName}->{$num}->start ;
					my $readStop = $self->{storeReads}->{allReads}->{$readName}->{$num}->end ;

			        unless (defined $self->{bedTree}->{$fam[0]}->{$readChr}){ # read Chr different than copies Chr
						print $list $read."\n" ;
						$alreadyPassed->{$read} = 1 ;
						next ; # because $bedTree->{$readChr} is not defined in this case
			        };
					
			        my $overlapping = $self->{bedTree}->{$fam[0]}->{$readChr}->fetch( $readStart, $readStop );
			        unless(defined $overlapping and @$overlapping){ # read Chr is not overlapping copies
						print $list $read."\n" ;
						$alreadyPassed->{$read} = 1 ;
			        }
				}
			}
			else{ # the other mate isn't mapped => keep it.
				print $list $read."\n" ;
				$alreadyPassed->{$read} = 1 ;
			}
		}
		delete($self->{targetFeat}->{$teid}->{reads}) ;
	}
	$self->{storeReads} = () ;
}
sub pairedEndGrab {
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: pairedEndGrab !\n" ;

	my $self = shift;
	my $sam ;
	my @features;

	my $fileName = $self->{option}->{unmappedFastqFile}->{prefix}.".matePairedEnd.bam" ;
	my @bamutils = ($self->{program}->{bamutils}, 'filter ', $self->{fileName}->{mappingOutput}->{bam} , $fileName, '-whitelist ', $self->{option}->{unmappedFastqFile}->{prefix}.'.fishReads.list') ;
	print STDERR join ( " ", @bamutils),"\n" ;
	system(join" ", @bamutils) ;
	my @index = ($self->{program}->{samtools},'index', $fileName) ;
	print STDERR join ( " ", @index),"\n" ;
	system(join" ", @index) ;
	$sam = Bio::DB::Sam->new(-bam => $fileName);
	@features = $sam->features("match") ; # output all matching reads.

	for my $a (@features){
	 	if(not defined $self->{storeReads}->{allReads}->{$a->display_name}){
			$self->{storeReads}->{allReads}->{$a->display_name}->{1} = $a ;
 		}
 		elsif(defined $self->{storeReads}->{allReads}->{$a->display_name}){
 			my $num = scalar (keys %{$self->{storeReads}->{allReads}->{$a->display_name}})+1 ;
			$self->{storeReads}->{allReads}->{$a->display_name}->{$num} = $a ;
 		}
	}
}
sub selectSplitReads_v2{
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: selectSplitReads !\n" ;

	my $self = shift;
	$self->{sam} = Bio::DB::Sam->new(-bam   => $self->{fileName}->{mappingOutput}->{bam} , 
		                             -fasta => $self->{fileName}->{refFasta}->{file});

	# select the bam file header to later print bam.
	if(not defined $self->{printBam}->{header}){
		my $bam = Bio::DB::Bam->open($self->{fileName}->{mappingOutput}->{bam}) ;
		$self->{printBam}->{header} = $bam->header() ;
	}

	my $list ;
	my $passed ;
	open($list, '>', $self->{option}->{unmappedFastqFile}->{prefix}.".fishReads.list") or die("DIE: Could not open file.\n");
	foreach my $teid (sort keys %{$self->{targetFeat}}){
		print STDERR $teid."\r" ;
		my $teObj = $self->{targetFeat}->{$teid}->{te} ;
		my @fam = $self->{targetFeat}->{$teid}->{'te'}->get_tag_values('fam') ;
		##### for 5prime of the TEs.
		my $position = $teObj->seq_id.":".($teObj->start)."-".($teObj->start+$self->{option}->{maxReadsLength}) ;
		$self->{sam}->fetch($position,
	      sub {
				my $a = shift;
	       		my ($readId, $halfRead, $size, $splitPosition, $break, $fullReadLength) = $self->spliceReadsId($a->display_name) ;
            	if($halfRead == 1){ # for h1
	          		if($a->strand == "-1"){ # if h1 and strand = + => useless

						# selected mNeoReadId
          				for(my $i = $break ; $i<$break+10 ; $i++){
			               	my $mNeoReadId = $readId."_".$i."_".$fullReadLength."/2" ;
          					if(not defined $passed->{$mNeoReadId}){
				               	print $list $mNeoReadId."\n" ;
          					}
          					$passed->{$mNeoReadId} = 1 ;
          				}

		               	# load $self->{nbTimeMapped} for multi-mapped reads
		               	if(defined $self->{nbTimeMapped}->{mTE}->{$a->display_name}){
		               		$self->{nbTimeMapped}->{mTE}->{$a->display_name} += 1 ;
		               	}
		               	else{
		               		$self->{nbTimeMapped}->{mTE}->{$a->display_name} = 1 ;
		               	}

		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{subPart}   = 1 ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{featureID} = $teid ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{fam}       = $fam[0] ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{obj}       = $a ;

		               	# exact insertion
		               	if($teObj->start == $a->start){ # read mapped exactly at the start
		               		if($teObj->strand == 1){
				               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mappedInStart} = 1 ;
		               		}
		               		elsif($teObj->strand == 1){
				               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mappedInStart} = 1 ;
		               		}
			               	$self->{mappedOnStartEnd}->{$readId}->{$splitPosition} = 1 ;
		               	}

		               	# determine mtExtremity
		               	if($teObj->strand == 1){
			            	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mtExtremity} = 5 ;
		               	}
		               	elsif($teObj->strand == -1){
			            	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mtExtremity} = 3 ;
		               	}
					}
	       		}
	       		elsif($halfRead == 2){ # for h2
	        		if($a->strand == "+1"){ # if h2 and strand = - => useless

						# selected mNeoReadId
          				for(my $i = $break ; $i>$break-10 ; $i--){
			               	my $mNeoReadId = $readId."_".$i."_".$fullReadLength."/1" ;
          					if(not defined $passed->{$mNeoReadId}){
				               	print $list $mNeoReadId."\n" ;
          					}
          					$passed->{$mNeoReadId} = 1 ;
          				}

		               	# load $self->{nbTimeMapped}
		               	if(defined $self->{nbTimeMapped}->{mTE}->{$a->display_name}){
		               		$self->{nbTimeMapped}->{mTE}->{$a->display_name} += 1 ;
		               	}
		               	else{
		               		$self->{nbTimeMapped}->{mTE}->{$a->display_name} = 1 ;
		               	}

		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{subPart}   = 2 ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{featureID} = $teid ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{fam}       = $fam[0] ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{obj}       = $a ;

		               	# exact insertion
		               	if($teObj->start == $a->start){ # read mapped exactly at the start
		               		if($teObj->strand == 1){
				               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mappedInStart} = 1 ;
		               		}
		               		elsif($teObj->strand == 1){
				               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mappedInEnd} = 1 ;
		               		}
			               	$self->{mappedOnStartEnd}->{$readId}->{$splitPosition} = 1 ;
		               	}

		               	if($teObj->strand == 1){
			               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mtExtremity} = 5 ;
		               	}
		               	elsif($teObj->strand == -1){
			               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mtExtremity} = 3 ;
		               	}
					}
	        	}
			});

			#### for 3prime of the TEs.
			$position = $teObj->seq_id.":".($teObj->end-$self->{option}->{maxReadsLength})."-".$teObj->end ;
			$self->{sam}->fetch($position,
		      sub {
		        my $a = shift;
	       		my ($readId, $halfRead, $size, $splitPosition, $break, $fullReadLength) = $self->spliceReadsId($a->display_name) ;
	        	if($halfRead == 1){ # for h1
	        		if($a->strand == "1"){ # if h1 and strand = - => useless

						# selected mNeoReadId
          				for(my $i = $break ; $i<$break+10 ; $i++){
			               	my $mNeoReadId = $readId."_".$i."_".$fullReadLength."/2" ;
          					if(not defined $passed->{$mNeoReadId}){
				               	print $list $mNeoReadId."\n" ;
          					}
          					$passed->{$mNeoReadId} = 1 ;
          				}

		               	# load $self->{nbTimeMapped}
		               	if(defined $self->{nbTimeMapped}->{mTE}->{$a->display_name}){
		               		$self->{nbTimeMapped}->{mTE}->{$a->display_name} += 1 ;
		               	}
		               	else{
		               		$self->{nbTimeMapped}->{mTE}->{$a->display_name} = 1 ;
		               	}

		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{subPart}   = 1 ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{featureID} = $teid ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{fam}       = $fam[0] ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{obj}       = $a ;

		               	# exact insertion
		               	if($teObj->end == $a->end){ # read mapped exactly at the start
		               		if($teObj->strand == 1){
				               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mappedInEnd} = 1 ;
		               		}
		               		elsif($teObj->strand == 1){
				               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mappedInStart} = 1 ;
		               		}
			               	$self->{mappedOnStartEnd}->{$readId}->{$splitPosition} = 1 ;
		               	}
		               	if($teObj->strand == 1){
			               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mtExtremity} = 3 ;
		               	}
		               	elsif($teObj->strand == -1){
			               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mtExtremity} = 5 ;
		               	}
					}
	        	}
	        	elsif($halfRead == 2){ # for h2
	        		if($a->strand == "-1"){ # if h2 and strand = + => useless

						# selected mNeoReadId
          				for(my $i = $break ; $i>$break-10 ; $i--){
			               	my $mNeoReadId = $readId."_".$i."_".$fullReadLength."/1" ;
          					if(not defined $passed->{$mNeoReadId}){
				               	print $list $mNeoReadId."\n" ;
          					}
          					$passed->{$mNeoReadId} = 1 ;
          				}

		               	# load $self->{nbTimeMapped}
		               	if(defined $self->{nbTimeMapped}->{mTE}->{$a->display_name}){
		               		$self->{nbTimeMapped}->{mTE}->{$a->display_name} += 1 ;
		               	}
		               	else{
		               		$self->{nbTimeMapped}->{mTE}->{$a->display_name} = 1 ;
		               	}

		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{subPart}   = 2 ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{featureID} = $teid ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{fam}       = $fam[0] ;
		               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{obj}       = $a ;

		               	# exact insertion
		               	if($teObj->end == $a->end){ # read mapped exactly at the start
		               		if($teObj->strand == 1){
				               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mappedInEnd} = 1 ;
		               		}
		               		elsif($teObj->strand == 1){
				               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mappedInStart} = 1 ;
		               		}
			               	$self->{mappedOnStartEnd}->{$readId}->{$splitPosition} = 1 ;
		               	}
		               	if($teObj->strand == 1){
			               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mtExtremity} = 3 ;
		               	}
		               	elsif($teObj->strand == -1){
			               	$self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$a->display_name}."--".$a->display_name}->{mtExtremity} = 5 ;
		               	}
					}
	        	}
	      });
	}
}
sub splitReadsGrab_v2{
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: splitReadsGrab !\n" ;

	my $self = shift;
	my $sam ;
	my @features;

	my $fileName = $self->{option}->{unmappedFastqFile}->{prefix}.".mateSplitReads.bam" ;
	if(-e $fileName){
		print STDERR "INFO $prog: File $fileName already exist => skip bamutils filter program.\n" ;
	}
	else{
		my @bamutils = ($self->{program}->{bamutils}, 'filter ', $self->{fileName}->{mappingOutput}->{bam}, $fileName, '-whitelist', $self->{option}->{unmappedFastqFile}->{prefix}.'.fishReads.list') ;
		print STDERR join ( " ", @bamutils),"\n" ;
		system(join" ", @bamutils) ;

		my @index = ($self->{program}->{samtools}, 'index', $fileName) ;
		print STDERR join ( " ", @index),"\n" ;
		system(join" ", @index) ;
	}
	
	$sam = Bio::DB::Sam->new(-bam => $fileName);
	@features = $sam->features("match") ; # output all matching reads.

	SPLITREAD:for my $a (@features) {

  		my ($readId, $halfRead, $size, $splitPosition) = $self->spliceReadsId($a->display_name) ;
		my $mTEhalfRead ;
		if($halfRead == 1){
			$mTEhalfRead = 2 ;
		}
		elsif($halfRead == 2){
			$mTEhalfRead = 1 ;
		}
   		# for splitReads having their mTE mapped on the edge of the TEs
   		if(defined $self->{mappedOnStartEnd}->{$readId}){
       		if(defined $self->{mappedOnStartEnd}->{$readId}->{$splitPosition}){
       			delete($self->{remember}->{$readId}) ;
				if(not defined $self->{nbTimeMapped}->{mNeo}->{$a->display_name}){
					$self->{nbTimeMapped}->{mNeo}->{$a->display_name} = 1 ;
				}
				else{
					$self->{nbTimeMapped}->{mNeo}->{$a->display_name} += 1 ;
				}
				$self->{reads}->{mNeo}->{$self->{nbTimeMapped}->{mNeo}->{$a->display_name}."--".$a->display_name}->{obj} = $a ;

				# remove all type of mTE reads
				my($split, $readLength) = ($splitPosition) =~/_([0-9]+)_([0-9]+)/ ;
				for(my $i = $self->{option}->{lengthPairedEnd}[0] ; $i<=$readLength-$self->{option}->{lengthPairedEnd}[0]+1 ; $i++){
					if($i == $split){ next ;} # keep the mTE read that is mapped to the Start or End
					if(defined $self->{nbTimeMapped}->{mTE}->{$readId."_".$i."_".$readLength."/".$mTEhalfRead}){
						for(my $j = 1 ; $j <= $self->{nbTimeMapped}->{mTE}->{$readId."_".$i."_".$readLength."/".$mTEhalfRead}; $j++){
							if(defined $self->{reads}->{mTE}->{$j."--".$readId."_".$i."_".$readLength."/".$mTEhalfRead}){
								# delete($self->{reads}->{mTE}->{$j."--".$readId."_".$i."_".$readLength."/".$mTEhalfRead}) ;
							}
						}
					}
				}
				next SPLITREAD ;
       		}
       		else{ # do not push the reads.
				# if one split-read variants already been part of $self->{reads}->{mNeo} no need to deal with other split-read variants
				foreach my $keys (sort keys %{$self->{mappedOnStartEnd}->{$readId}}){
					if(defined $self->{reads}->{mNeo}->{"1--".$readId.$keys."/1"}){
						next SPLITREAD ;
					}
					elsif(defined $self->{reads}->{mNeo}->{"1--".$readId.$keys."/2"}){
						next SPLITREAD ;
					}
				} ;

       			# push the longest reads
       			if(defined $self->{remember}->{$readId}){
	       			if($a->length > $self->{remember}->{$readId}->length){
	       				$self->{remember}->{$readId} = $a ;
	       			}
       			}
       			else{
       				$self->{remember}->{$readId} = $a ;
       			}
       			next SPLITREAD ;
       		}
       	}
       	else{ # select base on the size
	   		if(not defined $self->{selectBestSplit}->{$readId}->{$halfRead}){ # first time program read this readId
	       		$self->{selectBestSplit}->{$readId}->{$halfRead}->{size} = $size ;
	       		$self->{selectBestSplit}->{$readId}->{$halfRead}->{fullId} = $a->display_name ;
       			# store read
				if(not defined $self->{nbTimeMapped}->{mNeo}->{$a->display_name}){
					$self->{nbTimeMapped}->{mNeo}->{$a->display_name} = 1 ;
				}
				else{
					$self->{nbTimeMapped}->{mNeo}->{$a->display_name} += 1 ;
				}
				$self->{reads}->{mNeo}->{$self->{nbTimeMapped}->{mNeo}->{$a->display_name}."--".$a->display_name}->{obj} = $a ;
       		}
       		elsif(defined $self->{selectBestSplit}->{$readId}->{$halfRead}){
       			if($size < $self->{selectBestSplit}->{$readId}->{$halfRead}->{size}){ # new splitreads is smaller than previous -> dont store new split reads + delete hTE
       				# here delete mTE, used $self->{nbTimeMapped}->{$a->display_name} in a for loop
					# for(my $j = 1 ; $j <= $self->{nbTimeMapped}->{mTE}->{$readId.$splitPosition."/".$mTEhalfRead} ; $j++){
					# 	if(defined $self->{reads}->{mTE}->{$j."--".$readId.$splitPosition."/".$mTEhalfRead}){
							# delete($self->{reads}->{mTE}->{$j."--".$readId.$splitPosition."/".$mTEhalfRead}) ;
						# }
					# }
       			}
       			elsif($size > $self->{selectBestSplit}->{$readId}->{$halfRead}->{size}){ # new splitreads is longer than previous
					if(not defined $self->{nbTimeMapped}->{mNeo}->{$a->display_name}){
						$self->{nbTimeMapped}->{mNeo}->{$a->display_name} = 1 ;
					}
					else{
						$self->{nbTimeMapped}->{mNeo}->{$a->display_name} += 1 ;
					}
					$self->{reads}->{mNeo}->{$self->{nbTimeMapped}->{mNeo}->{$a->display_name}."--".$a->display_name}->{obj} = $a ;

       				# remove mNeo shorter version
					delete($self->{reads}->{mNeo}->{$self->{nbTimeMapped}->{mNeo}->{$self->{selectBestSplit}->{$readId}->{$halfRead}->{fullId}}."--".$self->{selectBestSplit}->{$readId}->{$halfRead}->{fullId}}) ;

        			# remove mTE
		       		my ($prefReadId, $prevHalfRead, $prevSize, $prevSplitPosition) = $self->spliceReadsId($self->{selectBestSplit}->{$readId}->{$halfRead}->{fullId}) ;
					# delete($self->{reads}->{mTE}->{$self->{nbTimeMapped}->{mTE}->{$prefReadId.$prevSplitPosition."/".$mTEhalfRead}."--".$prefReadId.$prevSplitPosition."/".$mTEhalfRead}) ; # this should not work
		       		
		       		# push curent mNeo read in selectBestSplit
		       		$self->{selectBestSplit}->{$readId}->{$halfRead}->{size} = $size ;
		       		$self->{selectBestSplit}->{$readId}->{$halfRead}->{fullId} = $a->display_name ;
       			}
       			elsif($size == $self->{selectBestSplit}->{$readId}->{$halfRead}->{size}){
					if(not defined $self->{nbTimeMapped}->{mNeo}->{$a->display_name}){
						$self->{nbTimeMapped}->{mNeo}->{$a->display_name} = 1 ;
					}
					else{
						$self->{nbTimeMapped}->{mNeo}->{$a->display_name} += 1 ;
					}
					$self->{reads}->{mNeo}->{$self->{nbTimeMapped}->{mNeo}->{$a->display_name}."--".$a->display_name}->{obj} = $a ;
       				# new size, new fullId
		       		$self->{selectBestSplit}->{$readId}->{$halfRead}->{size} = $size ;
		       		$self->{selectBestSplit}->{$readId}->{$halfRead}->{fullId} = $a->display_name ;
       			}
       		}
       	}
	}
	foreach my $read (sort keys %{$self->{remember}}){
		if(not defined $self->{reads}->{mNeo}->{"1--".$self->{remember}->{$read}->display_name}){
			$self->{reads}->{mNeo}->{"1--".$self->{remember}->{$read}->display_name}->{obj} = $self->{remember}->{$read} ;
			$self->{nbTimeMapped}->{mNeo}->{$self->{remember}->{$read}->display_name} = 1 ;
		}
		else{
			$self->{reads}->{mNeo}->{$self->{nbTimeMapped}->{$self->{remember}->{$read}->display_name}."--".$self->{remember}->{$read}->display_name}->{obj} = $self->{remember}->{$read} ;
			$self->{nbTimeMapped}->{mNeo}->{$self->{remember}->{$read}->display_name} += 1 ;
		}
	}
	delete($self->{selectBestSplit}) ;
	delete($self->{remember}) ;

	# foreach my $read (sort keys %{$self->{reads}->{mNeo}}){
	# 	my $mnObj = $self->{reads}->{mNeo}->{$read}->{obj} ;
	# 	if($mnObj->display_name =~/233949_Chr2:2860178-2860262/){
	# 		print Dumper $mnObj->display_name ;
	# 		print Dumper $mnObj->seq_id ;
	# 		print Dumper $mnObj->start ;
	# 		exit ;
	# 	}
	# }
	# exit ;
}
sub groupOverlapingmNeoReads{
	my $self = shift;
	my $datestring = localtime();
	my $nbReads = scalar(keys %{$self->{reads}->{mNeo}}) ;
	print STDERR "INFO $prog $datestring Run Module: groupOverlapingmNeoReads, $nbReads to process.\n" ;
	foreach my $mNeoReadId (
	sort {
		$self->{reads}->{mNeo}->{$a}->{obj}->seq_id cmp $self->{reads}->{mNeo}->{$b}->{obj}->seq_id or
		$self->{reads}->{mNeo}->{$a}->{obj}->start <=> $self->{reads}->{mNeo}->{$b}->{obj}->start
	}
	keys %{$self->{reads}->{mNeo}}){
  		my ($readId, $mNeoHalfRead, $size, $splitPosition, $break, $fullReadLength) = $self->spliceReadsId($self->{reads}->{mNeo}->{$mNeoReadId}->{obj}->display_name) ;
		my $mTEhalfRead ;
		if($mNeoHalfRead == 1){
			$mTEhalfRead = 2 ;
		}
		elsif($mNeoHalfRead == 2){
			$mTEhalfRead = 1 ;
		}
		my $mTEReadId = "1--".$readId.$splitPosition."/".$mTEhalfRead ; # here choose randomly 1-- could be immprove.

		my $readObj = $self->{reads}->{mNeo}->{$mNeoReadId}->{obj} ;
		if(not defined $self->{reads}->{mTE}->{$mTEReadId}){
			if($mNeoHalfRead == 1){
				for(my $i = $break ; $i<$break+10 ; $i++){
					my $alias_mTEReadId = $readId."_".$i."_".$fullReadLength."/2" ;
					if(defined $self->{reads}->{mTE}->{"1--".$alias_mTEReadId}){
						$self->{nbTimeMapped}->{mTE}->{$readId.$splitPosition."/".$mTEhalfRead} = $self->{nbTimeMapped}->{mTE}->{$alias_mTEReadId} ;
						if(not defined $self->{nbTimeMapped}->{mTE}->{$alias_mTEReadId}){
							next ;
						}
						for(my $j = 1 ; $j<= $self->{nbTimeMapped}->{mTE}->{$alias_mTEReadId} ; $j++){ # for mulitmapped reads
							$self->{reads}->{mTE}->{$j."--".$readId.$splitPosition."/".$mTEhalfRead} = $self->{reads}->{mTE}->{$j."--".$alias_mTEReadId} ;
						}
						last ;
					}
				}
			}
			elsif($mNeoHalfRead == 2){
				for(my $i = $break ; $i>$break-10 ; $i--){
					my $alias_mTEReadId = $readId."_".$i."_".$fullReadLength."/1" ;
					if(defined $self->{reads}->{mTE}->{"1--".$alias_mTEReadId}){
						$self->{nbTimeMapped}->{mTE}->{$readId.$splitPosition."/".$mTEhalfRead} = $self->{nbTimeMapped}->{mTE}->{$alias_mTEReadId} ;
						if(not defined $self->{nbTimeMapped}->{mTE}->{$alias_mTEReadId}){
							next ;
						}
						for(my $j = 1 ; $j<= $self->{nbTimeMapped}->{mTE}->{$alias_mTEReadId} ; $j++){ # for mulitmapped reads
							$self->{reads}->{mTE}->{$j."--".$readId.$splitPosition."/".$mTEhalfRead} = $self->{reads}->{mTE}->{$j."--".$alias_mTEReadId} ;
						}
						last ; 
					}
				}
			}
		}
		if(defined $self->{reads}->{mTE}->{$mTEReadId}){
			$self->groupOverlappingReads_v2($self->{reads}->{mNeo}->{$mNeoReadId}->{obj}, $self->{reads}->{mTE}->{$mTEReadId}->{featureID}, $mNeoReadId) ;
			$self->groupOverlappingReads_v2($self->{reads}->{mNeo}->{$mNeoReadId}->{obj}, $self->{reads}->{mTE}->{$mTEReadId}->{fam}, $mNeoReadId) ;
		}
	}
}
sub groupOverlappingReads_v2 {
	my $self = shift;
	my $readObj = shift ;
	my $featureID = shift ;
	my $readId = shift ;

	if(defined $self->{GroupPosition}->{$featureID}->{$readObj->seq_id}->{$readObj->start}){
		my $groupK = $self->{GroupPosition}->{$featureID}->{$readObj->seq_id}->{$readObj->start} ;
		my ($id, $subPart, $size, $splitPosition) = $self->spliceReadsId($readId) ;
		$self->{group}->{$featureID}->{$groupK}->{readsID}->{$readId}->{subPart} = $subPart ;
		$self->{group}->{$featureID}->{$groupK}->{readsID}->{$readId}->{size} = $size ;
		$self->{group}->{$featureID}->{$groupK}->{readsID}->{$readId}->{id} = $id.$splitPosition ;

		for(my $i = $self->{group}->{$featureID}->{$groupK}->{stop} ; $i <= $readObj->end ; $i++){
			$self->{GroupPosition}->{$featureID}->{$readObj->seq_id}->{$i} = $groupK ;
		}
		$self->{group}->{$featureID}->{$groupK}->{stop} = $readObj->end ; # new stop of the group
	}
	elsif(defined $self->{GroupPosition}->{$featureID}->{$readObj->seq_id}->{$readObj->end}){
		my $groupK = $self->{GroupPosition}->{$featureID}->{$readObj->seq_id}->{$readObj->end} ;
		my ($id, $subPart, $size, $splitPosition) = $self->spliceReadsId($readId) ;
		$self->{group}->{$featureID}->{$groupK}->{readsID}->{$readId}->{subPart} = $subPart ;
		$self->{group}->{$featureID}->{$groupK}->{readsID}->{$readId}->{size} = $size ;
		$self->{group}->{$featureID}->{$groupK}->{readsID}->{$readId}->{id} = $id.$splitPosition;

		for(my $i = $self->{group}->{$featureID}->{$groupK}->{start} ; $i >= $readObj->start ; $i--){
			$self->{GroupPosition}->{$featureID}->{$readObj->seq_id}->{$i} = $groupK ;
		}
		$self->{group}->{$featureID}->{$groupK}->{start} = $readObj->start ; # new start of the group
	}
	else{ # start or end of the read is not part of a group = defining a new group
		$self->{groupK}->{$featureID} += 1 ;
		for(my $i = $readObj->start ; $i <= $readObj->end+1 ; $i++){
			$self->{GroupPosition}->{$featureID}->{$readObj->seq_id}->{$i} = $self->{groupK}->{$featureID} ;
		}
		$self->{group}->{$featureID}->{$self->{groupK}->{$featureID}}->{start} = $readObj->start ;
		$self->{group}->{$featureID}->{$self->{groupK}->{$featureID}}->{stop} = $readObj->end ;
		$self->{group}->{$featureID}->{$self->{groupK}->{$featureID}}->{chr} = $readObj->seq_id ;
		my ($id, $subPart, $size, $splitPosition) = $self->spliceReadsId($readId) ;
		$self->{group}->{$featureID}->{$self->{groupK}->{$featureID}}->{readsID}->{$readId}->{subPart} = $subPart ;
		$self->{group}->{$featureID}->{$self->{groupK}->{$featureID}}->{readsID}->{$readId}->{size} = $size ;
		$self->{group}->{$featureID}->{$self->{groupK}->{$featureID}}->{readsID}->{$readId}->{id} = $id.$splitPosition;
	}
}
sub spliceReadsId{
	my $self = shift;
	my $readId = shift;
	my ($id, $subPart, $size, $splitPosition) ;
	if($readId =~/(.*)_([0-9]+)_([0-9]+)\/([1-2])$/){
		$id = $1 ;
		$splitPosition = "_".$2."_".$3 ;
		$subPart = $4 ;
		if($subPart == 1){ $size = $2 ; }
		elsif($subPart == 2){ $size = $3-$2+1 }
	}
	else{
		print Dumper $readId ; exit ;
		print STDERR "wrong reads ID\n" ;
	}
	return($id, $subPart, $size, $splitPosition,$2,$3) ;
}
sub neoInsertionFinder {
	my $datestring = localtime();
	print STDERR "INFO $prog $datestring Run Module: neoInsertionFinder !\n" ;
	my $self = shift;

	foreach my $featureID (sort keys %{$self->{group}}){
		foreach my $gk (sort keys %{$self->{group}->{$featureID}}){
			if(scalar keys %{$self->{group}->{$featureID}->{$gk}->{readsID}} < 4){ next ; } # if to small amount of reads in group next
			my $k = 0 ;
			my $countNbReads ;
			$self->{group}->{$featureID}->{$gk}->{nbReads}->{5} = 0 ;
			$self->{group}->{$featureID}->{$gk}->{nbReads}->{3} = 0 ;

			# if($self->{group}->{$featureID}->{$gk}->{start} > 2799388 and $self->{group}->{$featureID}->{$gk}->{stop} < 2799798 and $self->{group}->{$featureID}->{$gk}->{chr} eq "Chr2"){
			# 	print Dumper $self->{group}->{$featureID}->{$gk} ;
			# 	print Dumper $gk ;
			# 	print Dumper $featureID ;
			# }

			foreach my $mnReadId (
			sort {$self->{reads}->{mNeo}->{$a}->{obj}->start <=> $self->{reads}->{mNeo}->{$b}->{obj}->start }
			keys %{$self->{group}->{$featureID}->{$gk}->{readsID}}){

				# Get information for mNeo
				$k++ ;
				my $mnStrand        = $self->{reads}->{mNeo}->{$mnReadId}->{obj}->strand ; # mn = mate matching Neo 
				my $mnSubpart       = $self->{group}->{$featureID}->{$gk}->{readsID}->{$mnReadId}->{subPart} ;
				my $splitReadsId    = $self->{group}->{$featureID}->{$gk}->{readsID}->{$mnReadId}->{id} ;
				my $oldConcordanceF = $self->{group}->{$featureID}->{$gk}->{singleReadConcordance}->{F} ;
				my ($readId)        = $splitReadsId =~/^[0-9]+--(.*)_[0-9]+_[0-9]+$/ ;
				my $mnObj           = $self->{reads}->{mNeo}->{$mnReadId}->{obj} ;

				# Get information for mTE
				my $mtSubPart ; 
				my ($mtTEid) = $splitReadsId =~/\d+--(.*_[0-9]+_[0-9]+)$/ ;
				if($self->{group}->{$featureID}->{$gk}->{readsID}->{$mnReadId}->{subPart} == 1){
					$mtSubPart    = "2" ;
					$mtTEid      .="/2" ;
				}
				elsif($self->{group}->{$featureID}->{$gk}->{readsID}->{$mnReadId}->{subPart} == 2){
					$mtSubPart     = "1" ;
					$mtTEid       .="/1" ;
				}

				my $findSignle_mTE ;
				for(my $i = 1 ; $i<= $self->{nbTimeMapped}->{mTE}->{$mtTEid} ; $i++){ # for each single or multi mapped reads.
	    	       	if(not defined $self->{reads}->{mTE}->{$i."--".$mtTEid}){ next READ;}
	    	       	my $mtObj = $self->{reads}->{mTE}->{$i."--".$mtTEid}->{obj} ;
					my $mtExtremity  = $self->{reads}->{mTE}->{$i."--".$mtTEid}->{mtExtremity} ;
					my $teObj = $self->{targetFeat}->{$self->{reads}->{mTE}->{$i."--".$mtTEid}->{featureID}}->{'te'} ;
					if(not defined $mtExtremity){
						next ;
					}

					# Filter #1: Check Concordance: some combination are not possible, if concordance has not been proved the read is filtered out.
					my $concordance = $self->checkSingleReadConcordance($mnSubpart, $mnStrand, $mtSubPart, $mtObj->strand, $mtExtremity, $gk, $teObj->strand) ;
					if($concordance eq "FALSE"){ next ;} # next mtTEid
					my $distance ;
					if($teObj->strand == 1){
						if($mtExtremity == 5){
							$distance = abs($teObj->start - $mtObj->start) ;
						}
						elsif($mtExtremity == 3){
							$distance = abs($teObj->end - $mtObj->end) ;
						}
					}
					elsif($teObj->strand == -1){
						if($mtExtremity == 5){
							$distance = abs($teObj->end - $mtObj->end) ;
						}
						elsif($mtExtremity == 3){
							$distance = abs($teObj->start - $mtObj->start) ;
						}
					}
					$findSignle_mTE->{$mtExtremity}->{$i."--".$mtTEid} = $distance ;
				} # END of for multi mapped mtTEid reads

				# choose one single mTE read
				my $mtObj ;
				my $mtExtremity ;
				$mtTEid = '' ;
				if(defined $findSignle_mTE->{5} and defined $findSignle_mTE->{3}){
					if(scalar(keys %{$findSignle_mTE->{5}}) > scalar(keys %{$findSignle_mTE->{3}})){
						my @array = sort {$findSignle_mTE->{5}->{$a} <=> $findSignle_mTE->{5}->{$b}} keys %{$findSignle_mTE->{5}} ;
						$mtTEid = $array[0] ;
						$mtObj = $self->{reads}->{mTE}->{$array[0]}->{obj} ;
						$mtExtremity = 5 ;
					}
					elsif(scalar(keys %{$findSignle_mTE->{5}}) < scalar(keys %{$findSignle_mTE->{3}})){
						my @array = sort {$findSignle_mTE->{3}->{$a} <=> $findSignle_mTE->{3}->{$b}} keys %{$findSignle_mTE->{3}} ;
						$mtTEid = $array[0] ;
						$mtObj = $self->{reads}->{mTE}->{$array[0]}->{obj} ;
						$mtExtremity = 3 ;
					}
					elsif(scalar(keys %{$findSignle_mTE->{5}}) == scalar(keys %{$findSignle_mTE->{3}})){
						my @array5 = sort {$findSignle_mTE->{5}->{$a} <=> $findSignle_mTE->{5}->{$b}} keys %{$findSignle_mTE->{5}} ;
						my @array3 = sort {$findSignle_mTE->{3}->{$a} <=> $findSignle_mTE->{3}->{$b}} keys %{$findSignle_mTE->{3}} ;
						if($findSignle_mTE->{5}->{$array5[0]} < $findSignle_mTE->{3}->{$array3[0]}){ # keep smallest distance
							$mtTEid = $array5[0] ;
							$mtObj = $self->{reads}->{mTE}->{$array5[0]}->{obj} ;
							$mtExtremity = 5 ;
						}
						elsif($findSignle_mTE->{5}->{$array5[0]} > $findSignle_mTE->{3}->{$array3[0]}){
							$mtTEid = $array3[0] ;
							$mtObj = $self->{reads}->{mTE}->{$array3[0]}->{obj} ;
							$mtExtremity = 3 ;
						}
						else{ next ;}
					}
				}
				elsif(defined $findSignle_mTE->{5} and not defined $findSignle_mTE->{3}){
					my @array = sort {$findSignle_mTE->{5}->{$a} <=> $findSignle_mTE->{5}->{$b}} keys %{$findSignle_mTE->{5}} ;
					$mtTEid = $array[0] ;
					$mtObj = $self->{reads}->{mTE}->{$array[0]}->{obj} ;
					$mtExtremity = 5 ;
				}
				elsif(defined $findSignle_mTE->{3} and not defined $findSignle_mTE->{5} ){
					my @array = sort {$findSignle_mTE->{3}->{$a} <=> $findSignle_mTE->{3}->{$b}} keys %{$findSignle_mTE->{3}} ;
					$mtTEid = $array[0] ;
					$mtObj = $self->{reads}->{mTE}->{$array[0]}->{obj} ;
					$mtExtremity = 3 ;
				}

				# Get info TE
				my $teObj = $self->{targetFeat}->{$self->{reads}->{mTE}->{$mtTEid}->{featureID}}->{'te'} ;
				my @fam   = $self->{targetFeat}->{$self->{reads}->{mTE}->{$mtTEid}->{featureID}}->{'te'}->get_tag_values('fam') ;
				my @sF    = $self->{targetFeat}->{$self->{reads}->{mTE}->{$mtTEid}->{featureID}}->{'te'}->get_tag_values('sF') ;
				$self->{group}->{$featureID}->{$gk}->{teFam} = $fam[0] ;
				$self->{group}->{$featureID}->{$gk}->{tesF} = $sF[0] ;
				$self->{group}->{$featureID}->{$gk}->{teid}->{$self->{reads}->{mTE}->{$mtTEid}->{featureID}} = 1 ;

				# store correct Reads ID
				if($mtExtremity == 5){
					$self->{group}->{$featureID}->{$gk}->{readsList}->{neo}->{5}->{$mnObj->display_name} = 1;
					$self->{group}->{$featureID}->{$gk}->{readsList}->{te}->{5}->{$mtObj->display_name} = 1;
				}
				elsif($mtExtremity == 3){
					$self->{group}->{$featureID}->{$gk}->{readsList}->{neo}->{3}->{$mnObj->display_name} = 1;
					$self->{group}->{$featureID}->{$gk}->{readsList}->{te}->{3}->{$mtObj->display_name} = 1;
				}

				# estimated the length of the overlap between 5' and 3' reads
				$self->findOverlap($gk, $featureID, $mtExtremity, $mnObj) ;

				# find strand of the neoInsertion (here I could set up a detection of strand conflict)
				if(not defined $self->{group}->{$featureID}->{$gk}->{strand}){
					if($mtExtremity == 5){
						$self->{group}->{$featureID}->{$gk}->{strand} = 1 ;
					}
					elsif($mtExtremity == 3){ 
						$self->{group}->{$featureID}->{$gk}->{strand} = -1 ;
					}
				}

				if($self->{group}->{$featureID}->{$gk}->{strand} == 1){
					if($mtExtremity == 5){
						push(@{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}, $mnObj->end);
					}
					elsif($mtExtremity == 3){
						push(@{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}, $mnObj->start);
					}
				}
				elsif($self->{group}->{$featureID}->{$gk}->{strand} == -1){
					if($mtExtremity == 5){
						push(@{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}, $mnObj->start);
					}
					elsif($mtExtremity == 3){
						push(@{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}, $mnObj->end);
					}
				}

				# find position of putative TSD
				$self->findExactInsertion($featureID, $teObj, $gk, $mtTEid, $mnObj, $self->{group}->{$featureID}->{$gk}->{strand}, $mtExtremity, $mtObj) ;

				# analysis group concordance => check group conflict
				$self->loadGroupConflict($gk, $featureID, $mtExtremity) ;

				# if($mnObj->display_name =~/157529_Chr2:351115-351199/){
				# 	print Dumper $mnObj->display_name ;
				# 	print Dumper $mnObj->seq_id ;
				# 	print Dumper $mnObj->start ;
				# 	print Dumper $readId ;
				# 	print Dumper $self->{reads}->{mTE}->{$mtTEid}->{featureID} ;
				# }

				# keep reads for bam output
				$self->{printBam}->{reads}->{$mtObj->seq_id}->{$mtObj->start}->{$readId} = $mtObj ; # push mate mapped on TEs
				$self->{printBam}->{reads}->{$mnObj->seq_id}->{$mnObj->start}->{$readId} = $mnObj ;# push mNeo

				# count number of reads
				if(not defined $countNbReads->{$readId}){
					$self->{group}->{$featureID}->{$gk}->{nbReads}->{overall} +=1 ;
					if($mtExtremity == 5){
						$self->{group}->{$featureID}->{$gk}->{nbReads}->{5} +=1 ;
						$countNbReads->{$readId} = 5 ;
					}
					elsif($mtExtremity == 3){
						$self->{group}->{$featureID}->{$gk}->{nbReads}->{3} +=1 ;
						$countNbReads->{$readId} = 3 ;
					}
				}
				#report if group contain both 5 and 3 prime mtExtremity
				if($mtExtremity == 5){
					$self->{group}->{$featureID}->{$gk}->{typesOfMatch}->{5} = 5;
				}
				elsif($mtExtremity == 3){
					$self->{group}->{$featureID}->{$gk}->{typesOfMatch}->{3} = 3;
				}
			} # end of foreach $mateSplitReadsId

			# Find TSD
			my @array = sort {$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{$b} <=> $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{$a} } keys %{$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}} ;
			$self->{group}->{$featureID}->{$gk}->{exactSite}->{start} = $array[0] ;
			@array = sort {$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{$b} <=> $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{$a} } keys %{$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}} ;
			$self->{group}->{$featureID}->{$gk}->{exactSite}->{stop} = $array[0] ;

			# Filter #1: mNeoReads group should match both TE extremity 5 and 3

			# if($featureID eq "ATENSPM5" and $gk == 10){
			# 	print Dumper $self->{group}->{$featureID}->{$gk} ;
			# 	print Dumper "###########" ;
			# 	exit ;
			# }

			my $mtExtremity ;
			if(defined $self->{group}->{$featureID}->{$gk}->{typesOfMatch}->{5} and defined $self->{group}->{$featureID}->{$gk}->{typesOfMatch}->{3}){
				$mtExtremity = "5/3" ;
			}
			else{
				next ;
			}

			# Filter #2: Number of reads is the group
			if($self->{group}->{$featureID}->{$gk}->{nbReads}->{5} < 2 or $self->{group}->{$featureID}->{$gk}->{nbReads}->{3} < 2){
				if(not defined $debug){
					next ;
				}
			}
			elsif($self->{group}->{$featureID}->{$gk}->{nbReads}->{overall} < 5){
				if(not defined $debug){
					next ;
				}
			}

			# Filter #3: dispersion of position
			if(scalar(@{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}) < 5){ # calc dispersion
				my $dispersion = $self->{group}->{$featureID}->{$gk}->{pos}->{5}[scalar(@{$self->{group}->{$featureID}->{$gk}->{pos}->{5}})-1] - $self->{group}->{$featureID}->{$gk}->{pos}->{5}[0] ;
				if($dispersion>5){
					if(not defined $debug){
						next ;
					}
				}
				$self->{group}->{$featureID}->{$gk}->{dispersion}->{5} = $dispersion ;					
			}
			else{ # calc STDV
				my $statObj = Statistics::Descriptive::Full->new();
				$statObj->add_data(@{$self->{group}->{$featureID}->{$gk}->{pos}->{5}});
				if($statObj->standard_deviation>5){
					if(not defined $debug){
						next ;
					}
				}
				$self->{group}->{$featureID}->{$gk}->{standard_deviation}->{5} = $statObj->standard_deviation ;
			}
			if(scalar(@{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}) < 5){ # calc dispersion 3'
				my $dispersion = $self->{group}->{$featureID}->{$gk}->{pos}->{3}[scalar(@{$self->{group}->{$featureID}->{$gk}->{pos}->{3}})-1] - $self->{group}->{$featureID}->{$gk}->{pos}->{3}[0] ;
				if($dispersion>5){
					if(not defined $debug){
						next ;
					}
				}
				$self->{group}->{$featureID}->{$gk}->{dispersion}->{3} = $dispersion ;
			}
			else{ # calc STDV
				my $statObj = Statistics::Descriptive::Full->new();
				$statObj->add_data(@{$self->{group}->{$featureID}->{$gk}->{pos}->{3}});
				if($statObj->standard_deviation>5){
					if(not defined $debug){
						next ;
					}
				}
				$self->{group}->{$featureID}->{$gk}->{standard_deviation}->{3} = $statObj->standard_deviation ;
			}

			# Filter #4: check group Conflict only for LTR
			if($self->{group}->{$featureID}->{$gk}->{tesF} =~/LTR/){ 
				$self->checkGroupConflict($gk, $featureID) ;
				if(0 == $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{nbSwitch} % 2){
					if($self->{group}->{$featureID}->{$gk}->{groupConcordance}->{nbSwitch} > 2){
						if($self->{group}->{$featureID}->{$gk}->{tesF} =~/LTR/){
							if( $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{stdv} > $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{avg} ){
								$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{test} = "conflict" ;
							}
						}
						else{
							$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{test} = "conflict" ;
						}
					}
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{test} = "conflict" ;
				}
				if(defined $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{test}){
					next ;
				}
			}

			# Filter #5: based on overlap
			# 5.1 based on lenght ov: keep group having ov < 25
			if(defined $self->{group}->{$featureID}->{$gk}->{overlap}->{1} and defined $self->{group}->{$featureID}->{$gk}->{overlap}->{2}){
				$self->{group}->{$featureID}->{$gk}->{overlap}->{ov} = $self->{group}->{$featureID}->{$gk}->{overlap}->{1}-$self->{group}->{$featureID}->{$gk}->{overlap}->{2}+1 ;
				if($self->{group}->{$featureID}->{$gk}->{overlap}->{ov} > 25){
					if(not defined $debug){
						next ;
					}
				}
			}
			else{
				if(not defined $debug){
					next ;
				}
			}
			# 5.2 based on ratio grpLenght / ov > 0.5 
			my $grpLenght = $self->{group}->{$featureID}->{$gk}->{stop} - $self->{group}->{$featureID}->{$gk}->{start} + 1 ;
			$self->{group}->{$featureID}->{$gk}->{tmpRatio} = sprintf "%.2f", $self->{group}->{$featureID}->{$gk}->{overlap}->{ov}/$grpLenght ;
			if($self->{group}->{$featureID}->{$gk}->{overlap}->{ov}/$grpLenght > 0.5){
				if(not defined $debug){
					next ;
				}
			}
			# 5.3 # of ov reads
			$self->nbOvReads($featureID, $gk) ;
			$self->{group}->{$featureID}->{$gk}->{nbOvReads}->{5} = sprintf "%.2f", ($self->{group}->{$featureID}->{$gk}->{nbOvReads}->{5}/$self->{group}->{$featureID}->{$gk}->{nbReads}->{5})*100 ;
			$self->{group}->{$featureID}->{$gk}->{nbOvReads}->{3} = sprintf "%.2f", ($self->{group}->{$featureID}->{$gk}->{nbOvReads}->{3}/$self->{group}->{$featureID}->{$gk}->{nbReads}->{3})*100 ;
			if( ($self->{group}->{$featureID}->{$gk}->{nbOvReads}->{5}/$self->{group}->{$featureID}->{$gk}->{nbReads}->{5}) < 0.5 or ($self->{group}->{$featureID}->{$gk}->{nbOvReads}->{3}/$self->{group}->{$featureID}->{$gk}->{nbReads}->{3})<0.5){
				if(not defined $debug){
					next ;
				}
			}

			# Filter #6: coverage
			$self->nbReadsWhitinGroupPosition($gk, $featureID) ;
			if( ($self->{group}->{$featureID}->{$gk}->{covStep2}->{expected} * 2) < $self->{group}->{$featureID}->{$gk}->{covStep2}->{observed}){
				if(not defined $debug){
					next ;
				}
			}

			# Filter #7: group should not be overlaping any copie of the same family.
			my $test = $self->isInSameFamily($gk, $featureID) ;
			if($test eq "TRUE"){
				next ;
			}

			# if($featureID eq "AT5TE39825" and $gk == 6){
			# 	print Dumper $self->{group}->{$featureID}->{$gk} ;
			# 	exit ;
			# }

			my $kPrint = scalar(keys%{$self->{print}})+1 ;
			$self->{print}->{$kPrint}->{chr} = $self->{group}->{$featureID}->{$gk}->{chr} ;
			$self->{print}->{$kPrint}->{start} = $self->{group}->{$featureID}->{$gk}->{start} ;
			$self->{print}->{$kPrint}->{stop} = $self->{group}->{$featureID}->{$gk}->{stop} ;
			$self->{print}->{$kPrint}->{TSDstart} = $self->{group}->{$featureID}->{$gk}->{exactSite}->{start} ;
			$self->{print}->{$kPrint}->{TSDstop} = $self->{group}->{$featureID}->{$gk}->{exactSite}->{stop} ;
			$self->{print}->{$kPrint}->{strand} = $self->{group}->{$featureID}->{$gk}->{strand} ;
			$self->{print}->{$kPrint}->{teid} = $self->{group}->{$featureID}->{$gk}->{teid} ;
			$self->{print}->{$kPrint}->{featureID} = $featureID ;
			$self->{print}->{$kPrint}->{family} = $self->{group}->{$featureID}->{$gk}->{teFam} ;
			$self->{print}->{$kPrint}->{sF} = $self->{group}->{$featureID}->{$gk}->{tesF} ;
			$self->{print}->{$kPrint}->{mtExtremity} = $mtExtremity ;
			$self->{print}->{$kPrint}->{overlap} = $self->{group}->{$featureID}->{$gk}->{overlap}->{ov} ;
			$self->{print}->{$kPrint}->{nbOverallReads} = $self->{group}->{$featureID}->{$gk}->{nbReads}->{overall} ;
			$self->{print}->{$kPrint}->{nbR5} = $self->{group}->{$featureID}->{$gk}->{nbReads}->{5} ;
			$self->{print}->{$kPrint}->{nbR3} = $self->{group}->{$featureID}->{$gk}->{nbReads}->{3} ;
			$self->{print}->{$kPrint}->{readsList} = $self->{group}->{$featureID}->{$gk}->{readsList} ;

			$self->methylation($kPrint, $gk, $featureID) ;

		} # enf of foreach group ID
	} # end of foreach featureID
	delete($self->{group}) ;
	delete($self->{groupK}) ;
	delete($self->{GroupPosition});
	delete($self->{reads});
	delete($self->{nbTimeMapped});
	delete($self->{mappedOnStartEnd});
}
sub isInSameFamily{
	my $self = shift ;
	my ($gk, $featureID) = @_ ;
	my $chr   = $self->{group}->{$featureID}->{$gk}->{chr} ;
	my $start = $self->{group}->{$featureID}->{$gk}->{start} ;
	my $stop  = $self->{group}->{$featureID}->{$gk}->{stop} ;
	my $fam   = $self->{group}->{$featureID}->{$gk}->{teFam} ;
	my $test  = "FALSE" ;

    if(defined $self->{bedTree}->{$fam}->{$chr}){ # read Chr different than copies Chr
		my $overlapping = $self->{bedTree}->{$fam}->{$chr}->fetch( $start, $stop );
		if(defined $overlapping and @$overlapping){
			$test  = "TRUE" ;
		}
    }
    return($test) ;
}
sub loadGroupConflict{
	my $self = shift ;
	my ($gk, $featureID, $mtExtremity) = @_ ;
	my $k ;
	if(not defined $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}){
		$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}->{1}->{$mtExtremity} += 1 ;
		return(0) ;
	}
	else{
		$k = scalar(keys %{$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}}) + 1 ;
	}
	foreach my $lastExt (keys %{$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}->{$k-1}}){
		if($lastExt eq $mtExtremity){
			$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}->{$k-1}->{$mtExtremity} += 1 ;
		}
		else{
			$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}->{$k}->{$mtExtremity} += 1 ;
		}
	}
}
sub checkGroupConflict{
	my $self = shift ;
	my ($gk, $featureID) = @_ ;
	my $nbSwitch = scalar(keys %{$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}})  ;
	my @stat ;
	foreach my $key (keys %{$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}}){
		foreach my $ext (keys %{$self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}->{$key}}){
			push(@stat, $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{store}->{$key}->{$ext}) ;
		}
	}
	my $statObj = Statistics::Descriptive::Full->new();
	$statObj->add_data(@stat);
    $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{nbSwitch} = $nbSwitch;
    $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{stdv} = sprintf "%.8f", $statObj->standard_deviation;
    $self->{group}->{$featureID}->{$gk}->{groupConcordance}->{avg} = sprintf "%.8f", $statObj->mean;
}
sub checkSingleReadConcordance {
	my $self = shift ;
	my ($mnSubpart, $mnStrand, $mtSubPart, $mtStrand, $mtExtremity, $gk, $testrand) = @_ ;
	my $concordance = "TRUE" ;
	# $mtExtremity == 5
	if($testrand==1 and $mtExtremity==5 and $mtSubPart == 1 and $mtStrand == 1 and $mnSubpart == 2 and $mnStrand == 1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==1 and $mtExtremity==5 and $mtSubPart == 1 and $mtStrand == 1 and $mnSubpart == 2 and $mnStrand == -1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==1 and $mtExtremity==5 and $mtSubPart == 2 and $mtStrand == -1 and $mnSubpart == 1 and $mnStrand == 1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==1 and $mtExtremity==5 and $mtSubPart == 2 and $mtStrand == -1 and $mnSubpart == 1 and $mnStrand == -1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==-1 and $mtExtremity==5 and $mtSubPart == 2 and $mtStrand == 1 and $mnSubpart == 1 and $mnStrand == 1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==-1 and $mtExtremity==5 and $mtSubPart == 2 and $mtStrand == 1 and $mnSubpart == 1 and $mnStrand == -1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==-1 and $mtExtremity==5 and $mtSubPart == 1 and $mtStrand == -1 and $mnSubpart == 2 and $mnStrand == 1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==-1 and $mtExtremity==5 and $mtSubPart == 1 and $mtStrand == -1 and $mnSubpart == 2 and $mnStrand == -1){
		$concordance = "FALSE" ;
	}
	###3
	elsif($testrand==1 and $mtExtremity==3 and $mtSubPart == 2 and $mtStrand == 1 and $mnSubpart == 1 and $mnStrand == 1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==1 and $mtExtremity==3 and $mtSubPart == 2 and $mtStrand == 1 and $mnSubpart == 1 and $mnStrand == -1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==1 and $mtExtremity==3 and $mtSubPart == 1 and $mtStrand == -1 and $mnSubpart == 2 and $mnStrand == 1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==1 and $mtExtremity==3 and $mtSubPart == 1 and $mtStrand == -1 and $mnSubpart == 2 and $mnStrand == -1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==-1 and $mtExtremity==3 and $mtSubPart == 1 and $mtStrand == 1 and $mnSubpart == 2 and $mnStrand == 1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==-1 and $mtExtremity==3 and $mtSubPart == 1 and $mtStrand == 1 and $mnSubpart == 2 and $mnStrand == -1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==-1 and $mtExtremity==3 and $mtSubPart == 2 and $mtStrand == -1 and $mnSubpart == 1 and $mnStrand == 1){
		$concordance = "FALSE" ;
	}
	elsif($testrand==-1 and $mtExtremity==3 and $mtSubPart == 2 and $mtStrand == -1 and $mnSubpart == 1 and $mnStrand == -1){
		$concordance = "FALSE" ;
	}
	return($concordance) ;
}
sub findOverlap{
	my $self = shift ;
	my ($gk, $featureID, $mtExtremity, $mnObj) = @_ ;
	if(not defined $self->{group}->{$featureID}->{$gk}->{overlap}->{oldMatch}){
		$self->{group}->{$featureID}->{$gk}->{overlap}->{oldMatch} = $mtExtremity ;
		$self->{group}->{$featureID}->{$gk}->{overlap}->{1} = $mnObj->end+1 ;
	}
	else{
		if($self->{group}->{$featureID}->{$gk}->{overlap}->{oldMatch} == $mtExtremity
			and defined $self->{group}->{$featureID}->{$gk}->{overlap}->{1}
			and not defined $self->{group}->{$featureID}->{$gk}->{overlap}->{2}){
			if($self->{group}->{$featureID}->{$gk}->{overlap}->{1}<= $mnObj->end){
				$self->{group}->{$featureID}->{$gk}->{overlap}->{1} = $mnObj->end+1 ;
			}
		}
		elsif(not defined $self->{group}->{$featureID}->{$gk}->{overlap}->{2}){
			$self->{group}->{$featureID}->{$gk}->{overlap}->{2} = $mnObj->start ;
		}
		$self->{group}->{$featureID}->{$gk}->{overlap}->{oldMatch} = $mtExtremity ;
	}
}
sub nbOvReads{
	my $self = shift ;
	my ($featureID, $gk) = @_ ;
	if($self->{group}->{$featureID}->{$gk}->{strand} == 1){
		POS:foreach my $pos5 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}){
			foreach my $pos3 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}){
				if($pos3<=$pos5+1){
					$self->{group}->{$featureID}->{$gk}->{nbOvReads}->{5}++ ;
					next POS ;
				}
			}
		}
		POS:foreach my $pos3 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}){
			foreach my $pos5 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}){
				if($pos3<=$pos5+1){
					$self->{group}->{$featureID}->{$gk}->{nbOvReads}->{3}++ ;
					next POS ;
				}
			}
		}
	}
	elsif($self->{group}->{$featureID}->{$gk}->{strand} == -1){
		POS:foreach my $pos3 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}){
			foreach my $pos5 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}){
				if($pos5<=$pos3+1){
					$self->{group}->{$featureID}->{$gk}->{nbOvReads}->{3}++ ;
					next POS ;
				}
			}
		}
		POS:foreach my $pos5 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}){
			foreach my $pos3 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}){
				if($pos5<=$pos3+1){
					$self->{group}->{$featureID}->{$gk}->{nbOvReads}->{5}++ ;
					next POS ;
				}
			}
		}
	}
}
sub nbReadsWhitinGroupPosition{
	my $self  = shift ;
	my ($gk, $featureID) = @_ ;
	my $observed = 0 ;

	my $nbReads = 0 ;
	if($self->{group}->{$featureID}->{$gk}->{strand} == 1){
		foreach my $pos5 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}){
			$nbReads += $pos5 - $self->{group}->{$featureID}->{$gk}->{overlap}->{2} +1 ;
		}
		foreach my $pos3 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}){
			$nbReads += $self->{group}->{$featureID}->{$gk}->{overlap}->{1} - $pos3 +1 ;
		}
	}
	elsif($self->{group}->{$featureID}->{$gk}->{strand} == -1){
		foreach my $pos5 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{5}}){
			$nbReads += $self->{group}->{$featureID}->{$gk}->{overlap}->{1} - $pos5 + 1 ;
		}
		foreach my $pos3 (sort {$a <=> $b} @{$self->{group}->{$featureID}->{$gk}->{pos}->{3}}){
			$nbReads += $pos3 - $self->{group}->{$featureID}->{$gk}->{overlap}->{2} +1 ;
		}
	}

	my $segment = $self->{sam}->segment($self->{group}->{$featureID}->{$gk}->{chr},$self->{group}->{$featureID}->{$gk}->{overlap}->{2},$self->{group}->{$featureID}->{$gk}->{overlap}->{1});
	my ($coverage)       = $segment->features('coverage');
	my @data_points      = $coverage->coverage;
	my $statObj = Statistics::Descriptive::Full->new();
	$statObj->add_data(@data_points);

	$self->{group}->{$featureID}->{$gk}->{covStep2}->{expected} = $nbReads ;
	$self->{group}->{$featureID}->{$gk}->{covStep2}->{observed} = $statObj->mean ;
	return(0) ;
}
sub findExactInsertion{
	my $self  = shift ;
	my ($featureID, $teObj, $gk, $mtTEid, $mnObj, $groupStrand, $mtExtremity, $mtObj) = @_ ;

	if($groupStrand == 1){
		if(defined $self->{targetFeat}->{$featureID}->{reads}->{$mtTEid}->{mappedInStart}){
			if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{$mnObj->end}){
				$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{$mnObj->end} += 5 ;
			}
			else{
				$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{$mnObj->end} = 5 ;
			}
		}
		elsif(defined $self->{targetFeat}->{$featureID}->{reads}->{$mtTEid}->{mappedInEnd}){
			if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{$mnObj->start}){
				$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{$mnObj->start} += 5 ;
			}
			else{
				$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{$mnObj->start} = 5 ;
			}
		}
		else{
			$self->readThatDontMapAtTheTEExtremity($teObj, $gk, $featureID, $mtTEid, $mnObj, $groupStrand, $mtExtremity, $mtObj) ;
		}
	}
	elsif($groupStrand == -1){
		if(defined $self->{targetFeat}->{$featureID}->{reads}->{$mtTEid}->{mappedInStart}){
			if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{$mnObj->start}){
				$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{$mnObj->start} += 5 ;
			}
			else{
				$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{$mnObj->start} = 5 ;
			}
		}
		elsif(defined $self->{targetFeat}->{$featureID}->{reads}->{$mtTEid}->{mappedInEnd}){
			if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{$mnObj->end}){
				$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{$mnObj->end} += 5 ;
			}
			else{
				$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{$mnObj->end} = 5 ;
			}
		}
		else{
			$self->readThatDontMapAtTheTEExtremity($teObj, $gk, $featureID, $mtTEid, $mnObj, $groupStrand, $mtExtremity, $mtObj) ;
		}
	}
}
sub readThatDontMapAtTheTEExtremity { # if reads do not match at the exact end of the TE
	my $self  = shift ;
	my ($teObj, $gk, $featureID, $mtTEid, $mnObj, $groupStrand, $mtExtremity, $mtObj) = @_ ;

	if($teObj->strand == 1){
		if($mtExtremity == 5){
			if($mtObj->strand == 1 and $mnObj->strand == 1){ # checked
				my $distance = $teObj->start - $mtObj->start ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == 1 and $mnObj->strand == -1){ 
				my $distance = $teObj->start - $mtObj->start ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == -1 and $mnObj->strand == -1){ # redondante with $mtObj->strand == 1 and $mnObj->strand == -1
				my $distance = $teObj->start - $mtObj->start ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == -1 and $mnObj->strand == 1){ # redondante with $mtObj->strand == 1 and $mnObj->strand == 1
				my $distance = $teObj->start - $mtObj->start ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)} = 1 ;
				}
			}
		}
		elsif($mtExtremity == 3){
			my $distance = $teObj->end - $mtObj->end ;
			if($distance > 100 or $distance < -100){ return(1) ; }
			if($mtObj->strand == 1 and $mnObj->strand == 1){
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == 1 and $mnObj->strand == -1){ ### validated ###
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end-$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end-$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end-$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == -1 and $mnObj->strand == -1){
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == -1 and $mnObj->strand == 1){
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end-$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end-$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end-$distance)} = 1 ;
				}
			}
		}
	}
	if($teObj->strand == -1){
		if($mtExtremity == 5){
			if($mtObj->strand == 1 and $mnObj->strand == 1){ # checked
				my $distance = $teObj->end - $mtObj->end ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == 1 and $mnObj->strand == -1){ # checked
				my $distance = $teObj->end - $mtObj->end ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end-$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end-$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end-$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == -1 and $mnObj->strand == -1){ # redondante with $mtObj->strand == 1 and $mnObj->strand == -1
				my $distance = $teObj->end - $mtObj->end ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->start+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == -1 and $mnObj->strand == 1){ # redondante with $mtObj->strand == 1 and $mnObj->strand == 1
				my $distance = $teObj->end - $mtObj->end ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end-$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end-$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{5}->{($mnObj->end-$distance)} = 1 ;
				}
			}
		}
		elsif($mtExtremity == 3){
			if($mtObj->strand == 1 and $mnObj->strand == 1){ # checked
				my $distance = $teObj->start - $mtObj->start ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == 1 and $mnObj->strand == -1){ # checked
				my $distance = $teObj->start - $mtObj->start ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start-$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start-$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start-$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == -1 and $mnObj->strand == -1){ # checked
				my $distance = $teObj->start - $mtObj->start ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end+$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end+$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->end+$distance)} = 1 ;
				}
			}
			elsif($mtObj->strand == -1 and $mnObj->strand == 1){ # checked
				my $distance = $teObj->start - $mtObj->start ;
				if($distance > 100 or $distance < -100){ return(1) ; }
				if(defined $self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start-$distance)}){
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start-$distance)} += 1 ;
				}
				else{
					$self->{group}->{$featureID}->{$gk}->{exactSite}->{hist}->{3}->{($mnObj->start-$distance)} = 1 ;
				}
			}
		}
	}
}
sub print {
	my $self  = shift ;

	#---- Step 1: Select overlapping neoinsertion
	my $passedGrp ;
	foreach my $gk1 (sort {
	$self->{print}->{$a}->{chr} cmp $self->{print}->{$b}->{chr} or
	$self->{print}->{$a}->{start} <=> $self->{print}->{$b}->{start} or
	$self->{print}->{$a}->{stop} <=> $self->{print}->{$b}->{stop} 
	}
	keys %{$self->{print}}){
		foreach my $gk2 (sort{
		$self->{print}->{$a}->{chr} cmp $self->{print}->{$b}->{chr} or
		$self->{print}->{$a}->{start} <=> $self->{print}->{$b}->{start} or
		$self->{print}->{$a}->{stop} <=> $self->{print}->{$b}->{stop} 
		} keys %{$self->{print}}){
			if(defined $passedGrp->{$gk2}){ next ;}
			if($gk1 == $gk2){ next ;}
			if($self->{print}->{$gk1}->{sF} ne $self->{print}->{$gk2}->{sF}){ next ;}
			if($self->{print}->{$gk1}->{chr} ne $self->{print}->{$gk2}->{chr}){ next ;}

			# gk1 equal gk2 
			if( ($self->{print}->{$gk2}->{start} == $self->{print}->{$gk1}->{start}) and ($self->{print}->{$gk1}->{stop} == $self->{print}->{$gk2}->{stop})){
				$self->selectOverlappingNeoInsertion($gk1,$gk2) ;
			}
			# gk1 include in gk2 
			elsif( ($self->{print}->{$gk2}->{start} <= $self->{print}->{$gk1}->{start}) and ($self->{print}->{$gk1}->{stop} <= $self->{print}->{$gk2}->{stop})){
				$self->selectOverlappingNeoInsertion($gk1,$gk2) ;
			}
			# gk2 include in gk1 
			elsif($self->{print}->{$gk1}->{start} <= $self->{print}->{$gk2}->{start} and $self->{print}->{$gk2}->{stop} <= $self->{print}->{$gk1}->{stop}){
				$self->selectOverlappingNeoInsertion($gk1,$gk2) ;
			}
			# overlap 3' gk1 and 5' gk2 
			elsif($self->{print}->{$gk1}->{start} < $self->{print}->{$gk2}->{start} and $self->{print}->{$gk1}->{stop} > $self->{print}->{$gk2}->{start} and $self->{print}->{$gk1}->{stop} < $self->{print}->{$gk2}->{stop}){
				$self->selectOverlappingNeoInsertion($gk1,$gk2) ;
			}
		}
		$passedGrp->{$gk1} = 1 ;
	}

	#---- Step 2: Find neoinsertion having same set of reads (at repetitive position)
	GRP1:foreach my $gk1 (keys %{$self->{print}}){
		if(defined $self->{print}->{$gk1}->{overlappingGroup}){ next ; }
		my $nbSimilarReads ;
		foreach my $gk2 (keys %{$self->{print}}){
			if($gk1 == $gk2){ next ;}
			if(defined $self->{print}->{$gk2}->{overlappingGroup}){ next ; }
			if(defined $self->{print}->{$gk2}->{repetitivePos}){ next ;}
			$nbSimilarReads->{withGrp2} = 0 ;
			$nbSimilarReads->{inGrp1} = 0 ;
			$nbSimilarReads->{grp1} = $gk1 ;
			$nbSimilarReads->{grp2} = $gk2 ;
			foreach my $extremity (sort keys %{$self->{print}->{$gk1}->{readsList}->{neo}}){
				foreach my $readsID (sort keys %{$self->{print}->{$gk1}->{readsList}->{neo}->{$extremity}}){
					$nbSimilarReads->{inGrp1}++ ;
					if(defined $self->{print}->{$gk2}->{readsList}->{neo}->{$extremity}->{$readsID}){
						$nbSimilarReads->{withGrp2}++ ;
					}
				}
			}
			if($nbSimilarReads->{withGrp2}/$nbSimilarReads->{inGrp1} >= 0.8){
				$self->{print}->{$gk1}->{repetitivePos} = 1 ;
				$self->{print}->{$gk2}->{mappability} = "multi" ;
				next GRP1 ;
			}
		}
	}


	my $dir = cwd ;
	my $fileName = $dir.'/'.$self->{option}->{unmappedFastqFile}->{prefix}.".newInsertionSite.tab" ;
	open(COORD,'>', $fileName) or die "Can't create $fileName.";

	foreach my $gk (
		sort {
		$self->{print}->{$a}->{chr} cmp $self->{print}->{$b}->{chr} or
		$self->{print}->{$a}->{start} <=> $self->{print}->{$b}->{start} or 
		$self->{print}->{$a}->{stop} <=> $self->{print}->{$b}->{stop} 
		}
		keys %{$self->{print}}){
		if(defined $self->{print}->{$gk}->{overlappingGroup}){ next ; }
		if(defined $self->{print}->{$gk}->{repetitivePos}){ next ; }
		if(not defined $self->{print}->{$gk}->{mappability}){
			$self->{print}->{$gk}->{mappability} = "uniq" ;
		}
		( $self->{print}->{$gk}->{TSDstart}, $self->{print}->{$gk}->{TSDstop} ) = ( $self->{print}->{$gk}->{TSDstop}, $self->{print}->{$gk}->{TSDstart} ) if ( $self->{print}->{$gk}->{TSDstart} > $self->{print}->{$gk}->{TSDstop} );

		my @keys = sort keys (%{$self->{print}->{$gk}->{teid}}) ;
		my $teid = join("|", @keys) ;

		$self->methylationCalc($gk, "neo") ;
		$self->methylationCalc($gk, "te") ;

		print COORD join ("\t",
		$self->{print}->{$gk}->{chr},
		$self->{print}->{$gk}->{start},
		$self->{print}->{$gk}->{stop},
		"neo".$gk,
		$self->{print}->{$gk}->{mappability},
		$self->{print}->{$gk}->{strand},
		$self->{print}->{$gk}->{TSDstart},
		$self->{print}->{$gk}->{TSDstop},
		$self->{print}->{$gk}->{nbOverallReads},
		# $self->{print}->{$gk}->{nbR5},
		# $self->{print}->{$gk}->{nbR3},
		$self->{print}->{$gk}->{family},
		$teid
		),"\n" ;

	}
}
sub printBam{
	my $self = shift ;
	my $dir = cwd ;
	my $fileName = $dir.'/'.$self->{option}->{unmappedFastqFile}->{prefix}.".newInsertionSite.sam" ;
	open(SAM,'>', $fileName) or die "Can't create $fileName.";

	# print Header
	my $header = $self->{sam}->header ;
	my @seq_ids = $self->{sam}->seq_ids ;
	print SAM join("\t", "\@HD", "VN:1.0", "SO:coordinate"),"\n" ;
	foreach my $id (@seq_ids){
		print SAM join("\t", "\@SQ", "SN:".$id, "LN:".$self->{sam}->length($id)),"\n" ;
	}
	print SAM join("\t", "\@PG", $prog),"\n" ;

	foreach my $seqid (sort keys %{$self->{printBam}->{reads}}){
		foreach my $pos (sort keys %{$self->{printBam}->{reads}->{$seqid}}){
			foreach my $readId (sort keys %{$self->{printBam}->{reads}->{$seqid}->{$pos}}){

				print SAM join("\t",
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->qname,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->flag,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->seq_id,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->start,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->qual,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->cigar_str,
					"*",
					# $self->{printBam}->{reads}->{$kr}->mtid,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->mpos,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->isize,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->query->dna,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->qstring,
					$self->{printBam}->{reads}->{$seqid}->{$pos}->{$readId}->aux
					),"\n" ;
			}
		}
	}
	close(SAM);
}
sub selectOverlappingNeoInsertion{
	my $self = shift ;
	my ($gk1, $gk2) = @_ ;

	if(defined $self->{targetFeatTmp}->{$self->{print}->{$gk1}->{featureID}}->{'te'}){
		$self->{print}->{$gk1}->{overlappingGroup} = "YES" ;
		foreach my $tmpTeid (sort keys %{$self->{print}->{$gk1}->{teid}}){
			$self->{print}->{$gk2}->{teid}->{$tmpTeid} = 1 ;
		}
	}
	elsif(defined $self->{targetFeatTmp}->{$self->{print}->{$gk2}->{featureID}}->{'te'}){
		$self->{print}->{$gk2}->{overlappingGroup} = "YES" ;
		foreach my $tmpTeid (sort keys %{$self->{print}->{$gk2}->{teid}}){
			$self->{print}->{$gk1}->{teid}->{$tmpTeid} = 1 ;
		}
	}
	elsif(($self->{print}->{$gk1}->{TSDstart} ne "nd" and $self->{print}->{$gk1}->{TSDstop} ne "nd") and
		($self->{print}->{$gk2}->{TSDstart} ne "nd" or $self->{print}->{$gk2}->{TSDstop} ne "nd")){
		$self->{print}->{$gk2}->{overlappingGroup} = "YES" ;
		foreach my $tmpTeid (sort keys %{$self->{print}->{$gk2}->{teid}}){
			$self->{print}->{$gk1}->{teid}->{$tmpTeid} = 1 ;
		}
	}
	elsif(($self->{print}->{$gk2}->{TSDstart} ne "nd" and $self->{print}->{$gk2}->{TSDstop} ne "nd") and
		($self->{print}->{$gk1}->{TSDstart} ne "nd" or $self->{print}->{$gk1}->{TSDstop} ne "nd")){
		$self->{print}->{$gk1}->{overlappingGroup} = "YES" ;
		foreach my $tmpTeid (sort keys %{$self->{print}->{$gk1}->{teid}}){
			$self->{print}->{$gk2}->{teid}->{$tmpTeid} = 1 ;
		}
	}
	elsif($self->{print}->{$gk1}->{nbOverallReads} < $self->{print}->{$gk2}->{nbOverallReads}){
		$self->{print}->{$gk1}->{overlappingGroup} = "YES" ;
		foreach my $tmpTeid (sort keys %{$self->{print}->{$gk1}->{teid}}){
			$self->{print}->{$gk2}->{teid}->{$tmpTeid} = 1 ;
		}
	}
	elsif($self->{print}->{$gk1}->{nbOverallReads} > $self->{print}->{$gk2}->{nbOverallReads}){
		$self->{print}->{$gk2}->{overlappingGroup} = "YES" ;
		foreach my $tmpTeid (sort keys %{$self->{print}->{$gk2}->{teid}}){
			$self->{print}->{$gk1}->{teid}->{$tmpTeid} = 1 ;
		}
	}
	else{
		$self->{print}->{$gk2}->{overlappingGroup} = "YES" ;
		foreach my $tmpTeid (sort keys %{$self->{print}->{$gk2}->{teid}}){
			$self->{print}->{$gk1}->{teid}->{$tmpTeid} = 1 ;
		}
	}
}
sub methylation{
	# 5' neoinsertionPosition maxReadsLength->1
	# 3' neoinsertionPosition 1->maxReadsLength

	my $self = shift ;
	my ($kPrint, $gk, $featureID) = @_ ;

	####### At the insertion site
	# assay meth at insertion site. 5'
	if($self->{print}->{$kPrint}->{strand} == 1){
		my $metStart  = ($self->{print}->{$kPrint}->{TSDstart}-$self->{option}->{maxReadsLength}) ;
		my $metEnd    = $self->{print}->{$kPrint}->{TSDstart} ;
		my $region    = $self->{print}->{$kPrint}->{chr}.":".$metStart."-".$metEnd ;
        $self->{sam}->pileup(
            $region,
            sub {
                my ($seqid, $pos, $p ) = @_;
                if($pos < $metStart or $pos > $metEnd){ return(1) ;}
				my $nucNumber = abs((($pos-$metStart)-$self->{option}->{maxReadsLength})-1) ;
				my $refBase=$self->{sam}->segment($seqid,$pos,$pos)->dna;

				if($refBase eq "C"){
					# print Dumper $seqid."\t".$pos."\t".$nucNumber."\t".$refBase ;
					my $ctx = $self->{sam}->segment($seqid,$pos+1,$pos+2)->dna;
					my @array = split("", $ctx) ;
					if($array[0] eq "G"){ # identifying the context CG, CHG, CHH
						$self->metIdentification("neo", "CG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[0] ne "G" and $array[1] eq "G"){
						$self->metIdentification("neo", "CHG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[0] ne "G" and $array[1] ne "G"){
						$self->metIdentification("neo", "CHH", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
				}
				elsif($refBase eq "G"){
					my $ctx = $self->{sam}->segment($seqid,$pos-2,$pos-1)->dna;
					my @array = split("", $ctx) ;
					if($array[1] eq "C"){ # identifying the context CG, CHG, CHH
						$self->metIdentification("neo", "CG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[1] ne "C" and $array[0] eq "C"){
						$self->metIdentification("neo", "CHG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[1] ne "C" and $array[0] ne "C"){
						$self->metIdentification("neo", "CHH", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
				}
            }
        );
	}
	elsif($self->{print}->{$kPrint}->{strand} == -1){
		my $metStart = $self->{print}->{$kPrint}->{TSDstart} ;
		my $metEnd   = ($self->{print}->{$kPrint}->{TSDstart}+$self->{option}->{maxReadsLength}) ;
		my $region   = $self->{print}->{$kPrint}->{chr}.":".$metStart."-".$metEnd ;
        $self->{sam}->pileup(
            $region,
            sub {
                my ($seqid, $pos, $p ) = @_;
                if($pos < $metStart or $pos > $metEnd){ return(1) ;}
				my $nucNumber = ($self->{option}->{maxReadsLength}-($metEnd-$pos))+1 ;
				my $refBase=$self->{sam}->segment($seqid,$pos,$pos)->dna;

				if($refBase eq "C"){
					my $ctx = $self->{sam}->segment($seqid,$pos+1,$pos+2)->dna;
					my @array = split("", $ctx) ;
					if($array[0] eq "G"){ # identifying the context CG, CHG, CHH
						$self->metIdentification("neo", "CG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[0] ne "G" and $array[1] eq "G"){
						$self->metIdentification("neo", "CHG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[0] ne "G" and $array[1] ne "G"){
						$self->metIdentification("neo", "CHH", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
				}
				elsif($refBase eq "G"){
					my $ctx = $self->{sam}->segment($seqid,$pos-2,$pos-1)->dna;
					my @array = split("", $ctx) ;
					if($array[1] eq "C"){ # identifying the context CG, CHG, CHH
						$self->metIdentification("neo", "CG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[1] ne "C" and $array[0] eq "C"){
						$self->metIdentification("neo", "CHG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[1] ne "C" and $array[0] ne "C"){
						$self->metIdentification("neo", "CHH", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
				}
            }
        );
	}
	#assay meth at insertion site. 3'
	if($self->{print}->{$kPrint}->{strand} == 1){
		my $metStart  = $self->{print}->{$kPrint}->{TSDstop} ;
		my $metEnd    = ($self->{print}->{$kPrint}->{TSDstop}+$self->{option}->{maxReadsLength}) ;
		my $region    = $self->{print}->{$kPrint}->{chr}.":".$metStart."-".$metEnd ;
        $self->{sam}->pileup(
            $region,
            sub {
                my ($seqid, $pos, $p ) = @_;
                if($pos < $metStart or $pos > $metEnd){ return(1) ;}
				my $nucNumber = ($self->{option}->{maxReadsLength}-($metEnd-$pos))+1 ;
				my $refBase=$self->{sam}->segment($seqid,$pos,$pos)->dna ;

				if($refBase eq "C"){
					my $ctx = $self->{sam}->segment($seqid,$pos+1,$pos+2)->dna;
					my @array = split("", $ctx) ;
					if($array[0] eq "G"){ # identifying the context CG, CHG, CHH
						$self->metIdentification("neo", "CG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[0] ne "G" and $array[1] eq "G"){
						$self->metIdentification("neo", "CHG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[0] ne "G" and $array[1] ne "G"){
						$self->metIdentification("neo", "CHH", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
				}
				elsif($refBase eq "G"){
					my $ctx = $self->{sam}->segment($seqid,$pos-2,$pos-1)->dna;
					my @array = split("", $ctx) ;
					if($array[1] eq "C"){ # identifying the context CG, CHG, CHH
						$self->metIdentification("neo", "CG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[1] ne "C" and $array[0] eq "C"){
						$self->metIdentification("neo", "CHG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[1] ne "C" and $array[0] ne "C"){
						$self->metIdentification("neo", "CHH", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
				}
            }
        );
	}
	elsif($self->{print}->{$kPrint}->{strand} == -1){
		my $metStart  = ($self->{print}->{$kPrint}->{TSDstop}-$self->{option}->{maxReadsLength}) ;
		my $metEnd    = $self->{print}->{$kPrint}->{TSDstop} ;
		my $region    = $self->{print}->{$kPrint}->{chr}.":".$metStart."-".$metEnd ;
        $self->{sam}->pileup(
            $region,
            sub {
                my ($seqid, $pos, $p ) = @_;
                if($pos < $metStart or $pos > $metEnd){ return(1) ;}
				my $nucNumber = abs((($pos-$metStart)-$self->{option}->{maxReadsLength})-1) ;
				my $refBase=$self->{sam}->segment($seqid,$pos,$pos)->dna;

				if($refBase eq "C"){
					my $ctx = $self->{sam}->segment($seqid,$pos+1,$pos+2)->dna;
					my @array = split("", $ctx) ;
					if($array[0] eq "G"){ # identifying the context CG, CHG, CHH
						$self->metIdentification("neo", "CG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[0] ne "G" and $array[1] eq "G"){
						$self->metIdentification("neo", "CHG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[0] ne "G" and $array[1] ne "G"){
						$self->metIdentification("neo", "CHH", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
				}
				elsif($refBase eq "G"){
					my $ctx = $self->{sam}->segment($seqid,$pos-2,$pos-1)->dna;
					my @array = split("", $ctx) ;
					if($array[1] eq "C"){ # identifying the context CG, CHG, CHH
						$self->metIdentification("neo", "CG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[1] ne "C" and $array[0] eq "C"){
						$self->metIdentification("neo", "CHG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
					elsif($array[1] ne "C" and $array[0] ne "C"){
						$self->metIdentification("neo", "CHH", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
					}
				}
            }
        );
	}

	####### At the TE
	# assay meth at TE 5'
	foreach my $teid (sort keys %{$self->{print}->{$kPrint}->{teid}}){
		if($self->{targetFeat}->{$teid}->{te}->strand == 1){
			my $metStart  = $self->{targetFeat}->{$teid}->{te}->start ;
			my $metEnd    = $self->{targetFeat}->{$teid}->{te}->start+$self->{option}->{maxReadsLength} ;
			my $region    = $self->{targetFeat}->{$teid}->{te}->seq_id.":".$metStart."-".$metEnd ;
	        $self->{sam}->pileup(
	            $region,
	            sub {
	                my ($seqid, $pos, $p ) = @_;
	                if($pos < $metStart or $pos > $metEnd){ return(1) ;}
					my $nucNumber = ($self->{option}->{maxReadsLength}-($metEnd-$pos))+1 ;
					my $refBase = $self->{sam}->segment($seqid,$pos,$pos)->dna;
					if($refBase eq "C"){
						# print Dumper $seqid."\t".$pos."\t".$nucNumber."\t".$refBase ;
						my $ctx = $self->{sam}->segment($seqid,$pos+1,$pos+2)->dna;
						my @array = split("", $ctx) ;
						if($array[0] eq "G"){ # identifying the context CG, CHG, CHH
							$self->metIdentification("te", "CG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[0] ne "G" and $array[1] eq "G"){
							$self->metIdentification("te", "CHG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[0] ne "G" and $array[1] ne "G"){
							$self->metIdentification("te", "CHH", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
					}
					elsif($refBase eq "G"){
						my $ctx = $self->{sam}->segment($seqid,$pos-2,$pos-1)->dna;
						my @array = split("", $ctx) ;
						if($array[1] eq "C"){ # identifying the context CG, CHG, CHH
							$self->metIdentification("te", "CG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[1] ne "C" and $array[0] eq "C"){
							$self->metIdentification("te", "CHG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[1] ne "C" and $array[0] ne "C"){
							$self->metIdentification("te", "CHH", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
					}
	            }
	        );
		}
		elsif($self->{targetFeat}->{$teid}->{te}->strand == -1){
			my $metStart = $self->{targetFeat}->{$teid}->{te}->end-$self->{option}->{maxReadsLength} ;
			my $metEnd   = $self->{targetFeat}->{$teid}->{te}->end ;
			my $region   = $self->{targetFeat}->{$teid}->{te}->seq_id.":".$metStart."-".$metEnd ;
	        $self->{sam}->pileup(
	            $region,
	            sub {
	                my ($seqid, $pos, $p ) = @_;
	                if($pos < $metStart or $pos > $metEnd){ return(1) ;}
					my $nucNumber = abs((($pos-$metStart)-$self->{option}->{maxReadsLength})-1) ;
					my $refBase=$self->{sam}->segment($seqid,$pos,$pos)->dna;

					if($refBase eq "C"){
						my $ctx = $self->{sam}->segment($seqid,$pos+1,$pos+2)->dna;
						my @array = split("", $ctx) ;
						if($array[0] eq "G"){ # identifying the context CG, CHG, CHH
							$self->metIdentification("te", "CG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[0] ne "G" and $array[1] eq "G"){
							$self->metIdentification("te", "CHG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[0] ne "G" and $array[1] ne "G"){
							$self->metIdentification("te", "CHH", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
					}
					elsif($refBase eq "G"){
						my $ctx = $self->{sam}->segment($seqid,$pos-2,$pos-1)->dna;
						my @array = split("", $ctx) ;
						if($array[1] eq "C"){ # identifying the context CG, CHG, CHH
							$self->metIdentification("te", "CG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[1] ne "C" and $array[0] eq "C"){
							$self->metIdentification("te", "CHG", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[1] ne "C" and $array[0] ne "C"){
							$self->metIdentification("te", "CHH", "5", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
					}
	            }
	        );
		}
		#assay meth at insertion site. 3'
		if($self->{targetFeat}->{$teid}->{te}->strand == 1){
			my $metStart = $self->{targetFeat}->{$teid}->{te}->end-$self->{option}->{maxReadsLength} ;
			my $metEnd   = $self->{targetFeat}->{$teid}->{te}->end ;
			my $region   = $self->{targetFeat}->{$teid}->{te}->seq_id.":".$metStart."-".$metEnd ;
	        $self->{sam}->pileup(
	            $region,
	            sub {
	                my ($seqid, $pos, $p ) = @_;
	                if($pos < $metStart or $pos > $metEnd){ return(1) ;}
					my $nucNumber = abs((($pos-$metStart)-$self->{option}->{maxReadsLength})-1) ;
					my $refBase=$self->{sam}->segment($seqid,$pos,$pos)->dna;

					if($refBase eq "C"){
						my $ctx = $self->{sam}->segment($seqid,$pos+1,$pos+2)->dna;
						my @array = split("", $ctx) ;
						if($array[0] eq "G"){ # identifying the context CG, CHG, CHH
							$self->metIdentification("te", "CG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[0] ne "G" and $array[1] eq "G"){
							$self->metIdentification("te", "CHG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[0] ne "G" and $array[1] ne "G"){
							$self->metIdentification("te", "CHH", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
					}
					elsif($refBase eq "G"){
						my $ctx = $self->{sam}->segment($seqid,$pos-2,$pos-1)->dna;
						my @array = split("", $ctx) ;
						if($array[1] eq "C"){ # identifying the context CG, CHG, CHH
							$self->metIdentification("te", "CG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[1] ne "C" and $array[0] eq "C"){
							$self->metIdentification("te", "CHG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[1] ne "C" and $array[0] ne "C"){
							$self->metIdentification("te", "CHH", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
					}
	            }
	        );
		}
		elsif($self->{targetFeat}->{$teid}->{te}->strand == -1){
			my $metStart = $self->{targetFeat}->{$teid}->{te}->start ;
			my $metEnd   = $self->{targetFeat}->{$teid}->{te}->start+$self->{option}->{maxReadsLength} ;
			my $region   = $self->{targetFeat}->{$teid}->{te}->seq_id.":".$metStart."-".$metEnd ;
	        $self->{sam}->pileup(
	            $region,
	            sub {
	                my ($seqid, $pos, $p ) = @_;
	                if($pos < $metStart or $pos > $metEnd){ return(1) ;}
					my $nucNumber = ($self->{option}->{maxReadsLength}-($metEnd-$pos))+1 ;
					my $refBase=$self->{sam}->segment($seqid,$pos,$pos)->dna;

					if($refBase eq "C"){
						my $ctx = $self->{sam}->segment($seqid,$pos+1,$pos+2)->dna;
						my @array = split("", $ctx) ;
						if($array[0] eq "G"){ # identifying the context CG, CHG, CHH
							$self->metIdentification("te", "CG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[0] ne "G" and $array[1] eq "G"){
							$self->metIdentification("te", "CHG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[0] ne "G" and $array[1] ne "G"){
							$self->metIdentification("te", "CHH", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
					}
					elsif($refBase eq "G"){
						my $ctx = $self->{sam}->segment($seqid,$pos-2,$pos-1)->dna;
						my @array = split("", $ctx) ;
						if($array[1] eq "C"){ # identifying the context CG, CHG, CHH
							$self->metIdentification("te", "CG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[1] ne "C" and $array[0] eq "C"){
							$self->metIdentification("te", "CHG", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
						elsif($array[1] ne "C" and $array[0] ne "C"){
							$self->metIdentification("te", "CHH", "3", $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) ;
						}
					}
	            }
	        );
		}
	} # end of foreach teid
}
sub metIdentification {
	my $self = shift ;
	my ($area, $ctx, $mtExtremity, $nucNumber, $refBase, $p, $gk, $featureID, $pos, $kPrint) = @_ ; 
	for my $pileup(@{$p}){
		my $alignment = $pileup->alignment ;
		if(not defined $self->{group}->{$featureID}->{$gk}->{readsList}->{$area}->{$mtExtremity}->{$alignment->display_name}){ next ; } # read is not part of the group
		my $b = $pileup->b ;
		if($refBase eq "C" and $alignment->strand == 1){
			$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{ctx} = $ctx ;
			$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{pos} = $pos ;
			my $readBase=substr($b->qseq,$pileup->qpos,1);
			if($refBase eq "C" and $readBase eq "C"){ # cytosine is methylated
				$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{nbCm} += 1 ;
				$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{nbR} += 1 ;
			}
			elsif($refBase eq "C" and $readBase eq "T"){
				$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{nbR} += 1 ;
			}
		} 
		elsif($refBase eq "G" and $alignment->strand == -1){ # only consider reads in the - strand
			$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{ctx} = $ctx ;
			$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{pos} = $pos ;
			my $readBase=substr($b->qseq,$pileup->qpos,1);
			if($refBase eq "G" and $readBase eq "G"){ # cytosine is methylated
				$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{nbCm} += 1 ;
				$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{nbR} += 1 ;
			}
			elsif($refBase eq "G" and $readBase eq "A"){
				$self->{met}->{$area}->{$kPrint}->{$mtExtremity}->{$nucNumber}->{nbR} += 1 ;
			}
		}
	} # end of for pileup
}
sub methylationCalc {
	my $self = shift ;
	my ($gk, $area) = @_ ;

	my $teid ;
	foreach my $key (sort keys %{$self->{print}->{$gk}->{teid}}){
		$teid .= $key."|" ;
	}
	chop($teid);

	for(my $i = 1 ; $i<=$self->{option}->{maxReadsLength} ; $i++){
		if(defined $self->{met}->{$area}->{$gk}->{5}->{$i}){
			my $ctx = $self->{met}->{$area}->{$gk}->{5}->{$i}->{ctx} ;
			if(not defined $self->{met}->{$area}->{$gk}->{5}->{$i}->{nbR}){ next ;}
			if(not defined $self->{met}->{$area}->{$gk}->{5}->{$i}->{nbCm}){
				$self->{met}->{$area}->{$gk}->{5}->{$i}->{nbCm} = 0 ;
			}
            push(@{ $self->{windows}->{$ctx}->{ $self->{nc}->{$i} }->{$area}->{met}->{"5"}->{$gk}->{nbCm} }, $self->{met}->{$area}->{$gk}->{5}->{$i}->{nbCm});
            push(@{ $self->{windows}->{$ctx}->{ $self->{nc}->{$i} }->{$area}->{met}->{"5"}->{$gk}->{nbR} },  $self->{met}->{$area}->{$gk}->{5}->{$i}->{nbR});
            push(@{ $self->{windows}->{$ctx}->{ $self->{nc}->{$i} }->{$area}->{met}->{"8"}->{$gk}->{nbCm} }, $self->{met}->{$area}->{$gk}->{5}->{$i}->{nbCm});
            push(@{ $self->{windows}->{$ctx}->{ $self->{nc}->{$i} }->{$area}->{met}->{"8"}->{$gk}->{nbR} },  $self->{met}->{$area}->{$gk}->{5}->{$i}->{nbR});

			#### save row data
			my $tmp = join("\t", $self->{met}->{$area}->{$gk}->{5}->{$i}->{ctx}, $area, "5", $self->{met}->{$area}->{$gk}->{5}->{$i}->{nbCm}, $self->{met}->{$area}->{$gk}->{5}->{$i}->{nbR}, "neo".$gk, $teid)."\n" ;
			push(@{$self->{printRowDataMet}}, $tmp) ;
		}
		elsif(defined $self->{met}->{$area}->{$gk}->{3}->{$i}){
			my $ctx = $self->{met}->{$area}->{$gk}->{3}->{$i}->{ctx} ;
			if(not defined $self->{met}->{$area}->{$gk}->{3}->{$i}->{nbR}){ next ;}
			if(not defined $self->{met}->{$area}->{$gk}->{3}->{$i}->{nbCm}){
				$self->{met}->{$area}->{$gk}->{3}->{$i}->{nbCm} = 0 ;
			}
            push(@{ $self->{windows}->{$ctx}->{ $self->{nc}->{$i} }->{$area}->{met}->{"3"}->{$gk}->{nbCm} }, $self->{met}->{$area}->{$gk}->{3}->{$i}->{nbCm});
            push(@{ $self->{windows}->{$ctx}->{ $self->{nc}->{$i} }->{$area}->{met}->{"3"}->{$gk}->{nbR} },  $self->{met}->{$area}->{$gk}->{3}->{$i}->{nbR});
            push(@{ $self->{windows}->{$ctx}->{ $self->{nc}->{$i} }->{$area}->{met}->{"8"}->{$gk}->{nbCm} }, $self->{met}->{$area}->{$gk}->{3}->{$i}->{nbCm});
            push(@{ $self->{windows}->{$ctx}->{ $self->{nc}->{$i} }->{$area}->{met}->{"8"}->{$gk}->{nbR} },  $self->{met}->{$area}->{$gk}->{3}->{$i}->{nbR});

			#### save row data
			my $tmp = join("\t", $self->{met}->{$area}->{$gk}->{3}->{$i}->{ctx}, $area, "3", $self->{met}->{$area}->{$gk}->{3}->{$i}->{nbCm}, $self->{met}->{$area}->{$gk}->{3}->{$i}->{nbR}, "neo".$gk, $teid)."\n" ;
			push(@{$self->{printRowDataMet}}, $tmp) ;
		}
	}
}
sub makeWindows {
    my $self = shift;
    my $k    = 0;
    my $nc   = 0;
	my $remainder = $self->{option}->{maxReadsLength} % $windowSize;
	my $quotient = ($self->{option}->{maxReadsLength} - $remainder) / $windowSize;
	if($remainder != 0){ $quotient = $quotient+1 ; }

    for (my $i = 0 ; $i < $quotient; $i++)
    {	
    	++$k ;
        for ( my $j = 1 ; $j <= $windowSize ; $j++ ) {
            $self->{nc}->{ ++$nc } = $k;
        }
    }
}
sub calcMetPercentage {
    my $self = shift;
    foreach my $ctx (sort keys %{ $self->{windows} } ) {
	    foreach my $windK (sort { $a <=> $b } keys %{ $self->{windows}->{$ctx} } ) {
	   	    foreach my $area (sort keys %{$self->{windows}->{$ctx}->{$windK}} ){
		        if ( not defined $self->{windows}->{$ctx}->{$windK}->{$area}->{met} ) {
		            print STDERR "WARNING $prog: Not Cytosine in the window: ".$self->{windows}->{$ctx}->{$windK}->{$area}->{start}."-".$self->{windows}->{$ctx}->{$windK}->{$area}->{stop}."\n" ;
		            next ;
		        }
		        foreach my $extremity (sort keys %{$self->{windows}->{$ctx}->{$windK}->{$area}->{met}}){
			        foreach my $gk (sort keys %{$self->{windows}->{$ctx}->{$windK}->{$area}->{met}->{$extremity}}){
			            my $statNbCm = Statistics::Descriptive::Full->new();
			            $statNbCm->add_data(@{$self->{windows}->{$ctx}->{$windK}->{$area}->{met}->{$extremity}->{$gk}->{nbCm}});
			            my $statNbR = Statistics::Descriptive::Full->new();
			            $statNbR->add_data(@{$self->{windows}->{$ctx}->{$windK}->{$area}->{met}->{$extremity}->{$gk}->{nbR}});
			            my $percMet = ( $statNbCm->sum / $statNbR->sum ) * 100;
			            push(@{ $self->{windows}->{$ctx}->{$windK}->{$area}->{cov}->{$extremity} }, $percMet);
			        }
		        }
		    }
	    }
	}
 }
sub calcMeanAndCi {
    my $self = shift;
    my $displayType = "smooth" ;
    foreach my $ctx ( sort keys %{ $self->{windows} } ) {
	    foreach my $windK ( sort { $a <=> $b } keys %{ $self->{windows}->{$ctx} } ) {
			foreach my $area(sort keys %{$self->{windows}->{$ctx}->{$windK}} ){
		        foreach my $extremity ( sort keys %{ $self->{windows}->{$ctx}->{$windK}->{$area}->{cov} } )
		        {
		            if($displayMet eq "metaplot"){
			            my $statObj = Statistics::Descriptive::Full->new();
			            $statObj->add_data( @{ $self->{windows}->{$ctx}->{$windK}->{$area}->{cov}->{$extremity}} );
			            my $count  = sprintf "%.8f", $statObj->count;
			            my $sum    = sprintf "%.8f", $statObj->sum;
			            my $median = sprintf "%.8f", $statObj->median;
			            my $stddev = sprintf "%.8f", $statObj->standard_deviation;
			            my $min    = sprintf "%.8f", $statObj->min;
			            my $max    = sprintf "%.8f", $statObj->max;
			            $self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{mean} = sprintf "%.8f", $statObj->mean;
			            $self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{ci} = sprintf "%.8f", ( 1.960 * ( $stddev / sqrt($count) ) );
			            $self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{count} = sprintf "%.0f", $count ;
		            }
  		            elsif($displayMet eq "smooth"){
						@{$self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{smooth}} = @{ $self->{windows}->{$ctx}->{$windK}->{$area}->{cov}->{$extremity} } ;
  		            }
		        }
		    }
		}
	}
}
sub printMetaplot{
	my $self = shift ;
	my $dir = cwd ;
	my $fileName = $dir.'/'.$self->{option}->{unmappedFastqFile}->{prefix}.".met.meta.tab" ;
	open(META,'>', $fileName) or die "Can't create $fileName.";

	my $print ;
	# reorganized
	foreach my $ctx(sort keys %{$self->{printMet}} ){
		foreach my $area(sort keys %{$self->{printMet}->{$ctx}} ){
			foreach my $extremity(sort {$b<=>$a} keys %{$self->{printMet}->{$ctx}->{$area}}){
				foreach my $windK(sort {$a<=>$b} keys %{$self->{printMet}->{$ctx}->{$area}->{$extremity}}){
					my $bin = $windK ;
					if($extremity == 5 and $area eq "neo"){
						$bin = (-1*$windK)+1 ;
					}
					elsif($extremity == 5 and $area eq "te"){
						$bin = $windK  ;
					}
					elsif($extremity == 3 and $area eq "te"){
						$bin = ($windK-scalar(keys %{$self->{printMet}->{$ctx}->{$area}->{$extremity}}))-1 ;
					}
					elsif($extremity == 3 and $area eq "neo"){
						$bin = $windK - 1 ;
					}
					elsif($extremity == 8 and $area eq "neo"){
						$bin = (-1*$windK)+1 ;
					}
					elsif($extremity == 8 and $area eq "te"){
						$bin = $windK ;
					}
					if(not defined $self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{smooth}){
						$print->{$ctx}->{$area}->{$extremity}->{$bin}->{metaplot} = join ("\t",
							$ctx,
							$area,
							$extremity,
							$bin,
							$self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{mean},
							$self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{ci}
						) ;
					}
					elsif(defined $self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{smooth}){
						my $k = 0 ;
						foreach (@{$self->{printMet}->{$ctx}->{$area}->{$extremity}->{$windK}->{smooth}}){
							$print->{$ctx}->{$area}->{$extremity}->{$bin}->{smooth}->{++$k} = join ("\t",
								$ctx,
								$area,
								$extremity,
								$bin,
								$_
							) ;
						}
					}
				}
			}
		}
	}
	foreach my $ctx (sort keys %{$print}){
		foreach my $area(sort keys %{$print->{$ctx}}){
			foreach my $extremity(sort {$b<=>$a} keys %{$print->{$ctx}->{$area}}){
				foreach my $windK(sort {$a<=>$b} keys %{$print->{$ctx}->{$area}->{$extremity}}){

					if(defined $print->{$ctx}->{$area}->{$extremity}->{$windK}->{metaplot}){
						print META $print->{$ctx}->{$area}->{$extremity}->{$windK}->{metaplot}."\n" ;
					}
					elsif(defined $print->{$ctx}->{$area}->{$extremity}->{$windK}->{smooth}){
						foreach my $k (sort {$a <=> $b} keys %{$print->{$ctx}->{$area}->{$extremity}->{$windK}->{smooth}}){
							print META $print->{$ctx}->{$area}->{$extremity}->{$windK}->{smooth}->{$k}."\n" ;
						}
					}
				}
			}
		}
	}
}
sub printRowDataMet{
	my $self = shift ;
	my $dir = cwd ;
	my $fileName = $dir.'/'.$self->{option}->{unmappedFastqFile}->{prefix}.".met.row.tab" ;
	open(TAB,'>', $fileName) or die "Can't create $fileName.";
	foreach(@{$self->{printRowDataMet}}){
		print TAB $_ ;
	}
}
sub sF_TSDsize {
	my $self = shift ;
	$self->{tsd}->{"LTR/Copia"} = 5 ;
	$self->{tsd}->{"LTR/Gypsy"} = 5 ;
	$self->{tsd}->{"LINE/L1"} = 5 ;
	$self->{tsd}->{"SINE"} = 5 ;
	$self->{tsd}->{"DNA/En-Spm"} = 3 ;
	$self->{tsd}->{"DNA/MuDR"} = 9 ;
	$self->{tsd}->{"RC/Helitron"} = 0 ;

# DNA
# DNA/HAT
# DNA/Harbinger
# DNA/Mariner
# DNA/Pogo
# DNA/Tc1
# RathE1_cons
# RathE2_cons
# RathE3_cons
# Unassigned
}

sub help {
my $prog = basename($0) ;
print STDERR <<EOF ;
### HELP: $prog ###
AUTHOR:     Josquin DARON, Slotkin Lab, Ohio State University
VERSION:    $VERSION -- $lastmodif

PURPOSE: Identify non-reference TE insertion sites and their methylation level.

USAGE: $prog —l [max read length] -gff <gff> -t <target> —ref <fasta> -un <fastq>

        <gff>     TE annotation should be given in gff3 format.
                  For TE annotated features, column 9 should have following list of tags:
                  ID (teid), sF (superfamily name), fam (family name).
                  For LTR annotated features, they should be referred as LTR5 or LTR3 in column 3,
                  column 9 should have tags Parent (teid).
            
        <target>  list of TEid of interest

        <fasta>   FASTA formated (.fa, .fna or fasta) genome file.

        <fastq>   FASTQ file of reads that failed to map the reference genome (unmapped reads)

    OPTIONS
      EpiTEome Specific Options:
        -chop [integer] : read ends length of chopped (defaut 25,30,40).
                          Use of several length will improve epiTEome sensitivity. 
        -b    [integer] : number of TE per batch (defaut 5000).
        -w    [integer] : window size for methylation metaplot analysis.

      Alignment Options:
        -E    [integer] : segemehl max evalue (default:5)
        -p    [integer] : number of threads use in segemehl (defaut 1).
                          All other portions are single-threaded.

    OUTPUT
      epiTEome output 4 different files such as .newInsertionSite.tab, .newInsertionSite.sam, .met.meta.tab and .met.row.tab


EOF
exit(1) ;
}
