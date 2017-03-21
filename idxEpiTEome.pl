#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename ;
use Bio::Tools::GFF;
use Data::Dumper ;
use File::Which;
use Cwd;

my $VERSION = "1.0" ;
my $lastmodif = "2012-01-12" ;
my $prog = basename($0) ;

my $gffFile ;
my $target ;
my $readLength ;
my $refFasta ;
my $help ;

&GetOptions(
	"gff:s"     => \$gffFile,
	"t:s"       => \$target,
	"fasta:s"   => \$refFasta,
	"l:s"       => \$readLength,
	"h"         => \$help
);

my $self = {};
bless $self;
$self->checkPrograms();
$self->setOptions();
$self->readTargetFile();
$self->readGffFile();
$self->maskFastaIndex();

$help and &help;

# &main( \@ARGV );

sub setOptions {
	my $self = shift;
	if (defined ($gffFile)){
		unless (-e $gffFile){ print STDERR "DIE $prog: Could not find file: " . $gffFile . "\n" and die ; }
		$self->{option}->{gff} = $gffFile ;
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -gff.\n" and &help ;
	}
	if (defined ($target)){
		unless (-e $target){ print STDERR "DIE $prog: Could not find file: " . $target . "\n" and die ; }
		$self->{option}->{targetList} = $target ;
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -t.\n" and &help ;
	}
	if(defined $readLength){
		$self->{option}->{maxReadsLength} = $readLength ;
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -l.\n" and &help ;
	}
	if (defined ($refFasta)){
		unless (-e $refFasta){ print STDERR "DIE $prog: Could not find file: " . $refFasta . "\n" and die ; }
		if($refFasta =~/(.*).fasta/ or $refFasta =~/(.*).fna/ or $refFasta =~/(.*).fa/){
			$self->{fileName}->{refFasta}->{file} = $refFasta ;
			$self->{fileName}->{refFasta}->{prefix} = $1 ;
		}
		else{
			print STDERR "ERROR $prog: File give in -fasta option should have extention fasta, fa or fna.\n" and &help ;
		}
	}
	else{
		print STDERR "ERROR $prog: Need to specify option -fasta.\n" and &help ;
	}
}
sub checkPrograms{
	my $self = shift;
	$self->{program}->{rm} = which 'rm';
	if(not defined $self->{program}->{rm}){ print STDERR "DIE $prog: Missing program: rm.\n" and die ; }
	$self->{program}->{maskFastaFromBed} = which 'maskFastaFromBed';
	if(not defined $self->{program}->{maskFastaFromBed}){ print STDERR "DIE $prog: Missing program: maskFastaFromBed.\n" and die ; }
	$self->{program}->{segemehl} = which 'segemehl.x';
	if(not defined $self->{program}->{segemehl}){ print STDERR "DIE $prog: Missing program: segemehl.x.\n" and die ; }
	print STDERR "INFO $prog: All programs have been found succesfully.\n" ;
}
sub readTargetFile {
	my $self = shift;
	open("T", $self->{option}->{targetList}) or print STDERR "Die: Could not open file: " . $self->{option}->{targetList} . "\n" and exit ;
	while(<T>){
		chomp $_ ;
		$self->{targetList}->{$_} =  1 ;
	}
}
sub readGffFile {
	my $datestring = localtime();
	print STDERR "INFO $prog: Run Module: $datestring readGffFile !\n" ;
	my $self = shift;
	my $gffio = Bio::Tools::GFF->new(-file => $self->{option}->{gff}, -gff_version => 3);
    my $feature;
    my $checkConcodanceGffTarget = "FALSE" ;

    # loop over the input stream
    while($feature = $gffio->next_feature()) {
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

    		if(defined $self->{targetList}->{$id[0]}){
    			if($feature->end-$feature->start+1>$self->{option}->{maxReadsLength}){
	    			$checkConcodanceGffTarget = "TRUE" ;
	    			$self->{targetFeatTmp}->{$id[0]}->{$feature->primary_tag} = $feature ;
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

	print STDERR "INFO $prog: Reading GFF file ... done !\n" ;
	if($checkConcodanceGffTarget eq "FALSE"){
		print STDERR "DIE $prog: Features ID in target file (-t) are not present in gff file. Die.\n" and die ;
	}
	delete($self->{targetList}) ;
	return(1) ;
}
sub maskFastaIndex{
	my $datestring = localtime();
	print STDERR "INFO $prog: Run Module: $datestring maskFastaIndex !\n" ;
	my $self = shift;
	my $dir = cwd ;
	my $bedFileOut = $dir."/prepRefSeq.bed" ;
    open(BED, '>', $bedFileOut) or die("Could not open file $bedFileOut.\n");
	foreach my $feature(sort keys %{$self->{targetFeatTmp}}){
		if(defined $self->{targetFeatTmp}->{$feature}->{'LTR5'}){
			if($self->{targetFeatTmp}->{$feature}->{'LTR5'}->length > 300){
				# + strand
				if($self->{targetFeatTmp}->{$feature}->{'te'}->strand == "1"){
					print BED $self->{targetFeatTmp}->{$feature}->{'te'}->seq_id."\t".($self->{targetFeatTmp}->{$feature}->{'LTR5'}->end-150)."\t".$self->{targetFeatTmp}->{$feature}->{'LTR5'}->end."\t".$feature."_LTR5\n" ;
				}
				# - strand
				if($self->{targetFeatTmp}->{$feature}->{'te'}->strand == "-1"){
					print BED $self->{targetFeatTmp}->{$feature}->{'te'}->seq_id."\t".$self->{targetFeatTmp}->{$feature}->{'LTR5'}->start."\t".($self->{targetFeatTmp}->{$feature}->{'LTR5'}->start+150)."\t".$feature."_LTR5\n" ;
				}
			}
			elsif($self->{targetFeatTmp}->{$feature}->{'LTR5'}->length > 100){
				# + strand
				if($self->{targetFeatTmp}->{$feature}->{'te'}->strand == "1"){
					print BED $self->{targetFeatTmp}->{$feature}->{'te'}->seq_id."\t".($self->{targetFeatTmp}->{$feature}->{'LTR5'}->end-50)."\t".$self->{targetFeatTmp}->{$feature}->{'LTR5'}->end."\t".$feature."_LTR5\n" ;
				}
				# - strand
				if($self->{targetFeatTmp}->{$feature}->{'te'}->strand == "-1"){
					print BED $self->{targetFeatTmp}->{$feature}->{'te'}->seq_id."\t".$self->{targetFeatTmp}->{$feature}->{'LTR5'}->start."\t".($self->{targetFeatTmp}->{$feature}->{'LTR5'}->start+50)."\t".$feature."_LTR5\n" ;
				}
			}
		}
		if(defined $self->{targetFeatTmp}->{$feature}->{'LTR3'}){
			if($self->{targetFeatTmp}->{$feature}->{'LTR3'}->length > 300){
				# + strand
				if($self->{targetFeatTmp}->{$feature}->{'te'}->strand == "1"){
					print BED $self->{targetFeatTmp}->{$feature}->{'te'}->seq_id."\t".$self->{targetFeatTmp}->{$feature}->{'LTR3'}->start."\t".($self->{targetFeatTmp}->{$feature}->{'LTR3'}->start+150)."\t".$feature."_LTR3\n" ;
				}
				# - strand
				if($self->{targetFeatTmp}->{$feature}->{'te'}->strand == "-1"){
					print BED $self->{targetFeatTmp}->{$feature}->{'te'}->seq_id."\t".($self->{targetFeatTmp}->{$feature}->{'LTR3'}->end-150)."\t".$self->{targetFeatTmp}->{$feature}->{'LTR3'}->end."\t".$feature."_LTR3\n" ;
				}
			}
			elsif($self->{targetFeatTmp}->{$feature}->{'LTR3'}->length > 100){
				# + strand
				if($self->{targetFeatTmp}->{$feature}->{'te'}->strand == "1"){
					print BED $self->{targetFeatTmp}->{$feature}->{'te'}->seq_id."\t".$self->{targetFeatTmp}->{$feature}->{'LTR3'}->start."\t".($self->{targetFeatTmp}->{$feature}->{'LTR3'}->start+50)."\t".$feature."_LTR3\n" ;
				}
				# - strand
				if($self->{targetFeatTmp}->{$feature}->{'te'}->strand == "-1"){
					print BED $self->{targetFeatTmp}->{$feature}->{'te'}->seq_id."\t".($self->{targetFeatTmp}->{$feature}->{'LTR3'}->end-50)."\t".$self->{targetFeatTmp}->{$feature}->{'LTR3'}->end."\t".$feature."_LTR3\n" ;
				}
			}
		}
	}
 	close BED ;

	my @maskFastaFromBed = ($self->{program}->{maskFastaFromBed}, '-fi ', $self->{fileName}->{refFasta}->{file}, '-bed', $bedFileOut, '-fo ', $self->{fileName}->{refFasta}->{prefix}.'.epiTEome.masked.fasta' ) ;
	print STDERR join ( " ", @maskFastaFromBed),"\n" ;
	system(join" ", @maskFastaFromBed) ;
	my @segIndex = ($self->{program}->{segemehl}, '--silent -x', $self->{fileName}->{refFasta}->{prefix}.'.epiTEome.masked.fasta.ctidx', '-y', $self->{fileName}->{refFasta}->{prefix}.'.epiTEome.masked.fasta.gaidx', '-d', $self->{fileName}->{refFasta}->{prefix}.'.epiTEome.masked.fasta', '-F 1 2> log') ;
	print STDERR join ( " ", @segIndex),"\n" ;
	system(join" ", @segIndex) ;

	my @rm = ("rm", $bedFileOut);
	system(join" ", @rm) ;

	print "All files have been created succesfully !\n" ;
}



sub help {
my $prog = basename($0);
print STDERR <<EOF ;
### $prog ###
#
# CREATED:    xxxxxxxxxxxxxxx
# LAST MODIF: $lastmodif
# AUTHOR:     Josquin DARON
# VERSION:    $VERSION
# PURPOSE:    
#             
#

USAGE:
       $prog  [OPTIONS]  FastaFile

       ### OPTIONS ###
       -gff file
       -t file
       -fasta reffile
       -l lenght
       -h|--help:   print this help

EOF

exit(1) ;
}

# END OF SCRIPT

