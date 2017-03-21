#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO ;
use Getopt::Long;
use File::Basename ;
use Data::Dumper ;


my $VERSION = "1.0" ;
my $lastmodif = "2017-03-02" ;
my $prog = basename($0) ;

my $outformat= "FASTA"; 
my $bedFile ;
my $fastadb ;
my $help ;

&GetOptions(
	"bed=s"   => \$bedFile,
	"db=s"    => \$fastadb,	
	"h"       => \$help
);

$help and &help;

&main( \@ARGV );

sub main {
        my $self = {};
        bless $self;
        my $datestring = localtime();
        print STDERR "INFO $prog: Running script $datestring\n" ;
        $self->tsdSize();
        $self->setOptions();
        $self->{inputFiles} = shift;
    foreach ( @{ $self->{inputFiles} } ) {
        $self->{_currentFile} = $_;
        unless ( -e $self->{_currentFile} ) {
            print STDERR "INFO $prog: Could not find file: " . $self->{_currentFile} . " -> next file\n";
            next;
        }
		my $finalSeq ;
		my $idFile ;
		my $description ;
		my $seqio  = Bio::SeqIO->new( -format => "FASTA" , -file => $self->{_currentFile});
		while (my $seq = $seqio->next_seq) {
			$idFile = $seq->display_id ;
			$description = $seq->desc ;

			my $kseq = 0;
			foreach my $chr (sort keys(%{$self->{bedFile}})){
				foreach my $bedkey (sort {$a <=> $b} keys(%{$self->{bedFile}->{$chr}})){
					my $tsd = $self->{tsd}->{$self->{bedFile}->{$chr}->{$bedkey}->{sF}} ;
					if($bedkey == 1){
						$self->{subSeq}->{++$kseq} = $seq->subseq(1,$self->{bedFile}->{$chr}->{$bedkey}->{start}) ;
						$self->{subSeq}->{++$kseq} = $self->{fastadb}->{$self->{bedFile}->{$chr}->{$bedkey}->{id}}->seq ;
					}
					else{
						my $prevSt = $self->{bedFile}->{$chr}->{$bedkey-1}->{start}-$tsd+1 ;

						if($prevSt > $self->{bedFile}->{$chr}->{$bedkey}->{start}){
				            print STDERR "DIE $prog: Bed file should be sorted\n";
							exit ;
						}
						$self->{subSeq}->{++$kseq} = $seq->subseq($prevSt,$self->{bedFile}->{$chr}->{$bedkey}->{start}) ;
						$self->{subSeq}->{++$kseq} = $self->{fastadb}->{$self->{bedFile}->{$chr}->{$bedkey}->{id}}->seq ;
					}
					if(not defined $self->{bedFile}->{$chr}->{$bedkey+1}){ # last feat, add the end of the seq
						$self->{subSeq}->{++$kseq} = $seq->subseq($self->{bedFile}->{$chr}->{$bedkey}->{start}-$tsd+1,$seq->length) ;
					}					
				}
			}
		}
		foreach my $key (sort {$a <=> $b} keys %{$self->{subSeq}}){
			$finalSeq .= $self->{subSeq}->{$key} ;
		}

		my $finalSeqObj = Bio::Seq ->new(
			                        '-seq'   => $finalSeq, 
									'-id'    => $idFile, 
									'-desc'  => $description) ;
		my $out_seqio_obj = Bio::SeqIO ->new('-fh' => \*STDOUT, '-format' => "FASTA") ;
		$out_seqio_obj ->write_seq ($finalSeqObj) ;

    } # END OF foreach inputFiles
}
sub setOptions {
    my $self = shift;
    if(defined $bedFile){
		$self->readBedFile() ;
    }
    else{
		print STDERR "ERROR $prog: Need to specify either option -bed.\n" and &help ;
    }
    if(defined $fastadb){
		$self->readFastadb() ;
    }
    else{
		print STDERR "ERROR $prog: Need to specify either option -db.\n" and &help ;
    }
}

sub readFastadb{
    my $self  = shift;
    unless (-e $fastadb) {print STDERR "ERROR $prog: Could not find file: $fastadb \n" ; &help ; }
	my $seqio  = Bio::SeqIO->new( -format => "FASTA" , -file => $fastadb);
	while (my $seq = $seqio->next_seq) {
		$self->{fastadb}->{$seq->display_id} = $seq ;
	}
}
sub readBedFile {
    my $self  = shift;
    my $count = 0;
    unless ( -e $bedFile ) {print STDERR "ERROR $prog: Could not find file: " . $bedFile . "\n" ; &help ; }

    open( "BED", $bedFile )
      or print STDERR "DIE $prog: Could not open file: " . $bedFile . "\n" and die;
    while (<BED>) {
        next if /^#/;
        chomp;
        my @line = split("\t");
        my ($chr, $start, $stop, $featid);
        if ( scalar(@line) < 3 ) {
            print STDERR "ERROR $prog: Bed file must have 3 column in line: $_ \n";
            die;
        }
        elsif(scalar(@line) == 3 ){
                ($chr, $start, $stop) = @line;
                $featid = $chr.":".$start."-".$stop ;
        }
        else{
                        ($chr, $start, $stop, $featid) = @line;
        }


        if ( $start !~ /^\d+$/ or $stop !~ /^\d+$/ ) {
            print STDERR "ERROR $prog: Warning in bed file: start or stop should be a number in line: $_ \n";
            next;
        }

        ( $start, $stop ) = ( $stop, $start ) if ( $start > $stop );
        $start += 1;    # put start in base 1
        $count += 1 ;

        $self->{bedFile}->{$chr}->{$count} = { start => $start, stop => $stop, id => $featid, fam => $line[4], sF => $line[5]};
    }
    print STDERR "Bed file loaded, $count lines\n";
    if ( defined( $bedFile ) ) {
        print STDERR "INFO $prog: read tab File $bedFile done !\n";
    }
}

sub tsdSize {
    my $self  = shift;
	$self->{tsd}->{"LTR/Copia"} = 5 ;
	$self->{tsd}->{"LTR/Gypsy"} = 5 ;
	$self->{tsd}->{"RC/Helitron"} = 0 ;
	$self->{tsd}->{"DNA/En-Spm"} = 3 ;
	$self->{tsd}->{"DNA/MuDR"} = 9 ;
	$self->{tsd}->{"Unassigned"} = 5 ;
	$self->{tsd}->{"DNA/HAT"} = 5 ;
	$self->{tsd}->{"LINE/L1"} = 9 ;

	$self->{tsd}->{"RLC"} = 5 ;
	$self->{tsd}->{"RLG"} = 5 ;
	$self->{tsd}->{"DHH"} = 0 ;
	$self->{tsd}->{"RST"} = 9 ;
	$self->{tsd}->{"DTC"} = 3 ;
	$self->{tsd}->{"DTM"} = 9 ;
	$self->{tsd}->{"DTA"} = 5 ;
	$self->{tsd}->{"DTH"} = 5 ;
}


#=============================================================================
# HELP
#=============================================================================
sub help {
my $prog = basename($0);
print STDERR <<EOF ;
### $prog ###
#
# CREATED:    2017-01-04
# LAST MODIF: $lastmodif
# AUTHOR:     Josquin DARON (Ohio State University, Slotkin Lab)
# VERSION:    $VERSION
# PURPOSE:    Insert sequences (TEs) and it corresponding TSD within single FASTA sequence.

USAGE: $prog -bed <bed> -db <FASTA> <FASTA>

       <FASTA> Single FASTA file such as one chr sequence.

       ### OPTIONS ###
       
       -bed <FILE> Bed file with customized format
                   Column 4: TE id of TE that will be inserted (and present in -db FASTA file) 
                   Column 6: superfamily id of TE that will be inserted (will decided TSD size).
                   Wicker classification accronymes is accepted (Wicker et al., 2007).

       -db <FILE>  Multi-FASTA file of TE sequences that will be inserted.

       -h|--help:   print this help

EOF

exit(1) ;
}

# END OF SCRIPT

