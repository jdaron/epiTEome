INSTALL

    Dependencies
        All dependencies must be executable and findable in the user's PATH.

        perl (version 5.x): Generally installed in linux and mac OS by default. Expected to be installed at /usr/bin/perl

        perl lib (Bio::SeqIO, Tools::GFF, DB::Sam, Set::IntervalTree, Statistics::Descriptive)

        samtools (version 1.3.1 or higher)

        bedtools (version v2.26.0 or higher)

        ngsutils (version 0.5.7 or higher)

        segemehl (version 0.2.0 or higher)

        gzip/bzip2 Generally installed in linux and mac OS by default.

    Install
        EpiTEome is a perl program that does not need to be compiled.
        Make sure it is executable. For convenience, can be added to your PATH.
        It expects your perl installation to be at /usr/bin/perl.

    Tips
        All libraries can be installed using `perl -MCPAN -e 'install Bio::SeqIO'`
        ngsutils can be installed using `pip install ngsutils`

    Test environment
        This release of epiTEome (v1) was tested on Mac OSX (10.11.6), perl
        5.18.2, samtools 1.3.1, bedtools v2.26.0, ngsutils 0.5.7, segemehl 0.2.0

USAGE

    INDEX: Reference fasta file should be indexed in the segemehl index format.
           idxEpiTEome.pl will masked the 3’ edge of the LTR5 and 5’ edge of the LTR3.

        Usage: idxEpiTEome.pl —l [max read length] -gff <gff3> -t <target> —ref <fasta>

        <gff3>    TE annotation in gff3 format.
            
        <target>  list of TEid of interest

        <ref>     FASTA formated (.fa, .fna or fasta) genome file.

        -l        Maximum reads length present in FASTQ file.


    EPITEOME: Identify new insertion sites and their methylation level.

        Usage: epiTEome.pl [options] -gff <gff3> -t <target> —ref <fasta> -un <fastq>

        <gff3>    TE annotation in gff3 format.
            
        <target>  TEid list of interest

        <ref>     FASTA formated (.fa, .fna or fasta) genome file.

        <un>      FASTQ file of reads that failed to map the reference genome (unmapped reads)

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

GFF3 INPUT FILE FORMAT

    TE annotated features (mandatory): 
        - Column 3 (type) should be referred as 'te'
        - Column 9 (attributes) should have following list of attributes: ID=teid, sF=superfamily_name, fam=family_name.
    LTR annotated features (optional):
         - Column 3 (type) should be referred as LTR5 or LTR3
         - Column 9 (attributes) should have tags Parent=teid.

OUTPUT FILE FORMAT

    .newInsertionSite.tab
        1.  chrom - name of the chromosome or scaffold
        2.  chromStart - start position of feature containing new insertion site (0 base)
        3.  chromEnd - end of the feature containing new insertion site
        4.  name - feature name
        5.  mapping type [uniq|multi] - Feature has been identify using split-reads that uniquely map the reference sequence (uniq)
                                        or map the reference sequence multiple time (multi).  
        6.  strand
        7.  tsdStart - Start of the TSD (target-site duplication)
        8.  tsdEnd - End of the TSD
        9.  nubReads - Number of split-reads aligned to identify this feature
        10. family - TE family name
        11. teid - teid name

    .met.row.tab
        1. methylation context [CG|CHG|CHH]
        2. location [neo|te] - neo (at flanking DNA at newinsertion site) or te (at TE)
        3. edge [5|3|8] - 5prime (5), 3prime (3), both (8)
        4. nbCm - number of cytosine methylated
        5. nbR - number of reads mapped
        6. name - feature name
        7. teid - teid name

    .met.meta.tab
        1. methylation context [CG|CHG|CHH]
        2. location [neo|te] - neo (at flanking DNA at newinsertion site) or te (at TE)
        3. edge: [5|3|8] - 5prime (5), 3prime (3), both (8)
        4. window id
        5. methylation level (%)
        6. confidence interval (95%)

TEST
    
    Test / demonstration data for epiTEome.
    - Step 1: Indexing reference file
	   $idxEpiTEome.pl -l 85 -gff tair10TEs.gff3 -t subteid.lst -fasta Chr2.fasta 

    - Step 2: Run epiTEome analysis
	   $epiTEome.pl -gff tair10TEs.gff3 -ref Chr2.epiTEome.masked.fasta -un unmapped.fastq -t teid.lst 

    - Output: 4 different files (unmapped.newInsertionSite.tab, unmapped.newInsertionSite.sam, unmapped.met.meta.tab and unmapped.met.row.tab) will be automatically generated by epiTEome. To check if epiTEome worked successfully, those files could be compare to reference output files present in the test folder (refOutput_*).
