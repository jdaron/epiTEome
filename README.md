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
        epiTEome assumes your perl installation to be at /usr/bin/perl.

    Tips
        All libraries can be installed using `perl -MCPAN -e 'install Bio::SeqIO'`
        ngsutils can be installed using `pip install ngsutils`

    Test environment
        This release of epiTEome (v1) was tested on Mac OSX (10.11.6), perl
        5.18.2, samtools 1.3.1, bedtools v2.26.0, ngsutils 0.5.7, segemehl 0.2.0

CITATION

    If you use epiTEome in your work, please cite one of the following:
    
USAGE

    INDEX: Reference fasta file should be indexed in the segemehl index format.
           Prior the indexing, idxEpiTEome.pl will mask the 3’ edge of the LTR5 
           and 5’ edge of the LTR3 to avoid multi-mapping read competition within 
           a single TE. Because LTRs of a single TE are identically duplicated at
           the TE insertion time, this light masking will prevent split-reads that
           could map at the junctions TE-flanking DNA to map inside the TE.

        Usage: idxEpiTEome.pl —l [max read length] -gff <gff3> -t <target> —ref <fasta>

        <gff3>    TE annotation in gff3 format.
            
        <target>  list of TEid of interest

        <ref>     FASTA formated (.fa, .fna or fasta) genome file.

        -l        Maximum read length present in FASTQ file.


    EPITEOME: Identify new TE insertion sites and quantify their methylation level from MethylC-seq datasets.

        Usage: epiTEome.pl [options] -gff <gff3> -t <target> —ref <fasta> -un <fastq>

        <gff3>    TE annotation in gff3 format.
            
        <target>  TEid list of interest

        <ref>     FASTA formated (.fa, .fna or fasta) genome file.

        <un>      FASTQ file of reads that failed to map to the reference genome (unmapped reads)

    OPTIONS
      EpiTEome Specific Options:
        -chop [integer] : read end length of chopped (defaut 25,30,40).
                          Usage of several lengths will improve epiTEome sensitivity. 
        -b    [integer] : number of TEs per batch (defaut 5000).
        -w    [integer] : window size for methylation metaplot analysis (defaut 10 bp)

      Alignment Options:
        -E    [integer] : segemehl max evalue (default:5)
        -p    [integer] : number of threads used in segemehl (defaut 1).
                          All other portions of epiTEome are single-threaded.

    OUTPUT
      epiTEome output 4 different files .newInsertionSite.tab, .newInsertionSite.sam, .met.meta.tab and .met.row.tab

GFF3 INPUT FILE FORMAT
     
    GFF3 input file follow the standard GFF3 format, except column 3 and 9 that have specific tags.
    TE annotated features (mandatory): 
        - Column 3 (type) should be referred to as 'te'
        - Column 9 (attributes) should have the following list of attributes: ID=teid, sF=superfamily_name, fam=family_name.
    LTR annotated features (optional):
         - Column 3 (type) should be referred to as LTR5 or LTR3
         - Column 9 (attributes) should have tag Parent=teid.

OUTPUT FILE FORMAT

    .newInsertionSite.tab: coordinate of non-reference TEs
        1.  chrom - name of the chromosome or scaffold
        2.  chromStart - start position of feature containing new insertion site (0 base)
        3.  chromEnd - end position of the feature containing new insertion site
        4.  name - feature name
        5.  mapping type [uniq|multi] - Feature has been identified using split-reads that uniquely map to the reference sequence (uniq)
        or map to the reference sequence multiple time (multi).  
        6.  strand
        7.  tsdStart - Start of the TSD (target-site duplication)
        8.  tsdEnd - End of the TSD
        9.  nubReads - Number of split-reads aligned to identify this feature
        10. family - TE family name
        11. teid - teid name

    .newInsertionSite.sam: standard sam aligment file diplaying split-reads aligment profile. Note that all 
                           split-reads susceptible to detect a non-reference TE will be store in this file,
                           allowing user to identify false negative predictions.

    .met.row.tab: methylation level at each cytosine position (used for barplot, Figure 4A)
        1. methylation context [CG|CHG|CHH]
        2. location [neo|te] - neo (at flanking DNA at newinsertion site) or te (at TE)
        3. edge [5|3|8] - 5prime (5), 3prime (3), both (8)
        4. nbCm - number of cytosines methylated
        5. nbR - number of reads mapped
        6. name - feature name
        7. teid - teid name

    .met.meta.tab: process methylation level for metaplot analysis (used for metaplot, Figure 4B)
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

    - Output: 4 different files (unmapped.newInsertionSite.tab, unmapped.newInsertionSite.sam, unmapped.met.meta.tab and unmapped.met.row.tab) will be automatically generated by epiTEome in the $CWD. To check whether epiTEome worked successfully, those files could be compared to reference output files present in the test folder (refOutput_*).
