#!/usr/bin/perl

# Culls contigs shorter than 200 bp for NCBI submission


## Read fasta file

  $/ = undef;   			 # read whole file, not line-by-line

  $Fasta = $ARGV[0];

  ($Genome_ID = $Fasta) =~ s/_.+//;

  ($FastaOut = $Fasta) =~ s/_nh/_final/;

  open(FASTAOUT, '>', $FastaOut) || die "Can't create output file\n";

  open(FASTA, "$Fasta");

  $wholeFile = <FASTA>;

  @SeqRecords = split(/>/, $wholeFile); # split records based on '>' symbol

  shift @SeqRecords;			# remove empty record at beginning	

  foreach $Record (@SeqRecords) {

    $SeqNo ++;

    ($Current_header, $Seq) = split(/\n/, $Record, 2);

    $Seq =~ s/[^A-Za-z]//g;


    # print out modified sequence record if > 200 bp in length

    print FASTAOUT ">$Genome_ID"."_contig$SeqNo\n$Seq\n" if length($Seq) >= 200;

    $Warning = 'Sequence omitted from final assembly (length < 200)';

    if(length($Seq) < 200) {

      print "$Genome_ID"."_contig$SeqNo\t$Warning\n"

    }

  }

  close FASTA;

  close FASTAOUT;

  ## print individual genome completion statement

## print overall run completion statement

print "Contig culling finished.\n";

print "\n########################################################\n\n"
