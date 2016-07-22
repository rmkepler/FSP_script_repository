#!/usr/bin/perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;

#usage genbank_taxonomy2csv.pl infile outfile

my $gbfile = shift;
my $outfile = shift;

my $in = Bio::SeqIO->new(-file => $gbfile);

open (my $out, ">", $outfile);

while (my $seq = $in->next_seq ) {
    my $accession = $seq->accession_number;
    my $dis_id = $seq->display_id;
    my @classification = $seq->species->classification;
    print $out $accession,",",$dis_id,",",$dis_id,"_bb", ",",join(",", @classification),"\n";
}

close $out;
__END__