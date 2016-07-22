#!/usr/bin/perl
use strict;
use warnings;
use Bio::Community;
use Bio::Community::IO;
use Cwd;
use File::Path 'make_path';
#use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::Tools::Run::Alignment::MAFFT;
#use Bio::SimpleAlign;
use Bio::AlignIO;

#
#This script must be executed from within an active macqiime environment.
#will run core_diversity_analyses.py, convert biom table to qiime format for every biom table in the folder
#OTU names will be extracted from qiime formatted table and used to pull and align representative sequences from ref_seq.fna
#Reference sequences are aligned with mafft and saved in phylip format

#usage ./special_taxa.pl ref_seqs.fna TAXONOMY.txt

my $rep = shift;
my $tax = shift;
my $dir = getcwd();
my $fasdb = Bio::DB::Fasta->new($rep) or die "$rep: $!";
my @ids = $fasdb->get_all_primary_ids;


#get the family name from a *.biom file, assumed to be the last taxonomy before the file extension
#make new directory with family name
#save qiime formatted otu table in new dir
opendir (DIR, $dir) or die $!;
while (my $file = readdir(DIR)) {
    next unless (-f "$dir/$file");
    if ($file =~/^.+___(\w+)\.biom$/) {
        my @seq_array;
        my %all_names;
        my $outdir = "$dir/output/$1_0513_aln";
        my $outfile = "$1_otu_table.txt"; 
        make_path $outdir;
        open (my $log_out, ">>", "output/alingment_log.txt");
        #system ("core_diversity_analyses.py --recover_from_failure -i $file -o $outdir -m nema_map_0510.txt -e 100 -p cd_params.txt -a --nonphylogenetic_diversity") == 0 or die "system core_diversity failed: $?\n";
        system ("biom convert -i $file -o $outdir/$outfile --to-tsv") == 0 or warn "system biom convert failed: $?\n";

	    #get the otus and extract them from the ref_seq.fna file.  Currently must be converted to qiime format (see above system call).
		my $in = Bio::Community::IO->new( 
		    -file => "$outdir/$outfile",
		    -format => 'qiime'
		);
        my $meta = $in->next_metacommunity;
        $in->close;
		while (my $community = $meta->next_community) {
		    my $comm_otus = $community->get_all_members(); #'$comm_otus' is a reference to an array of hash references containing community info
		    foreach my $otu (@$comm_otus) {
                $all_names{$$otu{id}}++;
		    }
		}
#		my @bad = ('New.CleanUp.ReferenceOTU394269', 'New.CleanUp.ReferenceOTU182375', 'New.CleanUp.ReferenceOTU58854', 'New.CleanUp.ReferenceOTU643381', 'New.CleanUp.ReferenceOTU199913');
#		foreach my $bad_name (@bad) {
#			if (exists $all_names{$bad_name}) {
#				delete $all_names{$bad_name};
#				}
#			}
        my @names = sort keys %all_names;

        #Get the sequence for each OTU, get data for the otu from the tax file, align with mafft, save alignment
        my $num_otu = 0;
        foreach my $id (@names) {
            $num_otu++;
            my $seq = $fasdb->get_Seq_by_id($id);
            push (@seq_array, $seq);
            open (my $taxfh, "<", $tax) or die "can't open $tax: $!\n";
            open (my $outtax, ">>", "$outdir/$1_taxonomy.txt");
            while (my $line = <$taxfh>){
                if ($line =~ /^$id\t(.+)\s\d.+/) {
                    my @kpcofgs = split (/;\s__/, $1);
                    print $outtax $id,"\t", join("\t", @kpcofgs),"\n" ;
                }
            }
            close $taxfh;
            close $outtax;
        }
        print $log_out $1,"\t",$num_otu,"\n";
        #if two or more sequences do alignment
        if ($num_otu >1) {
	        my $factory = Bio::Tools::Run::Alignment::MAFFT->new;
	        my $aln = $factory->align(\@seq_array);
	        my $out = Bio::AlignIO->new(-file => ">$outdir/$1.phy", -format => 'phylip', -longid => 1, -interleaved => 0); {
	            $out->write_aln($aln);
	        }
	    } 
        else {
            my $out = Bio::SeqIO->new(-file => ">$outdir/$1.fas", -format => 'fasta');{
                $out->write_seq($seq_array[0]);
            }
        }
    }
}

__END__
