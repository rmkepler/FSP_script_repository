#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Path 'make_path';
use Parallel::Loops;

#usage: ./script.pl MAPPING_FILE OUTPUT-FOLDER

#must be executed from an active macqiime environment in the directory containing Illumina output
#runs the join_paired_ends.py script for any *fastq.gz files matching criteria saving output to subdirectories matching the well number
#runs the convert_fastaqual_fastq.py script within each new directory
#parses master mapping file to individual mapping files in each new directory
#runs the split_libraries.py within each new directory

my $map = shift;
my $output = shift;
my $dir = getcwd();
my %well;
my $max_procs = 2;

#builds hash of well number keys and forward and reverse paired end file names
opendir (DIR, $dir) or die $!;
while (my $file = readdir(DIR)) {
    next unless (-f "$dir/$file");
    next if ($file =~ /^Undetermined.+\.fastq\.gz$/);
    #if ($file =~ /^.+_S(.+)_L001_R1_001\.fastq\.gz$/) {
    if ($file =~ /^.+-index(N(\d+))-(S(\d+))-.+-([A-H]{1}\d)_.+_L001_R1_001\.fastq\.gz$/) {
        $well{$5}{for} = $file;
        $well{$5}{index1} = $1;
        $well{$5}{index2} = $2;
        $well{$5}{index3} = $3;
        $well{$5}{index4} = $4;
    }
    #elsif ($file =~ /^.+_S(.+)_L001_R2_001\.fastq\.gz$/) {
    elsif ($file =~ /^.+-index(N(\d+))-(S(\d+))-.+-([A-H]{1}\d)_.+_L001_R2_001\.fastq\.gz$/) {
        $well{$5}{rev} = $file;
    }
}
closedir (DIR);

#This loop performs the work on *.fastq.gz files from Illumina
#The result is the creation of individual folders where qiime output is saved
#Performs tasks one at a time, with the final result contigs from paired end assembly

my $pl = Parallel::Loops->new($max_procs);
my @well_array = (sort keys %well);
$pl->foreach (\@well_array, sub{
    my $key = $_;
    if ($well{$key}{index2} <= 706) { #key less than or equal to 48, to filter half the wells.
    my ($fna, $qual, $map_head, $map_line);
    my $outdir = "$output/$key";
    my $outfile = "map_$key.txt";
    make_path $outdir  or die "couldn't create directory: $!\n";
    system ("join_paired_ends.py -f $well{$key}{for} -r $well{$key}{rev} -o $dir/$outdir/") == 0 or die "system failed: $?\n";
    #This secion converts from *.fastq to *.fna and *.qual files
    opendir (DIR, $outdir) or die "can't open directory: $!\n";
    while (my $jfile = readdir(DIR)) {
        next unless (-f "$outdir/$jfile");
        if ($jfile =~/fastqjoin.join.fastq/) {
            system("convert_fastaqual_fastq.py -c fastq_to_fastaqual -f $outdir/$jfile -o $outdir") == 0 or die "system failed: $?\n";
        }
    }
    closedir (DIR);

    #This section gets the file names for the *.fna and *.qual files     
    opendir (DIR, $outdir) or die "can't open directory: $!\n";
    while (my $sfile =readdir(DIR)) {
        next unless (-f "$outdir/$sfile");
        if ($sfile =~ /.+\.fna/) {
            $fna =$sfile;
        }
        elsif ($sfile =~ /.+\.qual/) {
            $qual =$sfile;
        }
    }
    closedir (DIR);

    #This section makes an individual mapping file specific to $key
    #runs the split_library script to perform demultiplexing
    open (IN, '<', $map) or die "can't open map file: $!\n";
    while (my $line = <IN>) {
        chomp $line;
        if ($line =~ /^#.+$/) {
            $map_head = $line;
        }
        elsif ($line =~ /^$key\s.+$/) {
            $map_line = $line;
        }
    }
    close IN;
    open (OUT, '>', "$outdir/$outfile");
    print OUT $map_head,"\n";
    print OUT $map_line,"\n"; 
    system ("split_libraries.py -m $outdir/$outfile -f $outdir/$fna -q $outdir/$qual -b 0 -z truncate_only -p -o $outdir") == 0 or die "system failed: $?\n";
}
});


__END__
