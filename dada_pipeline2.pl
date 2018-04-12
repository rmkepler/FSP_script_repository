use strict;
use warnings;
use File::Copy qw(copy);
use Parallel::ForkManager;
use Text::CSV;
use Cwd;
use File::Path qw(make_path remove_tree);
use Getopt::Long qw(GetOptions);


#must be executed in the directory containing Illumina output
#writes the first 12000 lines of raw miseq fastq.gz files to a single file to examine quality
my $usage = "perl $0 --plate SLURM_ARRAY_TASK_ID --outfolder OUTPUT-FOLDER --list NUMBERED-LIST-OF-DIRECTORIES --procs NUMBER-OF-PROCESSORS";

my $max_procs;
my $plate;
my $output;
my $list;
my $dir = getcwd();
my $workdir;
my $run;

GetOptions('plate=s' => \$plate, 'outfolder=s' => \$output, 'list=s' => \$list, 'procs=s' => \$max_procs) or die "Usage: $usage\n";

#open $list and scan for line starting with the number to find the directory
open (my $fh, '<', $list) or die "can't open $list: $!\n";
while (my $line = <$fh>) {
	chomp $line;
	my ($fnum, $fdir, $frun) = split(/\s+/, $line); #parses each line into three parts. The number of spaces separating each column doesn't matter
	if ($plate eq $fnum) {
		$workdir = $fdir;
		$run = $frun
		}	
}

close $fh;

my $outdir = "$dir/$output/$run";
my $cutfor = "$outdir/fdir";
my $cutrev = "$outdir/rdir";
my $cutout = "$outdir/cutout";
my $rcfor = "$outdir/rcfor";
my $rcrev = "$outdir/rcrev";
my $rcout = "$outdir/rcout";
make_path $outdir or die "couldn't create directory $outdir: $!\n";
make_path $cutfor or die "couldn't create directory $cutfor: $!\n";
make_path $cutrev or die "couldn't create directory $cutrev: $!\n";
make_path $cutout or die "couldn't create directory $cutout: $!\n";
make_path $rcfor or die "couldn't create directory $rcfor: $!\n";
make_path $rcrev or die "couldn't create directory $rcrev: $!\n";
make_path $rcout or die "couldn't create directory $rcout: $!\n";

my %well_hash;
opendir( my $work, "$workdir/Data/Intensities/BaseCalls/") or die $!;
while (my $fastq = readdir($work)) {
	next unless (-f "$workdir/Data/Intensities/BaseCalls/$fastq" && $fastq =~/^.+\.fastq\.gz$/);
	next if ($fastq =~ /^Undetermined.+\.fastq\.gz$/);
	if ($fastq =~ /^.+_(S{1}\d{1,2})_L001_R1_001\.fastq\.gz$/) {
    	$well_hash{$1}{for} = $fastq;
        }
    elsif ($fastq =~ /^.+_(S{1}\d{1,2})_L001_R2_001\.fastq\.gz$/) {
        $well_hash{$1}{rev} = $fastq;
        }
}
close $work;

my $pl = Parallel::ForkManager->new($max_procs);
my @well_keys = (sort keys %well_hash);
DATA_LOOP:
foreach my $well (@well_keys) {
	my $pid = $pl->start and next DATA_LOOP;
	system ("cutadapt -g ^CTTGGTCATTTAGAGGAAGTAA -G ^GCTGCGTTCTTCATCGATGC -q 22 --minimum-length 75 --discard-untrimmed -o $cutfor/$run.${well}_R1.fastq.gz -p $cutrev/$run.${well}_R2.fastq.gz $workdir/Data/Intensities/BaseCalls/$well_hash{$well}{for} $workdir/Data/Intensities/BaseCalls/$well_hash{$well}{rev} > $cutout/$run.${well}_out.txt"
			);
	system ("cutadapt -a GCATCGATGAAGAACGCAGC -A TTACTTCCTCTAAATGACCAAG --minimum-length 75 -o $rcfor/$run.${well}_R1.fastq.gz -p $rcrev/$run.${well}_R2.fastq.gz $cutfor/$run.${well}_R1.fastq.gz $cutrev/$run.${well}_R2.fastq.gz > $rcout/$run.${well}_out.txt"
    );
    $pl->finish;
}
$pl->wait_all_children;

remove_tree($cutfor, $cutrev);

system ("Rscript --vanilla dada_filter.R $outdir $rcfor $rcrev $run > $outdir/${run}_dada_filter.out 2>&1") == 0 or die "rscript system call failed: $?\n\n"

__END__
