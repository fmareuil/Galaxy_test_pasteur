#! /usr/bin/perl
# Version : 1.0
# setuppromoter.pl
# builds a fasta file of sequences from a reference genome using promoter region coordinates.
#
# The coordinate file can be in two formats:
#	chr:start-end
#	     or
#	rank	chr	start	end	strand	rest	(tab separated columns)
#

use Getopt::Std;
#use Bio::SeqIO;	# read fasta file directly is faster that Bio::SeqIO....

%seq = ();	# sequence cache

# this should probably be in a genomes_conf.pm file as it is also used by findpromoter.pl....
@genomes = ( {species => 'Hs', 
		Species => "Homo sapiens", 
		suffix => 'hs',
		annfile => "$anndir/human.gff3", 
		ann => {},
		genomefile=>"/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/hs_ncbi37.fa",
		blastdb => "/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/hs_ncbi37",
		genome => ()
		},
		{species => 'Mm', 
		Species => "Mus musculus", 
		suffix => 'mm',
		annfile => "$anndir/mouse.gff3", 
		ann => {},
		genomefile=>"/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/mm_ref_ncbi37.fa",
		blastdb=>"/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/mm_ref_ncbi37",
		genome => ()
		},
		{species => 'Md', 
		Species => "Monodelphis domestica", 
		suffix => 'md',
		annfile => "$anndir/opossum.gff3", 
		ann => {},
		genomefile=>"/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/mdm_ref.fa",
		blastdb=>"/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/mdm_ref",
		genome => ()
		},
);
$sidx = 1;	# default source genome is the mouse
$tidx = 0;	# default target genome is human...

@cols = (1,2,3);	# default columns for chromosome, start, end (in that order) (can be overiden by -C option)
@Species = ();

# make list of supported source/target species
for $elem (@genomes) {
	push @Species,$elem->{'species'}." (".$elem->{'Species'}.")"; # for example: "Hs (Homo sapiens)"
}

$e_value = "1e-10";
$blastcmd = 'blastall';
$pbblastcmd = 'pb blastall';

@options = (
["s","Source species", "species from which the promoters/adapter seq originate. Can be (".join(" | ",@Species).") [".$genomes[$sidx]->{'species'}."]"],
#["t","Target species", "species in which to search for promoters/adapter seq. Can be (".join(" | ",@Species).") [".$genomes[$tidx]->{'species'}."]"],
#["C","csv columns","for a CSV file, the columns where is the chromosome, start coordinate, and end coordinate. [". join(",",@cols)."]"],
["P","project name","Project name. The fasta file will be the project name with '.fa' appended. [$ARGV[0]]"],
#["b", undef,"Run the blast command."],
#["p", undef,"Run the blast command with Paracel Blast."],
#["e","e-value","Maximum BLAST e-value. [$e_value]"],
);
@otherOptions = (
	["coordinate.txt","Promoter region coordinates for the source species.\n\t2 formats possible:\n\t\tchr:start-end\n\t\tid chr start end strand rest"]
);


sub main::HELP_MESSAGE()
{
	($fh,@rest) = (@_);
	usage(undef);
}

$Getopt::Std::STANDARD_HELP_VERSION = 1;	# tell getopts to die if error

$opt = "h";

foreach $o (@options) {
	$opt .= $o->[0];
	$opt .= ':' if defined($o->[1]); # parameter takes an argument
};

@args = ();
foreach $arg (@ARGV) {
	$arg = "'$arg'" if $arg =~ m/\s/;
	push @args,$arg;
}
$cmdline = "$0 ".join(" ",@args);

getopts($opt);
usage(undef) if ($opt_h);

$sidx = findidx($opt_s) if $opt_s;
usage("Invalide source genome indicator ('$opt_s')") if $sidx == -1;		
print "Source genome: ".$genomes[$sidx]->{'Species'}."\n";
if ($opt_C) {
	@tcols = split/,\s*/,$opt_C;
	usage("Invalid number of columns: need 3 column numbers for chromosome, start, and end locations.") if @tcols ne 3;
	$msg = ""; $haserr = 0;
	for($i=0; $i < 3; ++$i) {
		unless ($tcols =~ m/^\d+$/) {
			$msg .= "Column ".($i+1)."is not valide (must be numeric)\n";
			$haserr = 1;
		}
		usage($msg) if $haserr;
	}
	@cols = @tcols;
}
usage("No coordinate file specified.") unless @ARGV == 1;

#@tcol = sort @cols;
#$maxcol = $tcol[2];
$maxcol = 5;
$| = 1;

#$in = Bio::SeqIO->new(-file=>$genomes[$sidx]->{'genomefile'}, -format => 'fasta');
#print "Read genome ".$genomes[$sdx]->{'genomefile'}.".. ";
#%seq= ();
#$chr = undef;
#while($seq = $in->next_seq()) {
#	($chr = $seq->id) =~ s/^chr//i;
#	print "$chr ";
#	$seq{$chr} = $seq->seq();
#}
#print "\n";
open SEQ,$genomes[$sidx]->{'genomefile'} or die "Cant open genome sequence file '".$genomes[$sidx]->{'genomefile'}."' : $!\n";
print "Read genome ".$genomes[$sidx]->{'genomefile'}.".. ";
while(<SEQ>) {
	chomp;
	if (m/^>(\w+)/) {
		$id = $1;
		print "$id ";
		$seq{$id} = "";
		next;
	}
	tr/AGTCXNagtcxn//cd;
	$seq{$id} .= uc $_;
}
close SEQ;
print "\n";

$coordFile = $ARGV[0];
($project = $coordFile) =~ s/\.(\w*)$//;
$project = $opt_P if $opt_P;

open IN,$coordFile or die "Cant open coordinate file '$coordFile' : $!\n";
$format = 1;		# column format by default
while(<IN>) {
	chomp;
	next if m/^#/;
	if (m/^\s*(\w+):(\d+)-(\d+)/) {
		$format = 0;
		print "Cordinate format is 'chr:start-end'\n";
		last;
	}
	(@flds) = split /\s/;
	usage("Not enough columns in data for maximum number of columns ($maxcol).") if @flds < $maxcol;
	print "Coordinate format is 'id chr start end strand'\n";
	last;
}
seek IN,0,0;

$nSeq = 0;
open OUT,">$project" or die "Cant create fasta file '$project' : $!\n";
while(<IN>) {
	chomp;
	next if m/^#/;
	++$nSeq;
	if ($format) {
		($id,$nchr,$start,$end,$strand,$rest) = split /\s/,$_,6;
		$id = $nSeq if ($id+0 == 0);
	} else {
		if (m/^\s*(\w+):(\d+)-(\d+)/) {
			($nchr,$start,$end) = ($1,$2,$3);
			$id = $nSeq;
			$strand = '+';
		} else {
			die ("Unable to parse coordinatesi ($_) at line $. of coordinate file '$coordFile'.\n");
		}
	}
	($chr = $nchr) =~ s/^chr//i;
	die ("Unknown chromosome '$nchr' at line $. in coordinate file '$coordFile'.\n") unless (exists($seq{$chr}));
	$id = sprintf("rank%5.5d",$id);
	print OUT ">$id $nchr|$start|$end|$strand\n";
	printSeq(*OUT{IO}, $seq{$chr},$start-1,$end-1,$strand);
}
close OUT;

sub printSeq {
	my ($io, $seq, $start, $end, $strand) = @_;

	# $start,$end are 0 based....
	($start,$end) = ($end, $start) if $start > $end;
	my $subseq = substr($seq, $start, $end-$start+1); 
	$subseq =~ tr/ATGCXNagtcxn//cd;
	if ($strand eq '-' || $strand == -1) {
		$subseq =~ tr/ATGCatgc/TACGtacg/;
		$subseq = reverse $subseq;
	}
	print $io formatSeq($subseq);
}
sub formatSeq {
	my ($subseq) = @_;
	my $len = length($subseq);
	my $idx = 0;
	my $str = "";
	while($idx < $len) {
		$need = $len-$idx;
		$need = 60 if $need > 60;
		$str .= substr($subseq,$idx,$need)."\n";
		$idx += $need;
	}
	$str
}
sub usage {
	$msg = shift;
	print "$msg\n" if defined $msg;
	print "Usage: $0";
	foreach $o (@options) {
		print " [-$o->[0]";
		print " $o->[1]" if defined $o->[1];
		print "]";
	}
	print " [-h]";
	foreach $o (@otherOptions) {
		print " $o->[0]";
	}
	print "\n"; # end of command line
	foreach $o (@options) {
		print "\t-$o->[0]\t$o->[2]\n";
	}
	print"\t-h\tThis useful help message\n";
	foreach $o (@otherOptions) {
		print "\t$o->[0]\t$o->[1]\n";
	}
	exit(1);
}
	
sub findidx {
	my $id = lc shift;
	for (my $idx = 0; $idx < @genomes; ++$idx) {
		return $idx if $genomes[$idx]->{'suffix'} eq $id;
	}
	return -1;
}



