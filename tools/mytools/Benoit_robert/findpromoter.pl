#! /usr/bin/perl 

# Version : 1.0
# program to search for a promoter motif and optionally, an adapter sequence from a source species on a target species
#	limits of the motif and adapter to it's nearest gene can be changed
#
# input files are blast output report of the sources species sequences blasted against the target genome
# the source species sequences must be fasta files, the fasta identifier must be in the form:
#		>id chr|start|last|strand
#			,where id is an identifier which appears on the report
#			chr|start|last|strand is the location/strand of the promoter region on the source species genome
# the 'chr' (chromosome) should be a short name like '1', '12', 'X', or 'Y' and must match the chromosome in the 
#		annotation GFF3 files
# the fasta id in the genome files should also be a short name for the chromosome and should also match the chromosome
#		name in the GFF3 annotation files
#
# the program produces a tabulation separated results file that contains:
#  the source species sequence identificator (from the fasta id from the source species sequences, the localisation of the region,
#  the location and sequence of the promoter, the location of the adapter, and the nearest gene. This information is optional for
#  the source species but included for the target species.
#
use Bio::SearchIO;
use Bio::SeqIO;
use Getopt::Std;
use Date::Format;

$minident = .90;
$minlen = 100;
$gaplimit = .20;		# max gap difference between two HSPs (in % > query length)

# defaults for mismatch from a consensus motif
$Promoter="TGGGTGGTC";					# Dai consensus motif
$Promotermismatch = 1;
$Promotermaxmm = 3;

# default for an array of Perl RE
#@Promoter = ("[TC]G[GC]GT[GA][GT].[ATC]");	# Dai motifs (Perl RE)
#$Promoter = "[CG][TC]G[GC]GT[GA][GT].[ATC][CT]";	# Dai + Vokes motifs (Perl RE)

# Dai + Vokes motifs (No Perl RE)
@Promoter = (qw(TGGGTGGAT CGGGTGGCA CGGGTGGTC GGGGTGGCA TGGGTGTTT TGCGTGGTA GGGGTGGGA TGGGTGGGT));

# by default, we use Perl RE
$PromoterUseMM = 0;

@Adpt = ("AATT");
$Adptmismatch = 0;
$Adptmaxmm = 1;

$PromoterAdptDist = 100;
$PromoterAdptMaxDist = 1000;

$annLimit = 10000;
$doAnn = 0;

$extraSeq = 0;
$maxExtraSeq = 10000;

%Genome = ();

$anndir = '/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/gtf';	# location of the GFF3 files
@genomes = ( {species => 'Hs', 
		Species => "Homo sapiens", 
		suffix => 'hs',
		annfile => "$anndir/human.gff3", 
		ann => {},
		genomefile=>"/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/hs_ncbi37.fa",
		genome => ()
		},
		{species => 'Mm', 
		Species => "Mus musculus", 
		suffix => 'mm',
		annfile => "$anndir/mouse.gff3", 
		ann => {},
		genomefile=>"/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/mm_ref_ncbi37.fa",
		genome => ()
		},
		{species => 'Md', 
		Species => "Monodelphis domestica", 
		suffix => 'md',
		annfile => "$anndir/opossum.gff3", 
		ann => {},
		genomefile=>"/pasteur/galaxy/galaxy-pasteur/projets/genomes/promoteur_ljm/mdm_ref.fa",
		genome => ()
		},
);
sub findidx {
	my $id = lc shift;
	for (my $idx = 0; $idx < @genomes; ++$idx) {
		return $idx if $genomes[$idx]->{'suffix'} eq $id;
	}
	return -1;
}
$suffix = "fa";	# sequence suffix

@Species = ();
# make list of supported source/target species
for $elem (@genomes) {
	push @Species,$elem->{'species'}." (".$elem->{'Species'}.")"; # for example: "Hs (Homo sapiens)"
}

use constant {
	DEB_HSP => 1,
	DEB_OTH => 2,
	DEB_INFO => 4
};
$debug = 0; #DEB_INFO;
$tidx = 0;
$sidx = 1;

@options = (
["s","Source species", "species from which the promoters/adapter seq originate. Can be (".join(" | ",@Species).") [".$genomes[$sidx]->{'species'}."]"],
["t","Target species", "species in which to search for promoters/adapter seq. Can be (".join(" | ",@Species).") [".$genomes[$tidx]->{'species'}."]"],
["G","Promoter sequence","are the Promoter sequences (separate by '|') or unique sequence with mismatches [".join("|",@Promoter)."]"],
["g","# Promoter mismatches" ,"number of Promoter mismatches [$Promotermismatch]"],
["M","Adpt sequence","is the Adpt sequence [".join("|",@Adpt)."]"],
["m","# Adpt mismatches","number of Adpt mismatches [$Adptmismatch]"],
["d","max distance Adpt/Promoter","maximum distance between Promoter motif and Adpt motif [$PromoterAdptDist]"],
["x",undef,"add Promoter/Adpt location from the mouse"],
["u",undef,"use ".($PromoterUseMM ? "Perl RE" : "mismatches"). " [".($PromoterUseMM ? "mismatches" : "Perl RE")."]"],
["a",undef,"list nearest mouse/human/opossum gene (within '-p' threshold) to Promoter sequence"],
["e","# extra sequence","add [$extraSeq] bp (maximum $maxExtraSeq) upstream/downstream of mapped hit"],
#["p","threshold","threshold in bp [$annLimit]"],
#["A","mouse GFF annotations","mouse annotations to use (in GFF format) [$ann]"],
#["H","human GFF annotations","human annotations to use (in GFF format) [$annhs]"],
#["P","opossum GFF annotations","opossum annotations to use (in GFF format) [$annmd]"],
["w","threshold","threshold for nearest gene (in bp) [$annLimit]"],
["i","blast identity","Blast hit minimum identity (%)[$minident]"],
["l","blast length", "Blast hit minimum length (bp) [$minlen]"],
["o","outfile","report output file"],
["S","source_seq.fa", "Source sequence file in FASTA format (FASTA file created by 'setuppromoter')"],
["D","debug level","Debug level 1: hsp parse 2:other info 4: chit/chat[$debug]"]
);
@otherOptions = (
	["blast.out ...", "output file(s) from Blast of source sequences to target genome"]
);

$opt = "h";

foreach $o (@options) {
	$opt .= $o->[0];
	$opt .= ':' if defined($o->[1]); # parameter takes an argument
};

#getopts('chG:g:M:m:aA:b:o:s:i:l:');
#die $usage if $opt_h;
@args = ();
foreach $arg (@ARGV) {
	$arg = "'$arg'" if $arg =~ m/\s/;
	push @args,$arg;
}
$cmdline = "$0 ".join(" ",@args);

getopts($opt);
usage(undef) if ($opt_h);

$debug = $opt_d if $opt_d;
$PromoterUseMM ^= 1 if $opt_u;
if ($PromoterUseMM) {
	$Promotermismatch = $opt_g if $opt_g;
	$Promoterlen = length($Promoter);
	$Promoterrev = reverse $Promoter;
	$Promoterrev =~ tr/[]ATGCatgc/][TACGtacg/;
} else {
	@Promoter = split /:/,$opt_G if $opt_G;
	@Promoterrev = ();
	foreach $gli3 (@Promoter) {
		$gli3rev = reverse $gli3;
		$gli3rev =~ tr/[]ATGCatgc/][TACGtacg/;
		push @Promoterrev,$gli3rev if $gli3rev ne $gli3;
	}
}

$Adpt = $opt_m if $opt_m;
@Adpt = split/:/,$opt_m if $opt_m;
@Adptrev = ();
foreach $msx (@Adpt) {
	$msxrev = reverse $msx;
	$msxrev =~ tr/[]ATGCatgc/][TACGtacg/;
	push @Adptrev,$msxrev if $msxrev ne $msx;
}
	
$Adptmismatch = $opt_M if $opt_M;
$PromoterAdptDist = $opt_d if $opt_d;

$sidx = findidx($opt_s) if $opt_s;
$tidx = findidx($opt_t) if $opt_t;
if ($sidx == -1 || $tidx == -1) {
	$msg = "";
	$msg .= "Invalide source genome indicator ('$opt_s')" if $sidx == -1;
	$msg .= ($tidx == -1 ? "\n" : "")."Invalide target genome indicator ('$opt_t')" if $tidx == -1;
	usage ($msg);
}
$doAnn = $opt_a if $opt_a;
#$ann = $opt_A if $opt_A;
#$annhs = $opt_H if $opt_H;
#$annmd = $opt_P if $opt_P;
$annLimit = $opt_w if $opt_w;
$minident = $opt_i if defined $opt_i;
$minident = $minident / 100. if $minident > 1;
$minlen = $opt_l if defined $opt_l;
$extraSeq = $opt_e if defined $opt_e && $opt_e > 0 && $opt_e < $maxExtraSeq;
$out = \*STDOUT;

$| = 1;

#if ($opt_S) {
#	$suffix = $opt_S;
#	$suffix =~ s/^\.//;
#}

$out->autoflush(1);	# flush the output file....

usage("Must use different source and target genomes (source and target species are '".$genomes[$sidx]->{'Species'}."')") if $sidx == $tidx;
usage("Promoter / Adpt distance must be between 1 and $PromoterAdptMaxDist") if $PromoterAdptDist < 1 || $PromoterAdptDist > $PromoterAdptMaxDist;
usage("Promoter mismatch must be between 0 and $Promotermaxmm") if $opt_g < 0 || $opt_g > $Promotermaxmm;
usage("Adpt mismatch must be between 0 and $Adptmaxmm") if $opt_g < 0 || $opt_g > $Adptmaxmm;
usage("blast process must be 1, 2, or 3") if defined($opt_b) && ($opt_b < 1 || $opt_b > 3);
usage("Missing blast results directory") unless @ARGV;
usage("percent identique option must be a percent between 0 and 100% (got $opt_i)") if $minident < 0 || $minident > 1;
usage("minimum length must be positive (got $opt_l)") if $minlen < 0;
usage("Missing source genome sequence file") unless defined($opt_S) && -r $opt_S;

# use if want to pass a directory instead of 1 or more blast results file
#@files = glob "$ARGV[0]/*$suffix";
#usage("No such file $ARGV[0]/*.$suffix") unless @files;

#@files = @ARGV;
$file = $ARGV[0];	# only want one file....
$err = 0;
$msg = "";
#foreach $file (@files) {
	unless (-r $file) {
		$msg .= ($err ? "\n" : "")."No such file '$file'";
		$err = 1;
	}
#}
usage($msg) if $err;

if ($opt_o) {
	open OUT, ">".$opt_o or die "Cant create '$opt_o' : $!\n";
	$out = \*OUT;

	@lt = localtime(time);
	$sspecies = $genomes[$sidx]->{'species'};
	$tspecies = $genomes[$tidx]->{'species'};
	print $out "# Run ".strftime("%d/%b/%Y %H:%M:%S",@lt)."\n";
	print $out "# Command Line : $cmdline\n";
	@hdr = ("#id","Src sp.","$sspecies Chr","$sspecies Start","$sspecies End","$sspecies Strand");
	push @hdr,"Promoter $sspecies loc", "Promoter $sspecies strand","Promoter $sspecies motif","Adpt $sspecies loc","Adpt $sspecies motif" if $opt_x;
	push @hdr,"Nearest $sspecies gene","Ensembl id","Gene start","Gene end","Gene strand","distance to gene" if $opt_a && $opt_x;
	push @hdr,"Promoter species","$tspecies chr","$tspecies loc","$tspecies strand", "$tspecies motif", "$tspecies Adpt loc","Adpt motif";
	push @hdr,"$tspecies narest gene","Ensembl id","Gene start","Gene end","Gene strand","distance to gene" if $opt_a;
	print $out join("\t",@hdr)."\n";
}

$anncnt = 0;
if ($opt_a) {
	readgff($genomes[$sidx]->{'annfile'},$genomes[$sidx]->{'ann'});
	readgff($genomes[$tidx]->{'annfile'},$genomes[$tidx]->{'ann'});
	$ann = $genomes[$tidx]->{'ann'};
	$annsrc = $genomes[$sidx]->{'ann'};
}
%coord = ();

# fetch the genome sequences

#$in = Bio::SeqIO->new(-file=>$genomes[$tidx]->{'genomefile'}, -format => 'fasta');
#print "Read genome ".$genomes[$tidx]->{'genomefile'}.".. ";
#%Genome= ();
#$chr = undef;
#while($seq = $in->next_seq()) {
#	($chr = $seq->id) =~ s/^chr//i;
#	print "$chr ";
#	$Genome{$chr} = $seq;
#}
open SEQ,$genomes[$tidx]->{'genomefile'} or die "Cant open genome sequence file '".$genomes[$tidx]->{'genomefile'}."' : $!\n";
print "Read genome ".$genomes[$tidx]->{'genomefile'}.".. ";
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


$, = "\t";
$\ = "\n";
foreach $blastfile (@files) {
	%coord = ();
	#($prefix = $blastfile) =~ s/\.\w*//;
	#print "File '$blastfile' prefix '$prefix'\n";
	#open IN,$file or die "Cant open '$file' : $!";
	#$fafile = undef;
	#foreach $ext (qw(fa fasta fast seq fsa nt aa fna faa)) {
	#	if (-r "$prefix.$ext") {
	#		$fafile = "$prefix.$ext";
	#		last;
	#	}
	#}
	$fafile = $opt_S;
	unless (defined $fafile) {
		print STDERR "$file: Cant determine name of fasta sequence file.";
		next;
	}
	$in = Bio::SeqIO->new(-file=>$fafile, -format=> 'fasta');
	print STDERR "seqfile $fafile";
	%qseq = ();
	# read query sequences
	while($seq = $in->next_seq()) {
		$id = $seq->id;
		warn "Already seen sequence '$id'" if (exists $qseq{$id});
		$qseq{$id} = $seq;
		if (!$c_opt && $seq->desc =~ m/(\w+)\|(\d+)\|(\d+)\|(.)/) {
			$coord{$id} = "$1|$2|$3|$4";	# chr | start | end | strand
		} else {
			print STDERR "Cant parse qseq (".$seq->id.") desc (".$seq->desc.")" unless $c_opt;
		}
	}

	print STDERR "blastfile $blastfile";
	$blast = Bio::SearchIO->new(-file=>$blastfile,-format=>'blast');
	while($result=$blast->next_result()) {
		$query_name = $result->query_name();
		$query_length = $result->query_length();
		print "query name '$query_name'" if $debug & DEB_INFO;
		print STDERR "query name '$query_name' sequence not found" unless exists $qseq{$query_name};
		if ($query_name eq 'rank02098') {
			print "Got it";
		}
		$found = 0;
		$hitcnt = 0;
		if($hit = $result->next_hit()) {
			$hit_name = $hit->name();
			#print STDERR "hit_name '$hit_name'";
			($tchr = $hit_name) =~ s/\s//g;
			unless (defined $tchr && exists $Genome{$tchr}) {
				print STDERR "Query name '$query_name': Hit '$hit_name' has invalid chromosome (".(defined($tchr) ? $tchr : "UNDEF").")";
				next;
			}

			if ($hit->hsps() > 1) {
				($qstart,$qend,$sstart,$send) = (undef,undef,undef,undef);

				@strand = ();
				%strand = ();
				foreach $hsp ($hit->hsps()) {
					$strand = $hsp->strand('hit');
					for ($idx = 0; $idx < @strand; ++$idx) {
						last if ($strand[$idx] == $strand);
					}
					push @strand,$strand if $idx >= @strand; 
					$strand{$strand} = [] unless defined $strand{$strand};
					push @{$strand{$strand}}, $hsp;
				}
				($rqstart,$rqend,$rsstart,$rsend,$rstrand) = (0,0,0,0,0);
				foreach $strand (@strand) {
					print "\tstrand $strand has ".scalar(@{$strand{$strand}})." HSPS" if $debug & OPT_HSP;
					if (@{$strand{$strand}} > 1) {
						$besthsp = $strand{$strand}->[0];
						($bqstart,$bqend,$bsstart,$bsend) = ($besthsp->start('query'),$besthsp->end('query'),
													$besthsp->start('hit'),$besthsp->end('hit'));
						@hsps = sort {$a->start('query') <=> $b->start('query')} @{$strand{$strand}};

						($qstart,$qend,$sstart,$send) = ($bqstart,$bqend,$bsstart,$bsend);
						foreach $hsp(@hsps) {
							($tqstart,$tqend,$tsstart,$tsend) = ($hsp->start('query'),$hsp->end('query'),
																$hsp->start('hit'),$hsp->end('hit'));
							next if $bqstart == $tqstart && $bqend == $tqend && $bsstart == $tsstart && $bsend == $tsend;
							if ($debug & OPT_HSP) {
								print "\t   HSP $strand q ($tqstart,$tqend) h ($tsstart,$tsend)";
								print "\t\tcur q ($qstart,$qend) h ($sstart,$send)";
							}
							if ($tqend < $qstart || $tqstart > $qend) {
								if ($tqend < $qstart) {
									$qdiff = $qstart - $tqend;
										$sdiff = $strand == 1 ? $sstart - $tsend : $tsstart - $send;
								} else {
									$qdiff = $tqstart - $qend;
									$sdiff = $strand == 1 ? $tsstart - $send : $sstart - $tsend;
								}
								if ($qdiff > 0 && $sdiff > 0 && $sdiff < ($qdiff * (1.0 + $gaplimit))) {
									$qstart = $tqstart if $tqend < $qstart;
									$qend = $tqend if $tqend > $qend;
									if ($strand == 1) {
										$sstart = $tsstart if $tsstart < $sstart;
										$send = $tsend if $tsend > $send;
									} else {
										$sstart = $tsstart if $tsstart < $sstart;
										$send = $tsend if $tsend > $send;
									}
									print "\tadd HSP $strand q ($qstart,$qend) h ($sstart,$send)" if $debug & OPT_HSP;
								} else {
									print "qdiff $qdiff sdiff $sdiff qdiff+gaplimit ".($qdiff* (1.0+$gaplimit)) if $debug & OPT_HSP;
								}
							}
							# deal with overlapping contig (but only if target coordinates are close...)
						}
					} else {
						$hsp = $strand{$strand}->[0];
						($qstart,$qend,$sstart,$send) = ($hsp->start('query'),$hsp->end('query'),
														$hsp->start('hit'),$hsp->end('hit'));
						print "\tstrand $strand 1HSP Qlen $query_length qstart,qend ($qstart..$qend) sstart,send ($sstart,$send)" if $debug & OPT_HSP;
					}
					($rqstart,$rqend,$rsstart,$rsend,$rstrand) = ($qstart,$qend,$sstart,$send,$strand) if $qend-$qstart > $rqend - $rqstart;
				}
			} else {
				$hsp =$hit->next_hsp();
				($rqstart,$rqend,$rsstart,$rsend,$rstrand) = ($hsp->start('query'),$hsp->end('query'),
															$hsp->start('hit'),$hsp->end('hit'),$hsp->strand('hit'));
				print "\t1HSP Qlen $query_length qstart,qend ($rqstart..$rqend) sstart,send ($rsstart,$rsend) strand $rstrand " if $debug & OPT_HSP;
					
			}
			$rqlen = $rqend-$rqstart+1;
			$rslen = $rsend-$rsstart+1;
			$pcdiff = $rqlen/$rslen;
			if ($rqlen/$query_length < .50) {
				print "\t\tReject: Query not long enough" if $debug & DEB_INFO;
			} else {
				print "\t\tAccept: qlen $query_length qstart,qend ($rqstart..$rqend) qlen $rqlen sstart,send ($rsstart,$rsend) slen $rslen strand $rstrand ".
					(sprintf("%3.2f %%",100.*$pcdiff))if $debug & DEB_INFO;
				print "\t\tLength < .80" if $pcdiff <= .80 && $debug & DEB_INFO;
				print "\t\tLength > 1.20" if $pcdiff >= 1.20 && $debug & DEB_INFO;
				print "\t\tBest hit Strand change" if $rstrand != $strand[0] && $debug & DEB_INFO;

					
				$qseq = uc $qseq{$query_name}->subseq($rqstart,$rqend);
				if ($extraSeq) {
					$rsstart -= $extraSeq;
					$rsstart = 1 if $rsstart < 1;
					$rsend += $extraSeq;
					$rsend = $Genome{$tchr}->length if $rsend > $Genome{$tchr}->length;
				}
				$tseq = uc $Genome{$tchr}->subseq($rsstart, $rsend);

				if (exists $coord{$query_name}) {
					$hasCoord = 1;
					($mmChr,$mmStart,$mmEnd,$mmStrand) = split (/\|/,$coord{$query_name},4);
				} else {
					$hasCoord = 0;
					($mmChr,$mmStart,$mmEnd,$mmStrand) = (undef,undef,undef,undef);
				}
				($mPromLoc,$mPromStrand,$mPromMatch,$mAdptLoc,$mAdptSeq) = (undef,undef,undef,undef,undef);
				($mEnsGene, $mEns, $mEnsStart, $mEnsEnd, $mEnsStrand, $mDiff) = (undef,undef,undef,undef,undef,undef);
				($EnsGene, $Ens, $EnsStart, $EnsEnd, $EnsStrand, $Diff) = (undef,undef,undef,undef,undef,undef);

				if ($opt_x) {	# want (mouse) query seq Promoter
					@qpos = findProm($qseq);
					print "findProm Mm ".scalar(@qpos) if $debug & DEB_INFO;
					($mLoc,$mStrand,$mMatch,$mAdptLoc,$mAdptSeq) = (undef,undef,undef,undef,undef);
					if (@qpos == 0) {
						++$stat{'MM_NOT_FOUND'};
					} else {
						if (@qpos == 1) {
							++$stat{'MM_FOUND_Promoter'};
							($mPromLoc,$mPromStrand,$PrommMatch) = (@{$qpos[0]});
							($mAdptLoc,$mAdptSeq) = findAdpt($qseq,$mPromLoc);
							$mAdptLoc += ($hasCoord ? ($mmStart-1) : 0) if defined $mAdptLoc;
							print "findAdpt Mm $mAdptLoc PromLoc $mPromLoc mStart ".(defined($mmStart) ? $mmStart : "UNDEF") if $debug & DEB_INFO;
							$mPromLoc += ($mmStart-1) if $hasCoord;
							($mEnsGene, $mEns, $mEnsStart, $mEnsEnd, $mEnsStrand, $mDiff) = nearestgene($mmChr,$mPromLoc,\%annsrc) if $opt_a && $hasCoord;
							($mLoc,$mStrand,$mMatch) = ($mPromLoc,strandSymbol($mPromStrand),$PrommMatch);

						} else {
							++$stat{'MM_FOUND_Multi_Promoter'};
							@mPromLoc = (); @mPromStrand=(); @mPromMatch=(); @mAdptLoc=(); @mAdptSeq=();
							foreach my $loc(@qpos) {
								my $mLoc = $loc->[0];
								($mAdptLoc,$mAdptSeq) = (undef,undef);
								($mAdptLoc,$mAdptSeq) = findAdpt($qseq,$mLoc);
								print "findAdpt(multi) Mm $mAdptLoc PromLoc $mLoc mStart ".(defined($mmStart) ? $mmStart : "UNDEF") if $debug & DEB_INFO;
								if (defined $mAdptLoc) {
									$mAdptLoc += ($hasCoord ? ($mmStart-1) : 0);
									push @mAdptLoc,$mAdptLoc;
									push @mAdptSeq,$mAdptSeq;
								}
								$mLoc += ($mmStart-1) if $hasCoord;
								push @mPromLoc,$mLoc;
								push @mPromStrand,strandSymbol($loc->[1]);
								push @mPromMatch,$loc->[2];
								$mDiff = $annLimit;
								($mEnsGene, $mEns, $mEnsStart, $mEnsEnd, $mEnsStrand, $mDiff) = (undef,undef,undef,undef,undef,undef);
								if ($opt_a && $hasCoord) {
									($tmEnsGene, $tmEns, $tmEnsStart, $tmEnsEnd, $tmEnsStrand, $tmDiff) = nearestgene($mmChr,$mPromLoc,\%annsrc);
									($mEnsGene, $mEns, $mEnsStart, $mEnsEnd, $mEnsStrand, $mDiff) = 
										($tmEnsGene, $tmEns, $tmEnsStart, $tmEnsEnd, $tmEnsStrand, $tmDiff)
											if !defined ($mEnsGene) || abs($tmDiff) < abs($mDiff);

								}
							}
							$mLoc = join(",",@mPromLoc);
							$mStrand = join(",",@mPromStrand);
							$mMatch = join(",",@mPromMatch);
							$mAdptLoc = join(",",@mAdptLoc);
							$mAdptSeq = join(",",@mAdptSeq);
						}
					}
				}
				@tpos = findProm($tseq);
				print "findProm ".scalar(@tpos) if $debug & DEB_INFO;
				if (@tpos == 0) {
					++$stat{'NOT_FOUND'};
				} else {
					if (@tpos == 1) {
						++$stat{'FOUND_Promoter'};
					} else {
						++$stat{'FOUND_Multi_Promoter'};
					}
					foreach $pos (@tpos) {
						($PromLoc, $PromStrand, $PromMatch) = (@{$pos});
						print "findAdpt PromLoc $PromLoc PromStrand $PromStrand PromMatch '$PromMatch'..." if $debug & DEB_INFO;
						($AdptLoc,$AdptSeq) = (undef,undef);
						($AdptLoc,$AdptSeq) = findAdpt($tseq,$PromLoc);
						if (defined $AdptLoc) {
							print "findAdpt rtn $AdptLoc" if $debug & DEB_INFO;
							$AdptLoc += ($hasCoord ? ($rsstart-1) : 0);
						}
						$PromLoc += ($PromStrand != -1 ? $rsstart : $rsend)-1 if $hasCoord;
						($EnsGene, $Ens, $EnsStart, $EnsEnd, $EnsStrand, $Diff) = nearestgene($tchr,$PromLoc, $ann) if $opt_a;

						# print the data...
						@cols = ($query_name,$genomes[$sidx]->{'species'},$mmChr,$mmStart,$mmEnd,$mmStrand);
						push @cols,$mLoc,$mStrand,$mMatch,$mAdptLoc,$mAdptSeq if $opt_x;
						push @cols,$mEnsGene, $mEns, $mEnsStart, $mEnsEnd, $mEnsStrand, $mDiff if $opt_a && $opt_x;
						push @cols,$genomes[$tidx]->{'species'},$tchr,$PromLoc,strandSymbol($PromStrand),$PromMatch,$AdptLoc,$AdptSeq;
						push @cols,$EnsGene, $Ens, $EnsStart, $EnsEnd, $EnsStrand, $Diff if $opt_a;
						print "printing....",@cols;
						print $out @cols;
					}
				}
			}
		} else {
			print STDERR "\t\tmissing hit" if $debug & DEB_HSP;
		}
	}
}
foreach $stat (sort keys %stat) {
	printf "%20s\t%5d\n",$stat,$stat{$stat};
}
exit(0);

sub nearestgene {
	my ($chr, $loc,$ann) = @_;
	return (undef,undef) unless defined $ann->{$chr} && ref ($ann->{$chr}) eq 'ARRAY' && scalar(@{$ann->{$chr}});

	my $lastann = undef;
	my $lastdiff = 0;
	my $annref;
	foreach $annref (@{$ann->{$chr}}) {
		my $isFwd = $annref->[2] eq '+';
		$diff = ($isFwd ? $annref->[0] : $annref->[1]) - $loc;
		unless (defined $lastann) {
			$lastann = $annref;
			$lastdiff = ($isFwd ? 1 : -1) * $diff;
			next;
		}
		if (abs($diff) < abs($lastdiff)) {
			$lastann = $annref;
			$lastdiff = ($isFwd ? 1 : -1) * $diff;
			next;
		} elsif (abs($diff) == abs($lastdiff)) {
			next;
		}
		#last;
	}
	return (undef,undef,undef,undef,undef,undef) if abs($lastdiff) > $annLimit;
	#return ($2,$1,$lastann->[0],$lastann->[1],$lastann->[2],$lastdiff) if $lastann->[3] =~ m/gene=(\S+);group=(\S+);*/;
	return ($2,$1,$lastann->[0],$lastann->[1],$lastann->[2],$lastdiff) if $lastann->[3] =~ m/ID=(\S+);Name=(\S+);*/;
	return ($1,$1,$lastann->[0],$lastann->[1],$lastann->[2],$lastdiff) if $lastann->[3] =~ m/ID=(\S+);*/;
	return ("$chr:$lastann->[0]-$lastann->[1]",undef,$lastann->[0],$lastann->[1],$lastann->[2],$lastdiff);
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
	print"\t-h\tThis help message\n";
	foreach $o (@otherOptions) {
		print "\t$o->[0]\t$o->[1]\n";
	}
	exit(1);
}
sub readgff {
	my ($annfile, $ann) = @_;

	$anncnt = 0;
	open GFF,$annfile or die "Cant open GFF file '$annfile' : $!";
	# read the GFF file 
	while(<GFF>) {
		next if m/^#/;
		next if m/^\s*$/;
		chomp;
		(@flds) = split /\t/,$_,-1;
		next if $flds[2] ne 'gene';	# only want gene
		$ann->{$flds[0]} = () unless exists $ann->{$flds[0]};
		push @{$ann->{$flds[0]}}, [$flds[3],$flds[4],$flds[6],$flds[8]];
		++$anncnt;
	}
	print STDERR "$annfile: read $anncnt GFF records...";
	foreach $key (sort chrtrie keys %{$ann}) {
		print STDERR "$key ";
		@tmp = sort {$a->[0] <=> $b->[0]} @{$ann->{$key}};
		$ann->{$key} = ();
		push @{$ann->{$key}},@tmp;
	}
	print STDERR "\n";
	close GFF;
}
sub chrtrie {
	return $a <=> $b if ($a =~ m/^\d+$/ && $b =~ m/^\d+$/);
	return $a cmp $b;
}
sub indexmm {
#
# find case_insensitive str in src allowing for mm mismatches and starting from pos
#
	my ($str,$src,$mm,$pos) = @_;
	return 0 unless length($str);
	return -1 unless length($src);
	$str = uc $str;
	$src = uc $src;
	$mm = 0 unless defined $mm || $mm >= 0;
	$pos = 0 unless defined $pos || $pos >= 0;
	return index($src,$str,$pos) if $mm == 0;
	$lsrc = length($src);
	$lstr = length($str);
	for (my $i = $pos; $i < $lsrc-$lstr+1;++$i) {
		my $it = $i;
		my $nmm = 0;
		my $fnd = 1;
		#print "$i: ";
		for (my $j = 0; $j < $lstr; ++$j,++$it) {
			if ($it > $lsrc) {
				$fnd = 0;
				last;
			}
			#print "[$j $it]";
			if (substr($str,$j,1) ne substr($src,$it,1)) {
				if (++$nmm > $mm) {
					$fnd = 0;
					last;
				}
			}else {
				#print substr($src,$it,1);
			}
		}
		return $i if $fnd;
		#print "\n";
	}
	return -1;
}
sub findProm {
	my $seq = shift;
	die "findProm must be passed a sequence\n" unless length($seq);
	my @pos = ();
	if ($PromoterUseMM) {
		my $lpos = 0;
		while(1) {
			my $pos = indexmm $Promoter, $seq, $Promotermismatch,$lpos;
			last if $pos < 0;
			$match = substr ($seq,$pos,$Promoterlen);
			push @pos,[$pos,1,$match];
			$lpos = $pos+1;
		}
		$lpos = 0;
		while(1) {
			my $pos=indexmm $Promoterrev,$seq,$Promotermismatch,$lpos;
			last if $pos < 0;
			my $match = substr ($seq,$pos,$Promoterlen);
			$match = reverse $match;
			$match =~ tr/AaTtGgCc/TtAaCcGg/;
			$pos += ($PromoterLen-1);
			push @pos,[$pos,-1,$match];
			$lpos = $pos+1;
		}
	} else {	# using an array of perl RE
		foreach my $gli3 (@Promoter) {
			while($seq =~ m/$gli3/ig) {
				my $match = $&;
				push @pos,[pos($seq)-length($&),1,$match];
			}
		}
		foreach $gli3rev (@Promoterrev) {
			while($seq =~ m/$gli3rev/ig) {
				my $match = reverse $&;
				$match =~ tr/AaTtGgCc/TtAaCcGg/;
				push @pos,[pos($seq)-1,-1,$match];
			}
		}
	}
	return @pos;
}
sub findAdpt {
	my ($seq, $loc) = @_;

	die "findAdpt must be passed a sequence and location\n" unless length($seq) &&  defined ($loc) && $loc >= 0;
	my $bestloc = undef;
	my $bestseq = undef;
	my $diff = undef;
	my $mpos = 0;
	my $msx;
	foreach $msx (@Adpt,@Adptrev) {
		$mpos = 0;
		while(1) {
			$mpos = indexmm $msx, $seq, $Adptmismatch,$mpos;
			if ($mpos >= 0) {
				if ((! defined($diff)) || abs($mpos-$loc) < $diff) {
					$diff = abs($mpos-$loc);
					$bestloc = $mpos;
					$bestseq = $msx;
				}
				++$mpos;
			} else {
				last;
			}
		}
	}
	return ($bestloc,$bestseq) if defined($diff) && $diff <= $PromoterAdptDist;
	return (undef,undef);
}
sub strandSymbol {
	my ($val) = @_;

	return '+' if $val == 1 || $val eq '+';
	return '-' if $val == -1 || $val == 0 || $val eq '-';
	print "strandSymbol bad value '$val'\n";
	return $val;
}
