#! /local/gensoft/adm/bin/perl -w


use strict;
my $current_id = "";
my $current_nb = 0;
my $previous_id = "";
my $previous_nb = 0;
my ($current_seq, $previous_seq) = ("","");
my $current_full = "";
my $previous_full = "";

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $se_outfile = $outfile . ".se.fasta";
my $pe_reads1_outfile = $outfile . ".pe_reads1.fasta";
my $pe_reads2_outfile = $outfile. ".pe_reads2.fasta";


open (F, $infile) or die;

open (SE, ">$se_outfile") or die;
open (PE1, ">$pe_reads1_outfile") or die;
open (PE2, ">$pe_reads2_outfile") or die;


while (<F>){
    chomp;

    if (m/\/1$/){
	$current_full = $_;
	($current_id, $current_nb) = m/(.*)\/([12]+)/;
	$current_seq = <F>;
	chomp($current_seq);
	
	if (($previous_nb == 1) || (($previous_nb == 2) && ($previous_id ne $current_id))){
	    print SE  $previous_full . "\n" . $previous_seq . "\n";
	}
    }

    else{
	if (m/\/2$/){
	    $current_full = $_;
	    ($current_id, $current_nb) = m/(.*)\/([12]+)/;
	    $current_seq= <F>;
	    chomp($current_seq);

	    if ($previous_nb == 1){
		if ($current_id eq $previous_id){
		    print PE1  $previous_full . "\n" . $previous_seq . "\n";
		    print PE2  $current_full . "\n" . $current_seq . "\n";
		    $current_nb = 0;
		}
		else{
		    print SE  $previous_full . "\n" . $previous_seq . "\n";
		}
	    }
	    else{
		if ($previous_nb == 2){
		    print SE  $previous_full . "\n" . $previous_seq . "\n";
		}
	    }
	    
	}

    }

    $previous_id = $current_id;
    $previous_nb = $current_nb;
    $previous_seq = $current_seq;
    $previous_full = $current_full;
}
if ($previous_nb != 0){
    print SE  $previous_full . "\n" . $previous_seq . "\n";
}
close(F);
close(SE);
close(PE1);
close(PE2);
