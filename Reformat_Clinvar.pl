#!/usr/bin/perl -w
#use strict;
#use warnings;
#use diagnostics;
use XML::Simple;
use Data::Dumper;
use feature 'say';
my $Variants		= "./clinvar.vcf";
my $Output_File		= "./clinvar_tmp.vcf";

my (@tmp1,%DATA1,@CLNALLE,@CLNSRC,@CLNORIGIN,@CLNSRCID,@CLNSIG,@CLNDSDB,@CLNDSDBID,@CLNDBN,$GENEINFO,$VC);


##Parse ClinVar VCF File
open (LIST, $Variants) || die "File not found\n";
     while (<LIST>) {
    chomp;
        @tmp1 = split(/\t/, $_);
        my @alts = split(/\,/,$tmp1[4]);
        my @DATA = split(/\;/, $tmp1[7]);
         foreach my $datum (@DATA) {
         my @current_datum = split(/\=/, $datum);
         if ($current_datum[0] eq "CLNALLE") {
         @CLNALLE = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "CLNSRC") {
         @CLNSRC = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "CLNORIGIN") {
         @CLNORIGIN = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "CLNSRCID") {
         @CLNSRCID = split(/\,/,$current_datum[1]);
        		}
        elsif  ($current_datum[0] eq "CLNSIG") {
         @CLNSIG = split(/\,/,$current_datum[1]);
        		}
        elsif  ($current_datum[0] eq "CLNDSDB") {
         @CLNDSDB = split(/\,/,$current_datum[1]);
        		}
        elsif  ($current_datum[0] eq "CLNDSDBID") {
         @CLNDSDBID = split(/\,/,$current_datum[1]);
        		}
        elsif  ($current_datum[0] eq "CLNDBN") {
         @CLNDBN = split(/\,/,$current_datum[1]);
        		}
        elsif  ($current_datum[0] eq "CLNHGVS") {
         @CLNHGVS = split(/\,/,$current_datum[1]);
        		}		
        elsif  ($current_datum[0] eq "GENEINFO") {
         $GENEINFO = $current_datum[1];
        		}
        elsif  ($current_datum[0] eq "VC") {
         $VC = $current_datum[1];
        		}												
        	}
        if (scalar(@alts) eq scalar(@CLNSRCID)) {
        foreach my $i (0 .. $#alts) {
		$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"Chr"} = $tmp1[0];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"Pos"} = $tmp1[1];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"dbSNP"} = $tmp1[2];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"Ref"} = $tmp1[3];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"Alt"} = $alts[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"GENEINFO"} = $GENEINFO;
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"VC"} = $VC;
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNALLE"} = $CLNALLE[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSRC"} = $CLNSRC[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNORIGIN"} = $CLNORIGIN[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSRCID"} = $CLNSRCID[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSIG"} = $CLNSIG[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDSDB"} = $CLNDSDB[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDSDBID"} = $CLNDSDBID[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDBN"} = $CLNDBN[$i];
	}
}	
       elsif((scalar(@CLNSRCID) eq 1) and ($CLNSRCID[0] ne ".")) {
       my @Possible_alts = split('>',$CLNHGVS[0]);
       if (scalar(@Possible_alts) > 1) {
       my $correct_alt = $Possible_alts[1];
       print $correct_alt."\n";
       foreach my $i (0 .. $#alts) {
       if ($alts[$i] eq $correct_alt) {
		$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"Chr"} = $tmp1[0];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"Pos"} = $tmp1[1];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"dbSNP"} = $tmp1[2];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"Ref"} = $tmp1[3];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"Alt"} = $alts[$i];
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"GENEINFO"} = $GENEINFO;
        $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"VC"} = $VC;
        if ($CLNALLE[$i]) {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNALLE"} = $CLNALLE[$i];} else {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNALLE"} = $CLNALLE[0];}
        if ($CLNSRC[$i]) {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSRC"} = $CLNSRC[$i];} else {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSRC"} = $CLNSRC[0];}
        if ($CLNORIGIN[$i]) {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNORIGIN"} = $CLNORIGIN[$i];} else {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNORIGIN"} = $CLNORIGIN[0];}
        if ($CLNSRCID[$i]) {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSRCID"} = $CLNSRCID[$i];} else {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSRCID"} = $CLNSRCID[0];}
        if ($CLNSIG[$i]) {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSIG"} = $CLNSIG[$i];} else {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNSIG"} = $CLNSIG[0];}
        if ($CLNDSDB[$i]) {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDSDB"} = $CLNDSDB[$i];} else {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDSDB"} = $CLNDSDB[0];}
        if ($CLNDSDBID[$i]) {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDSDBID"} = $CLNDSDBID[$i];} else {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDSDBID"} = $CLNDSDBID[0];}
        if ($CLNDBN[$i]) {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDBN"} = $CLNDBN[$i];} else {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"CLNDBN"} = $CLNDBN[0];}
				}
			}
		}
	}              	
};	
close(LIST);


##########  OUTPUT ###########################################################################################################################################
@VARs =keys(%DATA1);

open(MYOUTPUTFILE, ">$Output_File");
print MYOUTPUTFILE "Chr\tPos\tdbSNP\tRef\tAlt\tGENEINFO\tVC\tCLNALLE\tCLNSRC\tCLNORIGIN\tCLNSRCID\tCLNSIG\tCLNDSDB\tCLNDSDBID\tCLNDBN\n";

foreach (@VARs) {
if ($DATA1{$_}{"Chr"}) {print MYOUTPUTFILE $DATA1{$_}{"Chr"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Pos"}) {print MYOUTPUTFILE $DATA1{$_}{"Pos"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"dbSNP"}) {print MYOUTPUTFILE $DATA1{$_}{"dbSNP"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Ref"}) {print MYOUTPUTFILE $DATA1{$_}{"Ref"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Alt"}) {print MYOUTPUTFILE $DATA1{$_}{"Alt"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"GENEINFO"}) {print MYOUTPUTFILE "GENEINFO=".$DATA1{$_}{"GENEINFO"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"VC"}) {print MYOUTPUTFILE "VC=".$DATA1{$_}{"VC"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNALLE"}) {print MYOUTPUTFILE "CLNALLE=".$DATA1{$_}{"CLNALLE"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNSRC"}) {print MYOUTPUTFILE "CLNSRC=".$DATA1{$_}{"CLNSRC"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNORIGIN"}) {print MYOUTPUTFILE "CLNORIGIN=".$DATA1{$_}{"CLNORIGIN"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNSRCID"}) {print MYOUTPUTFILE "CLNSRCID=".$DATA1{$_}{"CLNSRCID"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNSIG"}) {print MYOUTPUTFILE "CLNSIG=".$DATA1{$_}{"CLNSIG"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNDSDB"}) {print MYOUTPUTFILE "CLNDSDB=".$DATA1{$_}{"CLNDSDB"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNDSDBID"}) {print MYOUTPUTFILE "CLNDSDBID=".$DATA1{$_}{"CLNDSDBID"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNDBN"}) {print MYOUTPUTFILE "CLNDBN=".$DATA1{$_}{"CLNDBN"}."\n";} else {print MYOUTPUTFILE "\n";}
}; 	
close(MYOUTPUTFILE);

#`grep 'CLNSRC=OMIM_Allelic_Variant' clinvar_tmp.vcf | grep 'CLNSIG=5\t' | grep 'CLNORIGIN=1\|CLNORIGIN=4\|CLNORIGIN=8\|CLNORIGIN=16\|CLNORIGIN=32\|CLNORIGIN=64\|CLNORIGIN=128' | grep 'MedGen:OMIM:Orphanet' > ClinVar_Reference_Final.txt`;

#`rm clinvar_tmp.vcf`;

