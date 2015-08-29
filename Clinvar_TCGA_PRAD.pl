#!/usr/bin/perl -w
#use strict;
#use warnings;
#use diagnostics;


my $Variants		= "/lustre/beagle2/preed/ClinVar_Reference_Final.txt";
my $TCGA_List		= "/lustre/beagle2/preed/File_Paths_PRAD.txt";
my $Output_File1	= "/lustre/beagle2/preed/TCGA_PRAD_Freq.txt";
my $Output_File2	= "/lustre/beagle2/preed/TCGA_PRAD_Samp.txt";
my (@tmp1,@FILES);
my (%DATA1,%SAMPLES);
my $chr = "chr";

##Parse ClinVar VCF File
open (LIST, $Variants) || die "File not found\n";
     while (<LIST>) {
    chomp;
        @tmp1 = split(/\t/, $_);      
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Chr"} = "$chr$tmp1[0]";
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Pos"} = $tmp1[1];
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Ref"} = $tmp1[3];
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Alt"} = $tmp1[4];
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"TCGA_Samples"} = 0;
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"TCGA_GT_Hets"} = 0;
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"TCGA_GT_Homs"} = 0;		
};
close(LIST);
print "1\n";

#############################################################################################################################################################
### TCGA Data ####

open (LIST, $TCGA_List) || die "File not found.1.\n";
	@FILES = <LIST>;
close(LIST);
 
foreach my $file (@FILES) {
chomp($file);
print "parsing file $file\n";

	my @Sample_data = split(/\//,$file);
    my $Sample = $Sample_data[(scalar(@Sample_data) - 1)];
     
open (LIST, $file) || die "File not found.3.\n";
     while (<LIST>) {
     chomp; 
         @tmp1 = split(/\t/, $_);
    if ($_ =~ /GT:GL:GOF:GQ:NR:NV/) {    
         my @alts = split(/\,/,$tmp1[4]);
         my @DATA = split(/\:/, $tmp1[9]);
         my @GT = split(/\,/, $DATA[0]);
         my @NR = split(/\,/, $DATA[4]);
         my @NV = split(/\,/, $DATA[5]);
        foreach my $i (0 .. $#alts) {       
        if ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}){
        if (($tmp1[6] eq "PASS") and ($NR[$i] >= 10)) {	
        push (@{$SAMPLES{"$Sample"}}, "$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]");
		$DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"TCGA_Samples"}++;
        if ($GT[$i] eq "0/1") {	 	
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"TCGA_GT_Hets"}++;
        push (@{$DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"TCGA_Het_Samples"}},"$Sample");
         	}
        elsif ($GT[$i] eq "1/1") {
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"TCGA_GT_Homs"}++;
        push (@{$DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"TCGA_Hom_Samples"}},"$Sample");
					}
				}		
			}
		}
	}
};		
close(LIST);
};

print "2\n";
##########  OUTPUT ###########################################################################################################################################
@VARs =keys(%DATA1);

open(MYOUTPUTFILE, ">$Output_File1");
print MYOUTPUTFILE "Chr\tPos\tRef\tAlt\tTCGA_Samples\tTCGA_GT_Hets\tTCGA_Het_Samples\tTCGA_GT_Homs\tTCGA_Hom_Samples\n";

foreach (@VARs) {
print MYOUTPUTFILE $DATA1{$_}{"Chr"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"Pos"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"Ref"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"Alt"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"TCGA_Samples"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"TCGA_GT_Hets"}."\t";

if ($DATA1{$_}{"TCGA_Het_Samples"}) {
my @array = @{$DATA1{$_}{"TCGA_Het_Samples"}};
foreach my $var (@array) {
print MYOUTPUTFILE $var.",";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

print MYOUTPUTFILE $DATA1{$_}{"TCGA_GT_Homs"}."\t";

if ($DATA1{$_}{"TCGA_Hom_Samples"}) {
my @array = @{$DATA1{$_}{"TCGA_Hom_Samples"}};
foreach my $var (@array) {
print MYOUTPUTFILE $var.",";
}
print MYOUTPUTFILE "\n";
}
else {print MYOUTPUTFILE "\n"};
};
close(MYOUTPUTFILE);



@VARs =keys(%SAMPLES);

open(MYOUTPUTFILE, ">$Output_File2");
print MYOUTPUTFILE "Sample\tVariants\n";
foreach (@VARs) {
print MYOUTPUTFILE $_."\t";
my @array = @{$SAMPLES{"$_"}};
foreach my $var (@array) {
print MYOUTPUTFILE $var.",";
}
print MYOUTPUTFILE "\n";
};

close(MYOUTPUTFILE);

