#!/usr/bin/perl -w
#use strict;
#use warnings;
#use diagnostics;


my $Variants		= "./ClinVar_Reference_Final.txt";
my $NDAR_DATA		= "./wigler-ndar-variants.PASS.vcf";
my $NDAR_HEAD		= "./wigler-ndar-variants.header";
my $Output_File1	= "./NDAR_wigler_Freq.txt";
my $Output_File2	= "./NDAR_wigler_Samp.txt";
my (@tmp1,%DATA1,%SAMPLES,@Sample);
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
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"NDAR_Samples"} = 0;
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"NDAR_GT_Hets"} = 0;
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"NDAR_GT_Homs"} = 0;		
};
close(LIST);
print "1\n";

#############################################################################################################################################################
### NDAR Data ####
open (LIST, $NDAR_HEAD) || die "File not found.2.\n";     
    while (<LIST>) {
    chomp;    
    @tmp1 = split(/\t/, $_);
	my $size = scalar(@tmp1);
    my $index = $size - 1;
    if ($_ =~ /#CHROM/) { @Sample = @tmp1[9..$index];}
    };
close(LIST);
print "2\n";

open (LIST, $NDAR_DATA) || die "File not found.2.\n";     
    while (<LIST>) {
    chomp;    
    @tmp1 = split(/\t/, $_);
    my @alts = split(/\,/,$tmp1[4]);
    my $size = scalar(@tmp1);
    my $index = $size - 1;
        if ($_ =~ /GT:AD:DP:GQ:PL/) {
        foreach my $i (0 .. $#alts) {
        	if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}){
         	my @Genotypes = @tmp1[9..$index];
         	my $GT_Het = 0;
         	my $GT_Hom = 0;
         	my $position = 0;
         	foreach my $GT (@Genotypes) {
         		if ($GT ne ".") {
         		my @data = split(':',$GT);
         			if (($tmp1[6] eq "PASS") and ($data[2] >= 10)) {
         			push (@{$SAMPLES{"$Sample[$position]"}}, "$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]");
					$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"NDAR_Samples"}++;
         				if ($data[0] eq "0/1") {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"NDAR_GT_Hets"}++;
         				push (@{$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"NDAR_Het_Samples"}},"$Sample[$position]");
         				}
         				elsif ($data[0] eq "1/1") {$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"NDAR_GT_Homs"}++;
         				push (@{$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"NDAR_Hom_Samples"}},"$Sample[$position]");	
         					}
         				}
         			}
         		$position++;
         		}						 	
			}
		}
	}
};	
close(LIST);
print "3\n";
##########  OUTPUT ###########################################################################################################################################
@VARs =keys(%DATA1);

open(MYOUTPUTFILE, ">$Output_File1");
print MYOUTPUTFILE "Chr\tPos\tRef\tAlt\tNDAR_Samples\tNDAR_GT_Hets\tNDAR_Het_Samples\tNDAR_GT_Homs\tNDAR_Hom_Samples\n";

foreach (@VARs) {
print MYOUTPUTFILE $DATA1{$_}{"Chr"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"Pos"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"Ref"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"Alt"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"NDAR_Samples"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"NDAR_GT_Hets"}."\t";

if ($DATA1{$_}{"NDAR_Het_Samples"}) {
my @array = @{$DATA1{$_}{"NDAR_Het_Samples"}};
foreach my $var (@array) {
print MYOUTPUTFILE $var.",";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

print MYOUTPUTFILE $DATA1{$_}{"NDAR_GT_Homs"}."\t";

if ($DATA1{$_}{"NDAR_Hom_Samples"}) {
my @array = @{$DATA1{$_}{"NDAR_Hom_Samples"}};
foreach my $var (@array) {
print MYOUTPUTFILE $var.",";
}
print MYOUTPUTFILE "\n";
}
else {print MYOUTPUTFILE "\n"};
}
close(MYOUTPUTFILE);


@VARs =keys(%SAMPLES);
open(MYOUTPUTFILE, ">$Output_File2");
print MYOUTPUTFILE "Sample\tVariants\n";
foreach (@VARs) {
print MYOUTPUTFILE $_."\t";
my @array = @{$SAMPLES{$_}};
foreach my $var (@array) {
print MYOUTPUTFILE $var.",";
}
print MYOUTPUTFILE "\n";
};
close(MYOUTPUTFILE);

