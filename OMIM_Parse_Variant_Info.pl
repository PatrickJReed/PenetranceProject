#!/usr/bin/perl -w
#use strict;
#use warnings;
#use diagnostics;
use lib qw(..);
use JSON qw( );

my $Variants		= "./ClinVar_Reference_Final.txt";
my $OMIM_List	= "./OMIM_Files.txt";
my $Output_File	= "./OMIM_CLNSRCID_Data.txt";
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
        @tmp2=split(/\=/,$tmp1[10]);
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"CLNSRCID"} = $tmp2[1];
};
close(LIST);
#print "1\n";

#############################################################################################################################################################
### OMIM Data ####
my %Variants;
open (LIST, $OMIM_List) || die "File not found.1.\n";
	@FILES = <LIST>;
close(LIST);

foreach my $file (@FILES) {
chomp($file);
my @Allele_Data = split(/[=&]+/,$file);
my $Sample = $Allele_Data[1];
#print $Sample."\n";
my $json_text = do {
   open(my $json_fh, "<:encoding(UTF-8)", $file)
      or die("Can't open \$file\": $!\n");
   local $/;
   <$json_fh>
};
my $AD = "autosomal dominant";
my $AR = "autosomal recessive";
my $HOM = "homozyg";
my $CH  = "compound heterozyg";
my $HET = "heterozyg";
my $DN  = "de novo";
my $CS  = " consanguin";
my $SP  = "sporadic";
my $XL  = "x-linked";
my $RF  = "risk factor";

my $json = JSON->new;
my $data = $json->decode($json_text);
my @variant_count = keys($data->{omim}->{entryList}[0]->{entry}->{allelicVariantList});
foreach my $var (@variant_count) {
my $Variant_Data = $data->{omim}->{entryList}[0]->{entry}->{allelicVariantList}[$var]->{allelicVariant}->{text};
 $Variant_Data =~ s/\R//g;
$Variants{$Sample}{$var}{"INFO"} = $Variant_Data;
$Variants{$Sample}{$var}{"Inheritance"} = "-";

if (index(lc($Variant_Data), lc($RF)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "RF";
}
if (index(lc($Variant_Data), lc($XL)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "XL";
}
if (index(lc($Variant_Data), lc($SP)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "AD";
}
if (index(lc($Variant_Data), lc($CS)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "AR";
}
if (index(lc($Variant_Data), lc($DN)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "AD";
}
if (index(lc($Variant_Data), lc($HET)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "AD";
}
if (index(lc($Variant_Data), lc($CH)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "CH";
}
if (index(lc($Variant_Data), lc($HOM)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "AR";
}
if (index(lc($Variant_Data), lc($AR)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "AR";
}
if (index(lc($Variant_Data), lc($AD)) != -1) {
    $Variants{$Sample}{$var}{"Inheritance"} = "AD";
    }
  };
};

my @VARs = keys(%DATA1);
foreach (@VARs) {
  my @ID = split(/[\.|]+/,$DATA1{$_}{"CLNSRCID"});
 #print $DATA1{$_}{"CLNSRCID"}."\n";
 #print "@ID"."\n";
  $DATA1{$_}{"OMIM_Variant_Info"} = $Variants{$ID[0]}{($ID[1]+=0)-1}{"INFO"};
  $DATA1{$_}{"OMIM_Inheritance"} = $Variants{$ID[0]}{($ID[1]+=0)-1}{"Inheritance"};
};

@VARs =keys(%DATA1);

open(MYOUTPUTFILE, ">$Output_File");
print MYOUTPUTFILE "CLNSRCID\tInheritance\tVariant_Info\n";
foreach (@VARs) {
print MYOUTPUTFILE $DATA1{$_}{"CLNSRCID"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"OMIM_Inheritance"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"OMIM_Variant_Info"}."\n";
};

close(MYOUTPUTFILE);
