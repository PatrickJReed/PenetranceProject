#!/usr/bin/perl -w
#use strict;
#use warnings;
#use diagnostics;
use XML::Simple;
use Data::Dumper;
use List::Util 'sum';
my $Variants		= "./ClinVar_Reference_Final.txt";
my $EXAC_DATA_All		= "./ExAC.r0.3.sites.PASS.vcf";
my $EXAC_DATA_nonTCGA		= "./ExAC.r0.3.nonTCGA.sites.vep.vcf";
my $EXAC_DATA_nonPSYCH		= "./ExAC.r0.3.nonpsych.sites.vcf";
my $KAV_DATA    = "./Kaviar-150810-Public-hg19.vcf";
my $HRC_DATA    = "./HRC.r1.GRCh37.autosomes.mac5.sites.vcf";
my $MedGen_OMIM		= "./MedGen_HPO_OMIM_Mapping.txt";
my $MG_ICD		    = "./MRCONSO_ENG_ICD10CM.RRF";
my $MG_REL			= "./MRREL.RRF";
my $OMIM_CLNSRCID_Data = "./OMIM_CLNSRCID_Data.txt";
my $Population_Data = "./Population_rawdata_2119.txt";
my $Output_File		= "./Mendelian_Data_Final.txt";
my (@tmp1,@tmp2,@VARs,@Samples,@HPO_CUIs,@RELAs,@RELAs1,@RELAs2,@EXAC_AC,@EXAC_AC_Het,@EXAC_AC_Hom,@EXAC_AF);
my (%DATA1,%DATA2,%DATA3_1,%DATA3_2,%DATA4,%DATA5,%NAMES,%SAMPLES,$NAME_1,$NAME_2,$NAME_3,%ORPHA_ages,%ORPHA_prev,%ORPHA_data,@CLNDSDB,@CLNDSDBID,@Disorders_data,%OMIM_Data);
my $chr = "chr";
my $xml_basic = new XML::Simple;
my $xml_ages = new XML::Simple;
my $xml_prev = new XML::Simple;
my $data_basic = $xml_basic->XMLin("./en_product1.xml");
my $data_ages = $xml_ages->XMLin("./en_product2_ages.xml");
my $data_prev = $xml_prev->XMLin("./en_product2_prev.xml");
my %population;
$population{"<1 / 1 000 000"} = 0.000001;
$population{"1-9 / 1 000 000"} = 0.000005;
$population{"1-9 / 100 000"} = 0.00005;
$population{"6-9 / 10 000"} = 0.00075;
$population{"1-5 / 10 000"} = 0.0003;
$population{">1 / 1000"}  = 0.001;
$population{"Unknown"}  = "-";

####Country/Region Population Data##############################
my %PopulationSize;
open (LIST, $Population_Data) || die "File not found\n";
     while (<LIST>) {
    chomp;
        @tmp1 = split( /\s{2,}/, $_);
        $PopulationSize{$tmp1[1]} = $tmp1[2];
};
close(LIST);


#Orpha_ICD10CM_###########################################################
 @Disorders_data = keys($data_basic->{DisorderList}->{Disorder});
foreach (@Disorders_data) {
my $OrphaNumber = $data_basic->{DisorderList}->{Disorder}->{$_}->{OrphaNumber};
my $Orpha_Name = $data_basic->{DisorderList}->{Disorder}->{$_}->{Name}->{content};
$ORPHA_data{"$OrphaNumber"}{"Orpha_Name"} = $Orpha_Name;
my $externalreferences =0;
if ($data_basic->{DisorderList}->{Disorder}->{$_}->{ExternalReferenceList}->{count}){ $externalreferences = $data_basic->{DisorderList}->{Disorder}->{$_}->{ExternalReferenceList}->{count};}
if ($externalreferences > 1) {
my @test = keys(%{ $data_basic->{DisorderList}->{Disorder}->{$_}->{ExternalReferenceList}->{ExternalReference} });
foreach my $t (@test) {
	if ($data_basic->{DisorderList}->{Disorder}->{$_}->{ExternalReferenceList}->{ExternalReference}{$t}->{Source} eq "ICD-10") {
		my $ICD_Code = $data_basic->{DisorderList}->{Disorder}->{$_}->{ExternalReferenceList}->{ExternalReference}{$t}->{Reference};
		push (@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}}, "$ICD_Code");
		@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}});
			}
		else {push (@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}}, "");
				}
			}
		}
elsif ($externalreferences == 1) {
	if ($data_basic->{DisorderList}->{Disorder}->{$_}->{ExternalReferenceList}->{ExternalReference}->{Source} eq "ICD-10") {
		my $ICD_Code = $data_basic->{DisorderList}->{Disorder}->{$_}->{ExternalReferenceList}->{ExternalReference}->{Reference};
		push (@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}}, "$ICD_Code");
		@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}});
			}
		else {push (@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}}, "");
		}
	}
elsif ($externalreferences == 0) {
		my $ICD_Code = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}}, $ICD_Code);
		@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}});
}

#Orpha_Inheritance###########################################################
my $TypeOfInheritances = 0;
if ($data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{count}) { $TypeOfInheritances = $data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{count};}
if ($TypeOfInheritances > 1) {
my @test = keys(%{ $data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{TypeOfInheritance} });
foreach my $t (@test) {
		my $inheritance = $data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{TypeOfInheritance}{$t}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}}, $inheritance);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}});
				}
			}
elsif ($TypeOfInheritances == 1) {
		my $inheritance = $data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{TypeOfInheritance}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}}, $inheritance);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}});
	}
elsif ($TypeOfInheritances == 0) {
		my $inheritance = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}}, $inheritance);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}});
	}

#Orpha_Age_Onset###########################################################
my $AgeOfOnsets = 0;
if ($data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{count}) { $AgeOfOnsets = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{count};}
if ($AgeOfOnsets > 1) {
my @test = keys(%{ $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{AverageAgeOfOnset} });
foreach my $t (@test) {
		my $onset = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{AverageAgeOfOnset}{$t}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}}, $onset);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}});
				}
			}
elsif ($AgeOfOnsets == 1) {
		my $onset = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{AverageAgeOfOnset}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}}, $onset);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}});
	}
elsif ($AgeOfOnsets == 0) {
		my $onset = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}}, $onset);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}});
	}
#Orpha_Age_Death###########################################################
my $AgeOfDeaths = 0;
if ($data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{count}) { $AgeOfDeaths = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{count};}
if ($AgeOfDeaths > 1) {
my @test = keys(%{ $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{AverageAgeOfDeath} });
foreach my $t (@test) {
		my $death = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{AverageAgeOfDeath}{$t}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}}, $death);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}});
				}
			}
elsif ($AgeOfDeaths == 1) {
		my $death = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{AverageAgeOfDeath}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}}, $death);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}});
	}
elsif ($AgeOfDeaths == 0) {
		my $death = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}}, $death);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}});
	}

#Orpha_Disorder_Type###########################################################
my $TypeOfDisorder = $data_ages->{DisorderList}->{Disorder}->{$_}->{DisorderType}->{Name}->{content};
$ORPHA_data{"$OrphaNumber"}{"Orpha_Type"} = $TypeOfDisorder;

#Orpha_Prevalence###########################################################
my $Prevelances = 0;
if ($data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{count}) { $Prevelances = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{count};}
if ($Prevelances > 1) {
my @test = keys(%{$data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence} });
foreach my $t (@test) {
		my $prevalence = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}{$t}->{PrevalenceClass}->{Name}->{content};
		my $type = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}{$t}->{PrevalenceType}->{Name}->{content};
		my $location = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}{$t}->{PrevalenceGeographic}->{Name}->{content};
		if ($data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}{$t}->{PrevalenceClass}->{Name}->{content}) {
		 if ($PopulationSize{"$location"}) {
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}}, "$PopulationSize{$location}:$population{$prevalence}");
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}});
					}
				}
			}
		}
elsif ($Prevelances == 1) {
		my $prevalence = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceClass}->{Name}->{content};
		my $type = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceType}->{Name}->{content};
		my $location = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceGeographic}->{Name}->{content};
		if ($data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceClass}->{Name}->{content}) {
		if ($PopulationSize{"$location"}) {
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}}, "$PopulationSize{$location}:$population{$prevalence}");
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}});
		}
	}
}
elsif ($Prevelances == 0) {
		my $type = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}}, $type);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}});
	}
};

print "1\n";
#############################################################################################################################################################

##Parse ClinVar VCF File
open (LIST, $Variants) || die "File not found\n";
     while (<LIST>) {
    chomp;
        @tmp1 = split(/\t/, $_);
		$DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Chr"} = "$chr$tmp1[0]";
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Pos"} = $tmp1[1];
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"dbSNP"} = $tmp1[2];
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Ref"} = $tmp1[3];
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Alt"} = $tmp1[4];
        @tmp2=split(/\=/,$tmp1[5]);
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"GENEINFO"} = $tmp2[1];
        @tmp2=split(/\=/,$tmp1[10]);
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"CLNSRCID"} = $tmp2[1];
        @tmp2=split(/\=/,$tmp1[14]);
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"CLNDBN"} = $tmp2[1];
        @tmp2=split(/\=/,$tmp1[12]);
        @CLNDSDB = split(/\:/,$tmp2[1]);
        @tmp2=split(/\=/,$tmp1[13]);
        @CLNDSDBID = split(/\:/,$tmp2[1]);
	my $position = 0;
	foreach (@CLNDSDB) {
	if ($CLNDSDB[$position] eq "MedGen") {
	$DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"MedGen"} = $CLNDSDBID[$position];
	}
	elsif ($CLNDSDB[$position] eq "OMIM") {
	$DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"OMIM"} = $CLNDSDBID[$position];
	}
	elsif ($CLNDSDB[$position] eq "Orphanet") {
	my @ORPHA_ID = split('ORPHA',$CLNDSDBID[$position]);
	$DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$tmp1[4]"}{"Orphanet"} = $ORPHA_ID[1];
	}
	$position++;
	}
};
close(LIST);
print "2\n";

################################################################################################################
##Parse omim allelic variant ID data (CLNSRCID) from files
open (LIST, $OMIM_CLNSRCID_Data) || die "File not found.2.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         $OMIM_Data{$tmp1[0]}{"Inheritance"} = $tmp1[1];
         $OMIM_Data{$tmp1[0]}{"INFO"} = $tmp1[2];

       };
close(LIST);



#############################################################################################################################################################
### Parse EXAC_All Data ####


open (LIST, $EXAC_DATA_All) || die "File not found.1.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         my @alts = split(/\,/,$tmp1[4]);
         my @DATA = split(/\|/, $tmp1[7]);
         my @Datum = split(/\;/,$DATA[0]);
         foreach my $datum (@Datum) {
         my @current_datum = split(/\=/, $datum);
         if ($current_datum[0] eq "AC") {
         @EXAC_AC = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AC_Het") {
         @EXAC_AC_Het = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AC_Hom") {
         @EXAC_AC_Hom = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AF") {
         @EXAC_AF = split(/\,/,$current_datum[1]);
        		}
        	}
        foreach my $i (0 .. $#alts) {
        if ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}){
		 $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_Chr"} = "chr$tmp1[0]";
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_Start"} = $tmp1[1];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_Ref"} = $tmp1[3];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_Alt"} = $alts[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_Quality"} = $tmp1[5];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_VQSR"} = $tmp1[6];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_Data"} = $DATA[0];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC_Het"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC_Hom"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AF"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC"} = $EXAC_AC[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC_Het"} = $EXAC_AC_Het[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC_Hom"} = $EXAC_AC_Hom[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AF"} = $EXAC_AF[$i];

        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_GT_Het"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC_Het"};
        unless ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC"} eq 0) {
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_GF_Het"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC_Het"} / ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC"} /  $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AF"});}
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_GT_Hom"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC_Hom"};
       unless ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC"} eq 0) {
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_GF_Hom"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC_Hom"} / ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AC"} /  $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_AF"});
        	}
        }
	}
};
close(LIST);

print "3\n";

### Parse EXAC_DATA_nonTCGA ####


open (LIST, $EXAC_DATA_nonTCGA) || die "File not found.1.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         my @alts = split(/\,/,$tmp1[4]);
         my @DATA = split(/\|/, $tmp1[7]);
         my @Datum = split(/\;/,$DATA[0]);
         foreach my $datum (@Datum) {
         my @current_datum = split(/\=/, $datum);
         if ($current_datum[0] eq "AC") {
         @EXAC_nonTCGA_AC = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AC_Het") {
         @EXAC_nonTCGA_AC_Het = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AC_Hom") {
         @EXAC_nonTCGA_AC_Hom = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AF") {
         @EXAC_nonTCGA_AF = split(/\,/,$current_datum[1]);
        		}
        	}
        foreach my $i (0 .. $#alts) {
        if ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}){
		 $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_Chr"} = "chr$tmp1[0]";
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_Start"} = $tmp1[1];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_Ref"} = $tmp1[3];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_Alt"} = $alts[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_Quality"} = $tmp1[5];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_VQSR"} = $tmp1[6];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_Data"} = $DATA[0];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC_Het"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC_Hom"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AF"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC"} = $EXAC_nonTCGA_AC[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC_Het"} = $EXAC_nonTCGA_AC_Het[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC_Hom"} = $EXAC_nonTCGA_AC_Hom[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AF"} = $EXAC_nonTCGA_AF[$i];

        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_GT_Het"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC_Het"};
        unless ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC"} eq 0) {
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_GF_Het"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC_Het"} / ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC"} /  $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AF"});
      }
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_GT_Hom"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC_Hom"};
       unless ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC"} eq 0) {
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_GF_Hom"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC_Hom"} / ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AC"} /  $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonTCGA_AF"});
        	}
        }
	}
};
close(LIST);

print "3_TCGA\n";

### Parse EXAC_DATA_nonPsych ####


open (LIST, $EXAC_DATA_nonPSYCH) || die "File not found.1.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         my @alts = split(/\,/,$tmp1[4]);
         my @DATA = split(/\|/, $tmp1[7]);
         my @Datum = split(/\;/,$DATA[0]);
         foreach my $datum (@Datum) {
         my @current_datum = split(/\=/, $datum);
         if ($current_datum[0] eq "AC") {
         @EXAC_nonPSYCH_AC = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AC_Het") {
         @EXAC_nonPSYCH_AC_Het = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AC_Hom") {
         @EXAC_nonPSYCH_AC_Hom = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AF") {
         @EXAC_nonPSYCH_AF = split(/\,/,$current_datum[1]);
        		}
        	}
        foreach my $i (0 .. $#alts) {
        if ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}){
		 $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_Chr"} = "chr$tmp1[0]";
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_Start"} = $tmp1[1];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_Ref"} = $tmp1[3];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_Alt"} = $alts[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_Quality"} = $tmp1[5];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_VQSR"} = $tmp1[6];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_Data"} = $DATA[0];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC_Het"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC_Hom"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AF"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC"} = $EXAC_nonPSYCH_AC[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC_Het"} = $EXAC_nonPSYCH_AC_Het[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC_Hom"} = $EXAC_nonPSYCH_AC_Hom[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AF"} = $EXAC_nonPSYCH_AF[$i];

        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_GT_Het"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC_Het"};
        unless ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC"} eq 0) {
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_GF_Het"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC_Het"} / ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC"} /  $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AF"});}
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_GT_Hom"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC_Hom"};
       unless ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC"} eq 0) {
        $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_GF_Hom"} = $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC_Hom"} / ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AC"} /  $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"EXAC_nonPSYCH_AF"});
        	}
        }
	}
};
close(LIST);

print "3_PSYCH\n";

#Parse Kaviar Vairnats
open (LIST, $KAV_DATA) || die "File not found.1.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         my @alts = split(/\,/,$tmp1[4]);
         my @DATA = split(/\|/, $tmp1[7]);
         my @Datum = split(/\;/,$DATA[0]);
         foreach my $datum (@Datum) {
         my @current_datum = split(/\=/, $datum);
         if ($current_datum[0] eq "AC") {
         @KAV_AC = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "DS") {
         @KAV_DS = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AF") {
         @KAV_AF = split(/\,/,$current_datum[1]);
        		}
        	}
        foreach my $i (0 .. $#alts) {
        if ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}){
		     $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_Chr"} = "chr$tmp1[0]";
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_Start"} = $tmp1[1];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_Ref"} = $tmp1[3];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_Alt"} = $alts[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_AC"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_AF"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_AC"} = $KAV_AC[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_AF"} = $KAV_AF[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"KAV_DataSource"} = $KAV_DS[$i];
        }
      }
    };
close(LIST);

print "3.2\n";



#Parse HRC Variants
open (LIST, $HRC_DATA) || die "File not found.1.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         my @alts = split(/\,/,$tmp1[4]);
         my @DATA = split(/\|/, $tmp1[7]);
         my @Datum = split(/\;/,$DATA[0]);
         foreach my $datum (@Datum) {
         my @current_datum = split(/\=/, $datum);
         if ($current_datum[0] eq "AC") {
         @HRC_AC = split(/\,/,$current_datum[1]);
         }
         elsif  ($current_datum[0] eq "AF") {
         @HRC_AF = split(/\,/,$current_datum[1]);
        		}
        	}
        foreach my $i (0 .. $#alts) {
        if ($DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}){
		     $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"HRC_Chr"} = "chr$tmp1[0]";
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"HRC_Start"} = $tmp1[1];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"HRC_Ref"} = $tmp1[3];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"HRC_Alt"} = $alts[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"HRC_AC"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"HRC_AF"} = 0;
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"HRC_AC"} = $HRC_AC[$i];
         $DATA1{"$chr$tmp1[0]\t$tmp1[1]\t$tmp1[3]\t$alts[$i]"}{"HRC_AF"} = $HRC_AF[$i];
        }
      }
    };
close(LIST);

print "3.3\n";

#############################################################################################################################################################
### Parse Exome Freq and Samp Files With Clinical metadata###
my $Unique_Individuals = "./unique_individuals.txt";
my %UnInd;
open (LIST, $Unique_Individuals) || die "File not found.2.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         $UnInd{"$tmp1[4]"}{"SampleID"} = $tmp1[0];
       };
close(LIST);

my $NDAR_Clinical1 = "./ndar_subject1.txt";
open (LIST, $NDAR_Clinical1) || die "File not found.2.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         $UnInd{"$tmp1[2]"}{"SampleID"} = $tmp1[2];
       };
close(LIST);



my $NDAR_Clinical2 = "./ndar_subject2.txt";
open (LIST, $NDAR_Clinical2) || die "File not found.2.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         $UnInd{"$tmp1[2]"}{"SampleID"} = $tmp1[2];
       };
close(LIST);


my $TCGA_Clinical = "./TCGA.June29.2015.clinical.patient.txt";
my %Clinical;
open (LIST, $TCGA_Clinical) || die "File not found.3.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         $Clinical{"$tmp1[2]"}{"Age"} = $tmp1[4]/-365;
         $Clinical{"$tmp1[2]"}{"Gender"} = $tmp1[10];
         $Clinical{"$tmp1[2]"}{"Race"} = $tmp1[19];
         $Clinical{"$tmp1[2]"}{"Disease"} = $tmp1[0];
       };
close(LIST);

my $NDAR_Clinical1 = "./ndar_subject1.txt";

open (LIST, $NDAR_Clinical1) || die "File not found.3.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         $Clinical{"$tmp1[2]"}{"Age"} = $tmp1[5]/12;
         $Clinical{"$tmp1[2]"}{"Gender"} = $tmp1[6];
         $Clinical{"$tmp1[2]"}{"Race"} = $tmp1[7];
         $Clinical{"$tmp1[2]"}{"Disease"} = $tmp1[9];
       };
close(LIST);

my $NDAR_Clinical2 = "./ndar_subject2.txt";
open (LIST, $NDAR_Clinical2) || die "File not found.3.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         $Clinical{"$tmp1[2]"}{"Age"} = $tmp1[5]/12;
         $Clinical{"$tmp1[2]"}{"Gender"} = $tmp1[6];
         $Clinical{"$tmp1[2]"}{"Race"} = $tmp1[7];
         $Clinical{"$tmp1[2]"}{"Disease"} = $tmp1[9];
       };
close(LIST);


my @Freq_Files = ("./NDAR_wigler_Freq.txt", "./TCGA_BLCA_Freq.txt","./TCGA_BRCA_Freq.txt", "./TCGA_COAD_Freq.txt", "./TCGA_DLBC_Freq.txt", "./TCGA_ESCA_Freq.txt", "./TCGA_GBM_Freq.txt", "./TCGA_HNSC_Freq.txt", "./TCGA_KIRC_Freq.txt", "./TCGA_KIRP_Freq.txt", "./TCGA_LGG_Freq.txt", "./TCGA_LIHC_Freq.txt", "./TCGA_LUAD_Freq.txt", "./TCGA_LUSC_Freq.txt", "./TCGA_OVCA_Freq.txt", "./TCGA_PAAD_Freq.txt", "./TCGA_PRAD_Freq.txt", "./TCGA_READ_Freq.txt", "./TCGA_SARC_Freq.txt", "./TCGA_STAD_Freq.txt", "./TCGA_THCA_Freq.txt", "./TCGA_UCEC_Freq.txt", "./TCGA_UVM_Freq.txt");


my @Samp_Files = ("./NDAR_wigler_Samp.txt", "./TCGA_BLCA_Samp.txt", "./TCGA_BRCA_Samp.txt", "./TCGA_COAD_Samp.txt", "./TCGA_DLBC_Samp.txt", "./TCGA_ESCA_Samp.txt", "./TCGA_GBM_Samp.txt", "./TCGA_HNSC_Samp.txt", "./TCGA_KIRC_Samp.txt", "./TCGA_KIRP_Samp.txt", "./TCGA_LGG_Samp.txt", "./TCGA_LIHC_Samp.txt", "./TCGA_LUAD_Samp.txt", "./TCGA_LUSC_Samp.txt", "./TCGA_OVCA_Samp.txt", "./TCGA_PAAD_Samp.txt", "./TCGA_PRAD_Samp.txt", "./TCGA_READ_Samp.txt", "./TCGA_SARC_Samp.txt", "./TCGA_STAD_Samp.txt", "./TCGA_THCA_Samp.txt", "./TCGA_UCEC_Samp.txt", "./TCGA_UVM_Samp.txt");

foreach my $file (@Freq_Files) {
chomp $file;
my @samples;
open (LIST, $file) || die "File not found.4.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
#         push(@{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_AC_Het"}},$tmp1[5]);
#         push(@{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_AC_Hom"}},$tmp1[7]);
         @samples_het = split(',',$tmp1[6]);
         @samples_hom = split(',',$tmp1[8]);
         #@samples = uniq(@samples);
        # @{$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples"}} = @samples;
        # print "@samples\n";
         foreach my $samps (@samples_het) {
        #   print "@samples_new"."\n\n\n\n\n\n";
            my $tmp_ID = $UnInd{$samps}{"SampleID"};
            if ($tmp_ID) {
        #      print $tmp_ID."\n";
          my $tmp_age = "NA";
            my $tmp_gender = "NA";
            my $tmp_race = "NA";
            my $tmp_disease = "NA";
            if ($Clinical{$tmp_ID}{"Age"}) {$tmp_age = $Clinical{$tmp_ID}{"Age"};}
            if ($Clinical{$tmp_ID}{"Gender"}) {$tmp_gender = $Clinical{$tmp_ID}{"Gender"};}
            if ($Clinical{$tmp_ID}{"Race"}) {$tmp_race= $Clinical{$tmp_ID}{"Race"};}
            if ($Clinical{$tmp_ID}{"Disease"}) {$tmp_disease= $Clinical{$tmp_ID}{"Disease"};}
            push(@{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples_Het"}}, $tmp_ID.':'.$tmp_age.':'.$tmp_gender.':'.$tmp_race.':'.$tmp_disease);
            push(@{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples_Het_Ages"}}, $tmp_age);
            @{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples_Het"}} = uniq(@{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples_Het"}});
          }
         }
         foreach my $samps (@samples_hom) {
        #   print "@samples_new"."\n\n\n\n\n\n";
            my $tmp_ID = $UnInd{$samps}{"SampleID"};
            if ($tmp_ID) {
        #      print $tmp_ID."\n";
          my $tmp_age = "NA";
            my $tmp_gender = "NA";
            my $tmp_race = "NA";
            my $tmp_disease = "NA";
            if ($Clinical{$tmp_ID}{"Age"}) {$tmp_age = $Clinical{$tmp_ID}{"Age"};}
            if ($Clinical{$tmp_ID}{"Gender"}) {$tmp_gender = $Clinical{$tmp_ID}{"Gender"};}
            if ($Clinical{$tmp_ID}{"Race"}) {$tmp_race= $Clinical{$tmp_ID}{"Race"};}
            if ($Clinical{$tmp_ID}{"Disease"}) {$tmp_disease= $Clinical{$tmp_ID}{"Disease"};}
            push(@{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples_Hom"}}, $tmp_ID.':'.$tmp_age.':'.$tmp_gender.':'.$tmp_race.':'.$tmp_disease);
            push(@{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples_Hom_Ages"}}, $tmp_age);
            @{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples_Hom"}} = uniq(@{$DATA2{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples_Hom"}});
          }
         }
         #print "@samples_new\n";
         #@samples_new = uniq(@samples_new);
         #push(@{$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples"}},@samples);
         #my @array_tmpxx = @{$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_Samples"}};
         #print $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Clinical_AC"}."\n";
         #print "@array_tmpxx"."\n";
         #@{$SAMPLE{"TCGA_NDAR"}{"$SampleName"}{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]"}{"Samples"}} = @samples_new;
        };
		close(LIST);
	};




# foreach my $file (@Samp_Files) {
# chomp($file);
# my @Sample_data = split(/\_/,$file);
# my $SampleName = "$Sample_data[1]";
# open (LIST, $file) || die "File not found.42.\n";
#      while (<LIST>) {
#      chomp;
#          @tmp1 = split(/\t/, $_);
# 		 $SAMPLE{$SampleName}{"Sample"} = $tmp1[0];
#          $SAMPLE{$SampleName}{"Variants"} = $tmp1[1];
#         };
# 		close(LIST);
# 	};

print "4\n";
#Parse MedGen_OMIM Relations
open (LIST, $MedGen_OMIM) || die "File not found.5.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\|/, $_);
         if ($tmp1[1] ne "MIM_number"){
#		 $DATA3{"$tmp1[1]"}{"OMIM_NAME"} = $tmp1[2];

		 $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"OMIM_CUI"} = $tmp1[0];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"MIM_number"} = $tmp1[1];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"OMIM_name"} = $tmp1[2];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"relationship"} = $tmp1[3];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"HPO_CUI"} = $tmp1[4];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"HPO_ID"} = $tmp1[5];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"HPO_name"} = $tmp1[6];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"MedGen_name"} = $tmp1[7];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"MedGen_source"} = $tmp1[8];
         $DATA3_1{"$tmp1[0]"}{"$tmp1[6]"}{"STY"} = $tmp1[9];

         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"OMIM_CUI"} = $tmp1[0];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"MIM_number"} = $tmp1[1];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"OMIM_name"} = $tmp1[2];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"relationship"} = $tmp1[3];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"HPO_CUI"} = $tmp1[4];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"HPO_ID"} = $tmp1[5];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"HPO_name"} = $tmp1[6];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"MedGen_name"} = $tmp1[7];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"MedGen_source"} = $tmp1[8];
         $DATA3_2{"$tmp1[1]"}{"$tmp1[6]"}{"STY"} = $tmp1[9];
      	}
};
close(LIST);


@VARs = keys(%DATA1);
foreach (@VARs) {
#	if ($DATA1{$_}{"OMIM"} and $DATA3{$DATA1{$_}{"OMIM"}}) {$DATA1{$_}{"OMIM_NAME"} = $DATA3{$DATA1{$_}{"OMIM"}}{"OMIM_NAME"};}
	if ($DATA1{$_}{"OMIM"} and $DATA3_2{$DATA1{$_}{"OMIM"}}) {
	$NAME_2 = $DATA1{$_}{"OMIM"};
	@RELAs = keys %{$DATA3_2{$NAME_2}};
	foreach $NAME_3 (@RELAs) {
		push (@{$DATA1{$_}{"HPO_Phenotypes"}}, $NAME_3);
		@{$DATA1{$_}{"HPO_Phenotypes"}} = uniq(@{$DATA1{$_}{"HPO_Phenotypes"}});
		}
	}
};


##Parse MedGen Ref Databases
open (LIST, $MG_REL) || die "File not found.55.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\|/, $_);
         if ($tmp1[7] eq "manifestation_of") {
		 $DATA4{"$tmp1[0]"}{"$tmp1[4]"}{"CUI1"} = $tmp1[0];
         $DATA4{"$tmp1[0]"}{"$tmp1[4]"}{"CUI2"} = $tmp1[4];
         $DATA4{"$tmp1[0]"}{"$tmp1[4]"}{"RELA"} = $tmp1[7];
	}
};
close(LIST);

open (LIST, $MG_ICD) || die "File not found.6.\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\|/, $_);
         $DATA5{"$tmp1[0]"}{"ICD-10"} = $tmp1[10];
};
close(LIST);

@VARs = keys(%DATA1);
foreach (@VARs) {
	if ($DATA1{$_}{"MedGen"} and $DATA4{$DATA1{$_}{"MedGen"}}) {
	$NAME_1 = $DATA1{$_}{"MedGen"};
	@RELAs = keys %{$DATA4{$NAME_1}};
	foreach $NAME_2 (@RELAs) {
		if ($DATA5{$NAME_2}{"ICD-10"}) {
		push (@{$DATA1{$_}{"MedGen_Relations"}}, $DATA5{$NAME_2}{"ICD-10"});
		@{$DATA1{$_}{"MedGen_Relations"}} = uniq(@{$DATA1{$_}{"MedGen_Relations"}});
			}
		}
	}
};

#convert CUIs to Names with MRCONSO
print "5\n";

@VARs = keys(%DATA1);
foreach (@VARs) {
	if ($DATA1{$_}{"Orphanet"} and $ORPHA_data{$DATA1{$_}{"Orphanet"}}) {
	 $NAME_1 = $DATA1{$_}{"Orphanet"};
	 @{$DATA1{$_}{"Orpha_Inheritance"}} = @{$ORPHA_data{$NAME_1}{"Orpha_Inheritance"}};
	 @{$DATA1{$_}{"Orpha_Onset"}} = @{$ORPHA_data{$NAME_1}{"Orpha_Onset"}};
	 @{$DATA1{$_}{"Orpha_Death"}} = @{$ORPHA_data{$NAME_1}{"Orpha_Death"}};
	 $DATA1{$_}{"Orpha_Type"} = $ORPHA_data{$NAME_1}{"Orpha_Type"};
	 $DATA1{$_}{"Orpha_Name"} = $ORPHA_data{$NAME_1}{"Orpha_Name"};
	 @{$DATA1{$_}{"Orpha_Prev"}} = @{$ORPHA_data{$NAME_1}{"Orpha_Prev"}};
	 @{$DATA1{$_}{"ICD-10_Codes"}} = @{$ORPHA_data{$NAME_1}{"ICD-10 Codes"}};
		}
	};

@VARs = keys(%DATA1);
foreach (@VARs) {
$DATA1{$_}{"OMIM_INFO"} = $OMIM_Data{$DATA1{$_}{"CLNSRCID"}}{"INFO"};
$DATA1{$_}{"OMIM_Inheritance"} = $OMIM_Data{$DATA1{$_}{"CLNSRCID"}}{"Inheritance"};
};


##########  OUTPUT ###########################################################################################################################################
@VARs =keys(%DATA1);
#@Samples = keys(%SAMPLE);
open(MYOUTPUTFILE, ">$Output_File");
print MYOUTPUTFILE "Chr\tPos\tdbSNP\tRef\tAlt\tCLNSRCID\tCLNDBN\tMedGen\tOMIM\tOrphanet\tExAC_AC\tExAC_AF\tEXAC_GT_Het\tEXAC_GF_Het\tEXAC_GT_Hom\tEXAC_GF_Hom\tExAC_nonTCGA_AC\tExAC_nonTCGA_AF\tEXAC_nonTCGA_GT_Het\tEXAC_nonTCGA_GF_Het\tEXAC_nonTCGA_GT_Hom\tEXAC_nonTCGA_GF_Hom\tExAC_nonPSYCH_AC\tExAC_nonPSYCH_AF\tEXAC_nonPSYCH_GT_Het\tEXAC_nonPSYCH_GF_Het\tEXAC_nonPSYCH_GT_Hom\tEXAC_nonPSYCH_GF_Hom\tKAV_AC\tKAV_AF\tKAV_DataSource\tHRC_AC\tHRC_AF\tOrpha_Inheritance\tOrpha_Age_Onset\tOrpha_Age_Death\tOrpha_Disorder_Type\tOrpha_Name\tOrpha_Prevalence\tOMIM_Inheritance\tOMIM_INFO\tICD-10\tHPO_Terms\tICD10CM_SubTerms\tTCGA_NDAR_AlleleCount_Het\tTCGA_NDAR_Het_Ages\tTCGA_NDAR_Samples_Het\tTCGA_NDAR_AlleleCount_Hom\tTCGA_NDAR_Hom_Ages\tTCGA_NDAR_Samples_Hom\n";
# foreach my $samp (@Samples) {
# #print MYOUTPUTFILE $samp."-Sample_Count\t";
# print MYOUTPUTFILE $samp."-GT_Hets\t";
# print MYOUTPUTFILE $samp."-Het_Samples\t";
# print MYOUTPUTFILE $samp."-GT_Homs\t";
# print MYOUTPUTFILE $samp."-Hom_Samples\t";
# };
#print MYOUTPUTFILE "\n";
foreach (@VARs) {
if ($DATA1{$_}{"Chr"}) {print MYOUTPUTFILE $DATA1{$_}{"Chr"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Pos"}) {print MYOUTPUTFILE $DATA1{$_}{"Pos"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"dbSNP"}) {print MYOUTPUTFILE $DATA1{$_}{"dbSNP"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Ref"}) {print MYOUTPUTFILE $DATA1{$_}{"Ref"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Alt"}) {print MYOUTPUTFILE $DATA1{$_}{"Alt"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNSRCID"}) {print MYOUTPUTFILE $DATA1{$_}{"CLNSRCID"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNDBN"}) {print MYOUTPUTFILE $DATA1{$_}{"CLNDBN"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"MedGen"}) {print MYOUTPUTFILE $DATA1{$_}{"MedGen"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"OMIM"}) {print MYOUTPUTFILE $DATA1{$_}{"OMIM"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Orphanet"}) {print MYOUTPUTFILE $DATA1{$_}{"Orphanet"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_AC"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_AC"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_AF"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_AF"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_GT_Het"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_GT_Het"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_GF_Het"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_GF_Het"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_GT_Hom"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_GT_Hom"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_GF_Hom"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_GF_Hom"}."\t";} else {print MYOUTPUTFILE "\t";}

if ($DATA1{$_}{"EXAC_nonTCGA_AC"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonTCGA_AC"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonTCGA_AF"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonTCGA_AF"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonTCGA_GT_Het"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonTCGA_GT_Het"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonTCGA_GF_Het"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonTCGA_GF_Het"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonTCGA_GT_Hom"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonTCGA_GT_Hom"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonTCGA_GF_Hom"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonTCGA_GF_Hom"}."\t";} else {print MYOUTPUTFILE "\t";}

if ($DATA1{$_}{"EXAC_nonPSYCH_AC"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonPSYCH_AC"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonPSYCH_AF"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonPSYCH_AF"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonPSYCH_GT_Het"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonPSYCH_GT_Het"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonPSYCH_GF_Het"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonPSYCH_GF_Het"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonPSYCH_GT_Hom"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonPSYCH_GT_Hom"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"EXAC_nonPSYCH_GF_Hom"}) {print MYOUTPUTFILE $DATA1{$_}{"EXAC_nonPSYCH_GF_Hom"}."\t";} else {print MYOUTPUTFILE "\t";}

if ($DATA1{$_}{"KAV_AC"}) {print MYOUTPUTFILE $DATA1{$_}{"KAV_AC"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"KAV_AF"}) {print MYOUTPUTFILE $DATA1{$_}{"KAV_AF"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"KAV_DataSource"}) {print MYOUTPUTFILE $DATA1{$_}{"KAV_DataSource"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"HRC_AC"}) {print MYOUTPUTFILE $DATA1{$_}{"HRC_AC"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"HRC_AF"}) {print MYOUTPUTFILE $DATA1{$_}{"HRC_AF"}."\t";} else {print MYOUTPUTFILE "\t";}



if ($DATA1{$_}{"Orpha_Inheritance"}) {
my @array1 = @{$DATA1{$_}{"Orpha_Inheritance"}};
foreach (@array1) {
print MYOUTPUTFILE $_." ";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

if ($DATA1{$_}{"Orpha_Onset"}) {
my @array1 = @{$DATA1{$_}{"Orpha_Onset"}};
foreach (@array1) {
print MYOUTPUTFILE $_." ";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

if ($DATA1{$_}{"Orpha_Death"}) {
my @array1 = @{$DATA1{$_}{"Orpha_Death"}};
foreach (@array1) {
print MYOUTPUTFILE $_." ";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

if ($DATA1{$_}{"Orpha_Type"}) {print MYOUTPUTFILE $DATA1{$_}{"Orpha_Type"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Orpha_Name"}) {print MYOUTPUTFILE $DATA1{$_}{"Orpha_Name"}."\t";} else {print MYOUTPUTFILE "\t";}

if ($DATA1{$_}{"Orpha_Prev"}) {
	if (@{$DATA1{$_}{"Orpha_Prev"}} eq "-") {print MYOUTPUTFILE @{$DATA1{$_}{"Orpha_Prev"}}."\t";}
	if (@{$DATA1{$_}{"Orpha_Prev"}} ne "-") {
my @array1 = @{$DATA1{$_}{"Orpha_Prev"}};
my $total_pop = 1;
my $affected  = 0;
foreach (@array1) {
 @tmp1 = split(/\:/,$_);
 if (($tmp1[0] ne "-") and ($tmp1[1] ne "-")) {
 $affected = $affected + ($tmp1[0] * $tmp1[1]);
 $total_pop = $total_pop + $tmp1[0];
 	}
 }
 my $frequency = $affected / $total_pop;
print MYOUTPUTFILE $frequency."\t";
	}
}
else {print MYOUTPUTFILE "\t"};
if ($DATA1{$_}{"OMIM_Inheritance"}) {print MYOUTPUTFILE $DATA1{$_}{"OMIM_Inheritance"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"OMIM_INFO"}) {print MYOUTPUTFILE $DATA1{$_}{"OMIM_INFO"}."\t";} else {print MYOUTPUTFILE "\t";}

if ($DATA1{$_}{"ICD-10_Codes"}) {
my @array1 = @{$DATA1{$_}{"ICD-10_Codes"}};
foreach (@array1) {
print MYOUTPUTFILE $_." ";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

if ($DATA1{$_}{"HPO_Phenotypes"}) {
my @array1 = @{$DATA1{$_}{"HPO_Phenotypes"}};
@array1 = uniq(@array1);
foreach (@array1) {
print MYOUTPUTFILE $_." ";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

if ($DATA1{$_}{"MedGen_Relations"}) {
my @array2 = @{$DATA1{$_}{"MedGen_Relations"}};
@array2 = uniq(@array2);
foreach (@array2) {
print MYOUTPUTFILE $_." ";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

#foreach my $samp (@Samples) {
#	if ($SAMPLE{$samp}{$_}{"Sample_Count"}) {print MYOUTPUTFILE $SAMPLE{$samp}{$_}{"Sample_Count"}."\t";} else {print MYOUTPUTFILE "\t";}
#    if ($SAMPLE{"TCGA_NDAR"}{$_}{"AlleleCount"}) {print MYOUTPUTFILE $SAMPLE{"TCGA_NDAR"}{$_}{"AlleleCount"}."\t";} else {print MYOUTPUTFILE "\t";}
#    if ($SAMPLE{"TCGA_NDAR"}{$_}{"Samples"}) {
#      my @array3 = @{$SAMPLE{"TCGA_NDAR"}{$_}{"Samples"}};
#      foreach (@array3) {
#      print MYOUTPUTFILE $_.",";
#      }
#      print MYOUTPUTFILE "\n";
#    }
#    else {
#      print MYOUTPUTFILE "\n";
#    }

    if ($DATA2{$_}{"Clinical_Samples_Het"}) {
      print MYOUTPUTFILE scalar(@{$DATA2{$_}{"Clinical_Samples_Het"}})."\t";} else {print MYOUTPUTFILE "\t";}
    if ($DATA2{$_}{"Clinical_Samples_Het_Ages"}) {
      my @array4 = @{$DATA2{$_}{"Clinical_Samples_Het_Ages"}};
      print MYOUTPUTFILE "@array4"."\t";
    }
    else {
      print MYOUTPUTFILE "\t";
    }
    if ($DATA2{$_}{"Clinical_Samples_Het"}) {
      my @array3 = @{$DATA2{$_}{"Clinical_Samples_Het"}};
      foreach (@array3) {
      print MYOUTPUTFILE $_.",";
      }
      print MYOUTPUTFILE "\t";
    }
    else {
      print MYOUTPUTFILE "\t";
    }
    if ($DATA2{$_}{"Clinical_Samples_Hom"}) {
      print MYOUTPUTFILE scalar(@{$DATA2{$_}{"Clinical_Samples_Hom"}})."\t";} else {print MYOUTPUTFILE "\t";}
    if ($DATA2{$_}{"Clinical_Samples_Hom_Ages"}) {
      my @array4 = @{$DATA2{$_}{"Clinical_Samples_Hom_Ages"}};
      print MYOUTPUTFILE "@array4"."\t";
    }
    else {
      print MYOUTPUTFILE "\t";
    }
    if ($DATA2{$_}{"Clinical_Samples_Hom"}) {
      my @array3 = @{$DATA2{$_}{"Clinical_Samples_Hom"}};
      foreach (@array3) {
      print MYOUTPUTFILE $_.",";
      }
      print MYOUTPUTFILE "\n";
    }
    else {
      print MYOUTPUTFILE "\n";
    }
    # if ($SAMPLE{$samp}{$_}{"GT_Homs"}) {print MYOUTPUTFILE $SAMPLE{$samp}{$_}{"GT_Homs"}."\t";} else {print MYOUTPUTFILE "\t";}
    # if ($SAMPLE{$samp}{$_}{"Hom_Samples"}) {
    #   my @array4 = @{$SAMPLE{$samp}{$_}{"Hom_Samples"}};
    #   foreach (@array4) {
    #   print MYOUTPUTFILE $_." ";
    #   }
    #   print MYOUTPUTFILE "\t"
    #   }
    #   else {
    #     print MYOUTPUTFILE "\t";
    #   }
#}
#print MYOUTPUTFILE "\n";
};

close(MYOUTPUTFILE);


sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
