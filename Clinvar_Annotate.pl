#!/usr/bin/perl -w
#use strict;
#use warnings;
#use diagnostics;
use XML::Simple;
use Data::Dumper;
use feature 'say';
my $Variants		= "./ClinVar_SNV_Pathogenic_Germline_Inherited_MOO.txt";
my $MedGen_OMIM		= "./MedGen_HPO_OMIM_Mapping.txt";
my $MG_ICD		    = "./MRCONSO_ENG_ICD10CM.RRF";
my $MG_REL			= "./MRREL.RRF";
my $Output_File		= "./ClinVar_Annotated_Final.txt";
my (@tmp1,@tmp2,@VARs,@HPO_CUIs,@RELAs,@RELAs1,@RELAs2);
my (%DATA1,%DATA2,%DATA3_1,%DATA3_2,%DATA4,%DATA5,%NAMES,$NAME_1,$NAME_2,$NAME_3,%ORPHA_ages,%ORPHA_prev,%ORPHA_data,@CLNDSDB,@CLNDSDBID,@Disorders_data);

my $xml_basic = new XML::Simple;
my $xml_ages = new XML::Simple;
my $xml_prev = new XML::Simple;
my $data_basic = $xml_basic->XMLin("./en_product1.xml");
my $data_ages = $xml_ages->XMLin("./en_product2_ages.xml");
my $data_prev = $xml_prev->XMLin("./en_product2_prev.xml");
#_ICD10CM_###########################################################
 @Disorders_data = keys($data_basic->{DisorderList}->{Disorder});
foreach (@Disorders_data) {
my $OrphaNumber = $data_basic->{DisorderList}->{Disorder}->{$_}->{OrphaNumber};
my $Orpha_Name = $data_basic->{DisorderList}->{Disorder}->{$_}->{Name}->{content};
$ORPHA_data{"$OrphaNumber"}{"Orpha_Name"} = $Orpha_Name;
#print "$Orpha_Name\n";
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
#		print "$ICD_Code";		
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
#		print "$ICD_Code";
		}
	}
elsif ($externalreferences == 0) {
		my $ICD_Code = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}}, $ICD_Code);
		@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"ICD-10 Codes"}});
#        print "$ICD_Code";	
}
#print "\n";

#Orpha_Inheritance###########################################################
my $TypeOfInheritances = 0;
if ($data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{count}) { $TypeOfInheritances = $data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{count};}
if ($TypeOfInheritances > 1) {
my @test = keys(%{ $data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{TypeOfInheritance} });
#print "@test";
foreach my $t (@test) {
		my $inheritance = $data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{TypeOfInheritance}{$t}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}}, $inheritance);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}});
#		print "$inheritance";		
				}
			}
elsif ($TypeOfInheritances == 1) {
		my $inheritance = $data_ages->{DisorderList}->{Disorder}->{$_}->{TypeOfInheritanceList}->{TypeOfInheritance}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}}, $inheritance);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}});
#        print "$inheritance";	
	}
elsif ($TypeOfInheritances == 0) {
		my $inheritance = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}}, $inheritance);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Inheritance"}});
#        print "$inheritance";	
	}	
#print "\n";

#Orpha_Age_Onset###########################################################
my $AgeOfOnsets = 0;
if ($data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{count}) { $AgeOfOnsets = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{count};}
if ($AgeOfOnsets > 1) {
my @test = keys(%{ $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{AverageAgeOfOnset} });
foreach my $t (@test) {
		my $onset = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{AverageAgeOfOnset}{$t}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}}, $onset);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}});
#		print "$onset";		
				}
			}
elsif ($AgeOfOnsets == 1) {
		my $onset = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfOnsetList}->{AverageAgeOfOnset}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}}, $onset);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}});
#		print "$onset";
	}
elsif ($AgeOfOnsets == 0) {
		my $onset = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}}, $onset);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Onset"}});
#		print "$onset";
	}
#print "\n";		
#Orpha_Age_Death###########################################################
my $AgeOfDeaths = 0;
if ($data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{count}) { $AgeOfDeaths = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{count};}
if ($AgeOfDeaths > 1) {
my @test = keys(%{ $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{AverageAgeOfDeath} });
foreach my $t (@test) {
		my $death = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{AverageAgeOfDeath}{$t}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}}, $death);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}});
#		print "$death";
				}
			}
elsif ($AgeOfDeaths == 1) {
		my $death = $data_ages->{DisorderList}->{Disorder}->{$_}->{AverageAgeOfDeathList}->{AverageAgeOfDeath}->{Name}->{content};
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}}, $death);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}});
#		print "$death";
	}
elsif ($AgeOfDeaths == 0) {
		my $death = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}}, $death);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Death"}});
#		print "$death";
	}
#print "\n";		
#Orpha_Disorder_Type###########################################################
my $TypeOfDisorder = $data_ages->{DisorderList}->{Disorder}->{$_}->{DisorderType}->{Name}->{content};
$ORPHA_data{"$OrphaNumber"}{"Orpha_Type"} = $TypeOfDisorder;
#print "$TypeOfDisorder";
#print "\n";
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
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}}, "$type: $prevalence: $location");
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}});
#		print "$type: $prevalence: $location";
				}
			}
		}	
elsif ($Prevelances == 1) {
		my $prevalence = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceClass}->{Name}->{content};
		my $type = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceType}->{Name}->{content};
		my $location = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceGeographic}->{Name}->{content};
		if ($data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceClass}->{Name}->{content}) {
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}}, "$type: $prevalence: $location");
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}});
#		print "$type: $prevalence: $location";
		}
	}
elsif ($Prevelances == 0) {
		my $type = "-";
		push (@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}}, $type);
		@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}} = uniq(@{$ORPHA_data{"$OrphaNumber"}{"Orpha_Prev"}});
#		print "$type";
	}
#print "\n";		
};

#foreach (@Disorders_prev) {
#my $OrphaNumber = $data_prev->{DisorderList}->{Disorder}->{$_}->{OrphaNumber};
#if ($data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}) {
#my @Orpha_prev_IDs = keys($data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence});
#		if ($Orpha_prev_IDs[0] eq ("PrevalenceGeographic" or "PrevalenceValidationStatus" or "PrevalenceQualification" or "ValMoy" or "id" or "Source" or "PrevalenceClass" or "PrevalenceType")) {
#	my $prevalence = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceClass}->{Name}->{content};
#	my $location = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{PrevalenceGeographic}->{Name}->{content};
#		if ($location and $prevalence) {
#	push (@{$ORPHA_prev{"$OrphaNumber"}{"Orpha_Prevalence"}}, "$location $prevalence");
#	@{$ORPHA_prev{"$OrphaNumber"}{"Orpha_Prevalence"}} = uniq(@{$ORPHA_prev{"$OrphaNumber"}{"Orpha_Prevalence"}});
#		}
#	}
#		else {
#	foreach $NAME_2 (@Orpha_prev_IDs) {
#	my $prevalence = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{$NAME_2}->{PrevalenceClass}->{Name}->{content};
#	my $location = $data_prev->{DisorderList}->{Disorder}->{$_}->{PrevalenceList}->{Prevalence}->{$NAME_2}->{PrevalenceGeographic}->{Name}->{content};
#		if ($location and $prevalence) {
#	push (@{$ORPHA_prev{"$OrphaNumber"}{"Orpha_Prevalence"}}, "$location $prevalence");
#	@{$ORPHA_prev{"$OrphaNumber"}{"Orpha_Prevalence"}} = uniq(@{$ORPHA_prev{"$OrphaNumber"}{"Orpha_Prevalence"}});
#				}
#			}
#		}
#	}
#};

##Parse ClinVar VCF File
open (LIST, $Variants) || die "File not found\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);         
		 $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Chr"} = $tmp1[0];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Pos"} = $tmp1[1];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"dbSNP"} = $tmp1[2];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Ref"} = $tmp1[3];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Alt"} = $tmp1[4];
		 my @DATA = split(/\;/, $tmp1[7]);
         if (scalar(@DATA) > 0) {
         foreach my $data (@DATA) {
         my @Current = split(/\=/, $data);
         if (scalar(@Current) > 0) {
         if ($Current[0] eq "CLNSRCID") {
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"CLNSRCID"} = $Current[1];
         		}
         elsif ($Current[0] eq "CLNDBN") {
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"CLNDBN"} = $Current[1];
         		}
         elsif ($Current[0] eq "CLNACC") {
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"CLNACC"} = $Current[1];
         		}
         elsif ($Current[0] eq "CLNDSDB") {
         @CLNDSDB = split(/\:/,$Current[1]);
         		}
         elsif ($Current[0] eq "CLNDSDBID") {
         @CLNDSDBID = split(/\:/,$Current[1]);
         		}						
			}
		}	
	}
	my $position = 0;
	foreach (@CLNDSDB) {
	if ($CLNDSDB[$position] eq "MedGen") { 
	$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"MedGen"} = $CLNDSDBID[$position];
	}
	elsif ($CLNDSDB[$position] eq "OMIM") { 
	$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"OMIM"} = $CLNDSDBID[$position];
	}
	elsif ($CLNDSDB[$position] eq "Orphanet") {
	my @ORPHA_ID = split('ORPHA',$CLNDSDBID[$position]); 
	$DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]\t$tmp1[3]\t$tmp1[4]"}{"Orphanet"} = $ORPHA_ID[1];
	}
	$position++;
		}		
};
close(LIST);
####################################

#Parse MedGen_OMIM Relations
open (LIST, $MedGen_OMIM) || die "File not found\n";
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\|/, $_);
         if ($tmp1[1] ne "MIM_number"){
		 $DATA3_1{"$tmp1[0]"}{"$tmp1[4]"}{"OMIM_CUI"} = $tmp1[0];
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
open (LIST, $MG_REL) || die "File not found\n";
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

open (LIST, $MG_ICD) || die "File not found\n";
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


@VARs =keys(%DATA1);
open(MYOUTPUTFILE, ">$Output_File");
print MYOUTPUTFILE "Chr\tPos\tdbSNP\tRef\tAlt\tCLNSRCID\tCLNDBN\tCLNACC\tMedGen\tOMIM\tOrphanet\tOrpha_Inheritance\tOrpha_Age_Onset\tOrpha_Age_Death\tOrpha_Disorder_Type\tOrpha_Name\tOrpha_Prevalence\tICD-10\tHPO_Terms\tICD10CM_SubTerms\n";
foreach (@VARs) {
if ($DATA1{$_}{"Chr"}) {print MYOUTPUTFILE $DATA1{$_}{"Chr"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Pos"}) {print MYOUTPUTFILE $DATA1{$_}{"Pos"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"dbSNP"}) {print MYOUTPUTFILE $DATA1{$_}{"dbSNP"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Ref"}) {print MYOUTPUTFILE $DATA1{$_}{"Ref"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Alt"}) {print MYOUTPUTFILE $DATA1{$_}{"Alt"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNSRCID"}) {print MYOUTPUTFILE $DATA1{$_}{"CLNSRCID"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNDBN"}) {print MYOUTPUTFILE $DATA1{$_}{"CLNDBN"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"CLNACC"}) {print MYOUTPUTFILE $DATA1{$_}{"CLNACC"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"MedGen"}) {print MYOUTPUTFILE $DATA1{$_}{"MedGen"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"OMIM"}) {print MYOUTPUTFILE $DATA1{$_}{"OMIM"}."\t";} else {print MYOUTPUTFILE "\t";}
if ($DATA1{$_}{"Orphanet"}) {print MYOUTPUTFILE $DATA1{$_}{"Orphanet"}."\t";} else {print MYOUTPUTFILE "\t";}

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
my @array1 = @{$DATA1{$_}{"Orpha_Prev"}};
foreach (@array1) {
print MYOUTPUTFILE $_." ";
}
print MYOUTPUTFILE "\t";
}
else {print MYOUTPUTFILE "\t"};

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
print MYOUTPUTFILE "\n";
}
else {print MYOUTPUTFILE "\n"};
  	};
  	
close(MYOUTPUTFILE);


sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}
