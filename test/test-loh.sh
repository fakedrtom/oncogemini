check() 
{
    if diff <( sort "$1" ) <( sort "$2" ); then
        echo ok
    else
        echo fail
        exit 1
    fi
}
#export -f check

###################################################################
# 1. Test basic functionality
###################################################################
printf "testing basic functionality...\n"

#printf "    loh.patientA...\n"
#echo "NameError: There are no normal samples; check the ped file for proper format and loading" > exp

#cancer_gemini loh \
#    --patient A \
#    --columns "chrom,start,end,ref,alt,gene" \
#    functional_tests_tools.db | tail -1 > obs
#check obs exp
#rm obs exp

printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
13	32912963	32912968	TGAAA	T	BRCA2	0.507246376812	0.986666666667	0.975	0.890625	0.875" > exp
cancer_gemini loh \
    --patient B \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    loh.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2" > exp
cancer_gemini loh \
    --patient C \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --maxNorm
###################################################################
printf "testing --maxNorm parameter...\n"
printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
cancer_gemini loh \
    --patient B \
    --maxNorm 0.507 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --minNorm
###################################################################
printf "testing --minNorm parameter...\n"
printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
cancer_gemini loh \
    --patient B \
    --minNorm 0.508 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 4. Test --minTumor
###################################################################
printf "testing --minTumor parameter...\n"
printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
cancer_gemini loh \
    --patient B \
    --minTumor 0.89 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "testing --minTumor parameter...\n"
printf "    loh.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	112176755	112176756	T	A	APC	0.578313253012	0.567901234568	0.794871794872" > exp
cancer_gemini loh \
    --patient C \
    --minTumor 0.55 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp
