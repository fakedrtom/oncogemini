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
printf "    unique.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A3" > exp
cancer_gemini unique \
    --patient A \
    --specific A3 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B1" > exp
cancer_gemini unique \
    --patient B \
    --specific B1 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404255319149
16	68842399	68842400	G	C	CDH1	0.588235294118" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --maxOthers
###################################################################
printf "testing --maxOthers parameter...\n"
printf "    unique.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1
16	68842399	68842400	G	C	CDH1	0.606060606061" > exp
cancer_gemini unique \
    --patient A \
    --specific A1 \
    --maxOthers 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B1
6	152419919	152419920	T	C	ESR1	0.552631578947" > exp
cancer_gemini unique \
    --patient B \
    --specific B1 \
    --maxOthers 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
16	68842399	68842400	G	C	CDH1	0.588235294118" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --maxOthers 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --increase
###################################################################
printf "testing --increase parameter...\n"
printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
16	68842399	68842400	G	C	CDH1	0.588235294118" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --increase 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp
