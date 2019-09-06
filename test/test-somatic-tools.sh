check() 
{
    if diff <( sort "$1" ) <( sort "$2" ); then
        echo ok
    else
        echo fail
        exit 1
    fi
}

###################################################################
# 1. Create somatic DB
###################################################################
printf "creating somatic database, somatic_tests_tools.db...\n"
cancer_gemini set_somatic \
    --minDP 10 \
    --tumAF 0.3 \
    somatic_tests_tools.db

printf "testing is_somatic_B...\n"
echo "6	152419919	152419920	T	C	ESR1
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3" > exp
cancer_gemini query \
    -q "select chrom,start,end,ref,alt,gene from variants where is_somatic_B == 1" \
    somatic_tests_tools.db > obs
check obs exp
rm obs exp

printf "testing is_somatic_C...\n"
echo "5	56168738	56168740	AC	A	MAP3K1
5	56183305	56183306	C	T	MAP3K1
6	152419921	152419922	T	A	ESR1
16	68842399	68842400	G	C	CDH1" > exp
cancer_gemini query \
    -q "select chrom,start,end,ref,alt,gene from variants where is_somatic_C == 1" \
    somatic_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --somatic_only bottleneck
###################################################################
printf "testing --somatic_only parameter...\n"
printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
cancer_gemini bottleneck \
    --patient B \
    --somatic_only \
    --columns "chrom,start,end,ref,alt,gene" \
    somatic_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
cancer_gemini bottleneck \
    --patient C \
    --somatic_only \
    --columns "chrom,start,end,ref,alt,gene" \
    somatic_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --somatic_only truncal
###################################################################
printf "testing --somatic_only parameter...\n"
printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
6	152419919	152419920	T	C	ESR1	0.0	0.552631578947	0.463414634146	0.383333333333	0.40625" > exp
cancer_gemini truncal \
    --patient B \
    --somatic_only \
    --columns "chrom,start,end,ref,alt,gene" \
    somatic_tests_tools.db > obs
check obs exp
rm obs exp

printf "    truncal.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.0	0.492753623188	0.476923076923
5	56183305	56183306	C	T	MAP3K1	0.0	0.415730337079	0.341772151899" > exp
cancer_gemini truncal \
    --patient C \
    --somatic_only \
    --columns "chrom,start,end,ref,alt,gene" \
    somatic_tests_tools.db > obs
check obs exp
rm obs exp

