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
# 1. Test --minDP
###################################################################
printf "testing --minDP parameter...\n"
printf "    bottleneck.patientB.minDP10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
cancer_gemini bottleneck \
    --patient B \
    --minDP 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB.minDP100...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863" > exp
cancer_gemini bottleneck \
    --patient B \
    --minDP 100 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    loh.patientB.minDP10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
13	32912963	32912968	TGAAA	T	BRCA2	0.507246376812	0.986666666667	0.975	0.890625	0.875" > exp
cancer_gemini loh \
    --patient B \
    --minDP 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    loh.patientB.minDP49...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
cancer_gemini loh \
    --patient B \
    --minDP 49 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    truncal.patientC.minDP10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.0	0.492753623188	0.476923076923
5	56183305	56183306	C	T	MAP3K1	0.0	0.415730337079	0.341772151899" > exp
cancer_gemini truncal \
    --patient C \
    --minDP 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    truncal.patientC.minDP79...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56183305	56183306	C	T	MAP3K1	0.0	0.415730337079	0.341772151899" > exp
cancer_gemini truncal \
    --patient C \
    --minDP 79 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientC.minDP10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404255319149
16	68842399	68842400	G	C	CDH1	0.588235294118" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --minDP 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientC.minDP35...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404255319149" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --minDP 35 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --minGQ
###################################################################
printf "testing --minGQ parameter...\n"

printf "    bottleneck.patientB.minGQ10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
cancer_gemini bottleneck \
    --patient B \
    --minGQ 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB.minGQ141...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863" > exp
cancer_gemini bottleneck \
    --patient B \
    --minGQ 141 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    loh.patientB.minGQ10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
13	32912963	32912968	TGAAA	T	BRCA2	0.507246376812	0.986666666667	0.975	0.890625	0.875" > exp
cancer_gemini loh \
    --patient B \
    --minGQ 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    loh.patientB.minGQ20...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
cancer_gemini loh \
    --patient B \
    --minGQ 20 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    truncal.patientB.minGQ10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287
6	152419919	152419920	T	C	ESR1	0.0	0.552631578947	0.463414634146	0.383333333333	0.40625" > exp
cancer_gemini truncal \
    --patient B \
    --minGQ 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    truncal.patientB.minGQ160...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
6	152419919	152419920	T	C	ESR1	0.0	0.552631578947	0.463414634146	0.383333333333	0.40625" > exp
cancer_gemini truncal \
    --patient B \
    --minGQ 160 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientC.minGQ10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404255319149
16	68842399	68842400	G	C	CDH1	0.588235294118" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --minGQ 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientC.minGQ144...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404255319149" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --minGQ 144 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp
 
###################################################################
# 3. Test --samples
###################################################################
printf "testing --samples parameter...\n"

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B1	alt_AF.B2	alt_AF.B3	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.00588235294118	0.00806451612903	0.267605633803	0.130861640431	-0.0370108061398	0.869627975964
6	152419921	152419922	T	A	ESR1	0.0	1.0	0.571428571429	0.285714285714	0.238095238095	0.569494797451
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0178571428571	0.35	0.175	-0.052380952381	0.88778411242" > exp
cancer_gemini bottleneck \
    --patient B \
    --samples B1,B2,B3 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B2	alt_AF.B3
13	32912963	32912968	TGAAA	T	BRCA2	0.507246376812	0.975	0.890625" > exp
cancer_gemini loh \
    --patient B \
    --samples B0,B2,B3 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B4
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.146496815287
6	152419919	152419920	T	C	ESR1	0.0	0.552631578947	0.40625" > exp
cancer_gemini truncal \
    --patient B \
    --samples B0,B1,B4 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C1
5	56168738	56168740	AC	A	MAP3K1	0.492753623188
5	56183305	56183306	C	T	MAP3K1	0.415730337079" > exp
cancer_gemini unique \
    --patient C \
    --samples C0,C1 \
    --specific C1 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 4. Test --purity
###################################################################
printf "testing --purity parameter...\n"

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	raw.alt_AF.B0	alt_AF.B1	raw.alt_AF.B1	alt_AF.B2	raw.alt_AF.B2	alt_AF.B3	raw.alt_AF.B3	alt_AF.B4	raw.alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.0	0.0588235294118	0.00588235294118	0.0806451612903	0.00806451612903	0.892018779343	0.267605633803	0.292993630573	0.146496815287	0.141918251108	-0.0189402820919	0.610347968039
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0	0.0	0.178571428571	0.0178571428571	1	0.35	0.528301886792	0.264150943396	0.205660377358	-0.0699460916442	0.762067229018" > exp
cancer_gemini bottleneck \
    --patient B \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	raw.alt_AF.B0	alt_AF.B1	raw.alt_AF.B1	alt_AF.B2	raw.alt_AF.B2	alt_AF.B3	raw.alt_AF.B3	alt_AF.B4	raw.alt_AF.B4
13	32912963	32912968	TGAAA	T	BRCA2	0.507246376812	0.507246376812	1	0.986666666667	1	0.975	1	0.890625	1	0.875" > exp
cancer_gemini loh \
    --patient B \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    truncal.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	raw.alt_AF.C0	alt_AF.C1	raw.alt_AF.C1	alt_AF.C2	raw.alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.0	0.0	1	0.492753623188	0.596153846154	0.476923076923
5	56183305	56183306	C	T	MAP3K1	0.0	0.0	1	0.415730337079	0.427215189873	0.341772151899" > exp
cancer_gemini truncal \
    --patient C \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2	raw.alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.505319148936	0.404255319149
16	68842399	68842400	G	C	CDH1	0.735294117647	0.588235294118" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 5. Test --cancers
###################################################################
printf "testing --cancers parameter...\n"

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863" > exp
cancer_gemini bottleneck \
    --patient B \
    --cancers AML,ST \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
cancer_gemini loh \
    --patient B \
    --cancers AML,ST \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287" > exp
cancer_gemini truncal \
    --patient B \
    --cancers AML,ST \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "     unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
16	68842399	68842400	G	C	CDH1	0.588235294118" > exp
cancer_gemini unique \
    --patient C \
    --specific C2 \
    --cancers AML,ST \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp
