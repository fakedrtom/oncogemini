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
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient A \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
cancer_gemini bottleneck \
    --patient B \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
cancer_gemini bottleneck \
    --patient C \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --minDP
###################################################################
printf "testing --minDP parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient A \
    --minDP 100 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863" > exp
cancer_gemini bottleneck \
    --patient B \
    --minDP 100 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient C \
    --minDP 100 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --minGQ
###################################################################
printf "testing --minGQ parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient A \
    --minGQ 99 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
cancer_gemini bottleneck \
    --patient B \
    --minGQ 99 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784" > exp
cancer_gemini bottleneck \
    --patient C \
    --minGQ 99 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 4. Test --minSlope
###################################################################
printf "testing --minSlope parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient A \
    --minSlope 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient B \
    --minSlope 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
cancer_gemini bottleneck \
    --patient C \
    --minSlope 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 5. Test --minR
###################################################################
printf "testing --minR parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient A \
    --minR 0.8 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
cancer_gemini bottleneck \
    --patient B \
    --minR 0.8 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
cancer_gemini bottleneck \
    --patient C \
    --minR 0.8 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 6. Test --minEnd 
###################################################################
printf "testing --minEnd parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient A \
    --minEnd 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
cancer_gemini bottleneck \
    --patient B \
    --minEnd 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
cancer_gemini bottleneck \
    --patient C \
    --minEnd 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 7. Test --endDiff
###################################################################
printf "testing --endDiff parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient A \
    --endDiff 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient B \
    --endDiff 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
cancer_gemini bottleneck \
    --patient C \
    --endDiff 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 8. Test --samples
###################################################################
printf "testing --samples parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A3	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.402985074627	0.47619047619	0.0732054015636	0.402985074627	1.0
5	56168738	56168740	AC	A	MAP3K1	0.433734939759	0.541176470588	0.107441530829	0.433734939759	1.0" > exp
cancer_gemini bottleneck \
    --patient A \
    --samples A1,A3 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

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

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C2	slope	intercept	r_value
5	56168738	56168740	AC	A	MAP3K1	0.0	0.476923076923	0.476923076923	0.0	1.0
5	56183305	56183306	C	T	MAP3K1	0.0	0.341772151899	0.341772151899	0.0	1.0
6	152419921	152419922	T	A	ESR1	0.0	0.404255319149	0.404255319149	0.0	1.0
16	68842399	68842400	G	C	CDH1	0.0	0.588235294118	0.588235294118	0.0	1.0" > exp
cancer_gemini bottleneck \
    --patient C \
    --samples C0,C2 \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

###################################################################
# 9. Test --purity
###################################################################
printf "testing --purity parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	raw.alt_AF.A1	alt_AF.A2	raw.alt_AF.A2	alt_AF.A3	raw.alt_AF.A3	slope	intercept	r_value" > exp
cancer_gemini bottleneck \
    --patient A \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp

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

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	raw.alt_AF.C0	alt_AF.C1	raw.alt_AF.C1	alt_AF.C2	raw.alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.0	0.0	0.505319148936	0.404255319149	0.252659574468	-0.084219858156	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.0	0.0	0.735294117647	0.588235294118	0.367647058824	-0.122549019608	0.866025403784" > exp
cancer_gemini bottleneck \
    --patient C \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    functional_tests_tools.db > obs
check obs exp
rm obs exp
