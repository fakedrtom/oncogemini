check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}
export -f check

###################################################################
# 1. Test basic functionality
###################################################################
printf "testing basic functionality...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
oncogemini bottleneck \
    --patient A \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
oncogemini bottleneck \
    --patient B \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
oncogemini bottleneck \
    --patient C \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --minSlope
###################################################################
printf "testing --minSlope parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
oncogemini bottleneck \
    --patient A \
    --minSlope 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value" > exp
oncogemini bottleneck \
    --patient B \
    --minSlope 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
oncogemini bottleneck \
    --patient C \
    --minSlope 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --minR
###################################################################
printf "testing --minR parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
oncogemini bottleneck \
    --patient A \
    --minR 0.8 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
oncogemini bottleneck \
    --patient B \
    --minR 0.8 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
oncogemini bottleneck \
    --patient C \
    --minR 0.8 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 4. Test --minEnd 
###################################################################
printf "testing --minEnd parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
oncogemini bottleneck \
    --patient A \
    --minEnd 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
oncogemini bottleneck \
    --patient B \
    --minEnd 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
oncogemini bottleneck \
    --patient C \
    --minEnd 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 5. Test --endDiff
###################################################################
printf "testing --endDiff parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
oncogemini bottleneck \
    --patient A \
    --endDiff 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value" > exp
oncogemini bottleneck \
    --patient B \
    --endDiff 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
oncogemini bottleneck \
    --patient C \
    --endDiff 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 6. Test --maxNorm
###################################################################
printf "testing --maxNorm parameter...\n"
printf "    bottleneck.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1	alt_AF.A2	alt_AF.A3	slope	intercept	r_value" > exp
oncogemini bottleneck \
    --patient A \
    --maxNorm 0.6 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287	0.0554716911435	-0.025333518655	0.73779842863
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.0	0.0	0.0178571428571	0.35	0.264150943396	0.0878301886792	-0.0492587601078	0.827307416385" > exp
oncogemini bottleneck \
    --patient B \
    --maxNorm 0.6 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
5	112176755	112176756	T	A	APC	0.578313253012	0.567901234568	0.794871794872	0.10827927093	0.538749489887	0.844996899579
6	152419921	152419922	T	A	ESR1	0.0	0.0	0.404255319149	0.202127659574	-0.0673758865248	0.866025403784
16	68842399	68842400	G	C	CDH1	0.0	0.0	0.588235294118	0.294117647059	-0.0980392156863	0.866025403784" > exp
oncogemini bottleneck \
    --patient C \
    --maxNorm 0.6 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp
