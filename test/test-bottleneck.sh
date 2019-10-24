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
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146	0.055	-0.025	0.738
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.000	0.018	0.350	0.264	0.088	-0.049	0.827" > exp
oncogemini bottleneck \
    --patient B \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.000	0.000	0.404	0.202	-0.067	0.866
16	68842399	68842400	G	C	CDH1	0.000	0.000	0.588	0.294	-0.098	0.866" > exp
oncogemini bottleneck \
    --patient C \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
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
16	68842399	68842400	G	C	CDH1	0.000	0.000	0.588	0.294	-0.098	0.866" > exp
oncogemini bottleneck \
    --patient C \
    --minSlope 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
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
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.000	0.018	0.350	0.264	0.088	-0.049	0.827" > exp
oncogemini bottleneck \
    --patient B \
    --minR 0.8 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.000	0.000	0.404	0.202	-0.067	0.866
16	68842399	68842400	G	C	CDH1	0.000	0.000	0.588	0.294	-0.098	0.866" > exp
oncogemini bottleneck \
    --patient C \
    --minR 0.8 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
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
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.000	0.018	0.350	0.264	0.088	-0.049	0.827" > exp
oncogemini bottleneck \
    --patient B \
    --minEnd 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
6	152419921	152419922	T	A	ESR1	0.000	0.000	0.404	0.202	-0.067	0.866
16	68842399	68842400	G	C	CDH1	0.000	0.000	0.588	0.294	-0.098	0.866" > exp
oncogemini bottleneck \
    --patient C \
    --minEnd 0.25 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
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
16	68842399	68842400	G	C	CDH1	0.000	0.000	0.588	0.294	-0.098	0.866" > exp
oncogemini bottleneck \
    --patient C \
    --endDiff 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
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
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146	0.055	-0.025	0.738
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.000	0.018	0.350	0.264	0.088	-0.049	0.827" > exp
oncogemini bottleneck \
    --patient B \
    --maxNorm 0.6 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    bottleneck.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2	slope	intercept	r_value
5	112176755	112176756	T	A	APC	0.578	0.568	0.795	0.108	0.539	0.845
6	152419921	152419922	T	A	ESR1	0.000	0.000	0.404	0.202	-0.067	0.866
16	68842399	68842400	G	C	CDH1	0.000	0.000	0.588	0.294	-0.098	0.866" > exp
oncogemini bottleneck \
    --patient C \
    --maxNorm 0.6 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp
