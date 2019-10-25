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
# 1. Create somatic DB
###################################################################
printf "creating somatic database, somatic_tests_tools.db...\n"
oncogemini set_somatic \
    --minDP 10 \
    --tumAF 0.3 \
    somatic_tests_tools.db

printf "testing is_somatic_B...\n"
echo "6	152419919	152419920	T	C	ESR1
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3" > exp
oncogemini query \
    -q "select chrom,start,end,ref,alt,gene from variants where is_somatic_B == 1" \
    somatic_tests_tools.db > obs
check obs exp
rm obs exp

printf "testing is_somatic_C...\n"
echo "5	56168738	56168740	AC	A	MAP3K1
5	56183305	56183306	C	T	MAP3K1
6	152419921	152419922	T	A	ESR1
16	68842399	68842400	G	C	CDH1" > exp
oncogemini query \
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
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.000	0.018	0.350	0.264	0.088	-0.049	0.827" > exp
oncogemini bottleneck \
    --patient B \
    --somatic_only \
    --columns "chrom,start,end,ref,alt,gene" \
    somatic_tests_tools.db | \
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
    --somatic_only \
    --columns "chrom,start,end,ref,alt,gene" \
    somatic_tests_tools.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --somatic_only truncal
###################################################################
printf "testing --somatic_only parameter...\n"
printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
6	152419919	152419920	T	C	ESR1	0.000	0.553	0.463	0.383	0.406" > exp
oncogemini truncal \
    --patient B \
    --somatic_only \
    --columns "chrom,start,end,ref,alt,gene" \
    somatic_tests_tools.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    truncal.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.000	0.493	0.477
5	56183305	56183306	C	T	MAP3K1	0.000	0.416	0.342" > exp
oncogemini truncal \
    --patient C \
    --somatic_only \
    --columns "chrom,start,end,ref,alt,gene" \
    somatic_tests_tools.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 8; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $9)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

