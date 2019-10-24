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
# 1. Test --minDP
###################################################################
printf "testing --minDP parameter...\n"
printf "    bottleneck.patientB.minDP10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146	0.055	-0.025	0.738
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.000	0.018	0.350	0.264	0.088	-0.049	0.827" > exp
oncogemini bottleneck \
    --patient B \
    --minDP 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB.minDP100...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146	0.055	-0.025	0.738" > exp
oncogemini bottleneck \
    --patient B \
    --minDP 100 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientB.minDP10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
13	32912963	32912968	TGAAA	T	BRCA2	0.507	0.987	0.975	0.891	0.875" > exp
oncogemini loh \
    --patient B \
    --minDP 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientB.minDP49...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
oncogemini loh \
    --patient B \
    --minDP 49 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    truncal.patientC.minDP10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.000	0.493	0.477
5	56183305	56183306	C	T	MAP3K1	0.000	0.416	0.342" > exp
oncogemini truncal \
    --patient C \
    --minDP 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 8; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $9)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    truncal.patientC.minDP79...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56183305	56183306	C	T	MAP3K1	0.000	0.416	0.342" > exp
oncogemini truncal \
    --patient C \
    --minDP 79 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 8; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $9)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    unique.patientC.minDP10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404
16	68842399	68842400	G	C	CDH1	0.588" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
    --minDP 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\n", $7)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    unique.patientC.minDP35...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
    --minDP 35 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\n", $7)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --minGQ
###################################################################
printf "testing --minGQ parameter...\n"

printf "    bottleneck.patientB.minGQ10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146	0.055	-0.025	0.738
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.000	0.018	0.350	0.264	0.088	-0.049	0.827" > exp
oncogemini bottleneck \
    --patient B \
    --minGQ 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    bottleneck.patientB.minGQ141...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146	0.055	-0.025	0.738" > exp
oncogemini bottleneck \
    --patient B \
    --minGQ 141 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientB.minGQ10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
13	32912963	32912968	TGAAA	T	BRCA2	0.507	0.987	0.975	0.891	0.875" > exp
oncogemini loh \
    --patient B \
    --minGQ 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientB.minGQ20...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
oncogemini loh \
    --patient B \
    --minGQ 20 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    truncal.patientB.minGQ10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146
6	152419919	152419920	T	C	ESR1	0.000	0.553	0.463	0.383	0.406" > exp
oncogemini truncal \
    --patient B \
    --minGQ 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    truncal.patientB.minGQ160...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
6	152419919	152419920	T	C	ESR1	0.000	0.553	0.463	0.383	0.406" > exp
oncogemini truncal \
    --patient B \
    --minGQ 160 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    unique.patientC.minGQ10...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404
16	68842399	68842400	G	C	CDH1	0.588" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
    --minGQ 10 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\n", $7)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    unique.patientC.minGQ144...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
    --minGQ 144 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\n", $7)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp
 
###################################################################
# 3. Test --samples
###################################################################
printf "testing --samples parameter...\n"

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B1	alt_AF.B2	alt_AF.B3	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.006	0.008	0.268	0.131	-0.037	0.870
6	152419921	152419922	T	A	ESR1	0.000	1.000	0.571	0.286	0.238	0.569
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.018	0.350	0.175	-0.052	0.888" > exp
oncogemini bottleneck \
    --patient B \
    --samples B1,B2,B3 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B2	alt_AF.B3
13	32912963	32912968	TGAAA	T	BRCA2	0.507	0.975	0.891" > exp
oncogemini loh \
    --patient B \
    --samples B0,B2,B3 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 8; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $9)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B4
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.146
6	152419919	152419920	T	C	ESR1	0.000	0.553	0.406" > exp
oncogemini truncal \
    --patient B \
    --samples B0,B1,B4 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 8; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $9)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C1
5	56168738	56168740	AC	A	MAP3K1	0.493
5	56183305	56183306	C	T	MAP3K1	0.416" > exp
oncogemini unique \
    --patient C \
    --samples C0,C1 \
    --specific C1 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\n", $7)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

###################################################################
# 4. Test --purity
###################################################################
printf "testing --purity parameter...\n"

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	raw.alt_AF.B0	alt_AF.B1	raw.alt_AF.B1	alt_AF.B2	raw.alt_AF.B2	alt_AF.B3	raw.alt_AF.B3	alt_AF.B4	raw.alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.000	0.000	0.009	0.006	0.011	0.008	0.335	0.268	0.293	0.146	0.091	-0.053	0.853
9	102595635	102595654	AACCTTCTCAGCCCTCTCC	A	NR4A3	0.000	0.000	0.000	0.000	0.023	0.018	0.438	0.350	0.528	0.264	0.149	-0.101	0.900" > exp
oncogemini bottleneck \
    --patient B \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 18; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $19)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	raw.alt_AF.B0	alt_AF.B1	raw.alt_AF.B1	alt_AF.B2	raw.alt_AF.B2	alt_AF.B3	raw.alt_AF.B3	alt_AF.B4	raw.alt_AF.B4
13	32912963	32912968	TGAAA	T	BRCA2	0.507	0.507	1.000	0.987	1.000	0.975	1.000	0.891	1.000	0.875" > exp
oncogemini loh \
    --patient B \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 15; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $16)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    truncal.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	raw.alt_AF.C0	alt_AF.C1	raw.alt_AF.C1	alt_AF.C2	raw.alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.000	0.000	1.000	0.493	0.596	0.477
5	56183305	56183306	C	T	MAP3K1	0.000	0.000	1.000	0.416	0.427	0.342" > exp
oncogemini truncal \
    --patient C \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 11; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $12)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2	raw.alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.505	0.404
16	68842399	68842400	G	C	CDH1	0.735	0.588" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
    --purity \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\t" "%0.3f\n", $7,$8)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

###################################################################
# 5. Test --cancers
###################################################################
printf "testing --cancers parameter...\n"

printf "    bottleneck.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4	slope	intercept	r_value
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146	0.055	-0.025	0.738" > exp
oncogemini bottleneck \
    --patient B \
    --cancers AML,ST \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 13; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $14)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
oncogemini loh \
    --patient B \
    --cancers AML,ST \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146" > exp
oncogemini truncal \
    --patient B \
    --cancers AML,ST \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "     unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
16	68842399	68842400	G	C	CDH1	0.588" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
    --cancers AML,ST \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\n", $7)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

###################################################################
# 6. Test expected errors
###################################################################
printf "testing expected errors...\n"
printf "    bottleneck.no_patient...\n"
echo "More than 1 patient is present, specify a patient_id with --patient" > exp
echo "$(oncogemini bottleneck \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp

printf "    bottleneck.wrong_patient...\n"
echo "Specified patient is not found, check the sample manifest file for available patient_ids" > exp
echo "$(oncogemini bottleneck \
    --patient Tom \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp

printf "    bottleneck.wrong_sample...\n"
echo "Specified samples, Tom, is not found" > exp
echo "$(oncogemini bottleneck \
    --patient B \
    --sample Tom \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp
