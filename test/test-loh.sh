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
printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
13	32912963	32912968	TGAAA	T	BRCA2	0.507	0.987	0.975	0.891	0.875" > exp
oncogemini loh \
    --patient B \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2" > exp
oncogemini loh \
    --patient C \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --maxNorm
###################################################################
printf "testing --maxNorm parameter...\n"
printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
oncogemini loh \
    --patient B \
    --maxNorm 0.507 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --minNorm
###################################################################
printf "testing --minNorm parameter...\n"
printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
oncogemini loh \
    --patient B \
    --minNorm 0.508 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 4. Test --minTumor
###################################################################
printf "testing --minTumor parameter...\n"
printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4" > exp
oncogemini loh \
    --patient B \
    --minTumor 0.89 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "testing --minTumor parameter...\n"
printf "    loh.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	112176755	112176756	T	A	APC	0.578	0.568	0.795" > exp
oncogemini loh \
    --patient C \
    --minTumor 0.55 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 8; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $9)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

###################################################################
# 5. Test --specific
###################################################################
printf "testing --specific parameter...\n"
printf "    loh.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A2	alt_AF.A3
5	112176755	112176756	T	A	APC	0.500	0.983
6	152419921	152419922	T	A	ESR1	0.407	1.000" > exp
oncogemini loh \
    --patient A \
    --specific A3 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\t" "%0.3f\n", $7,$8)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B3	alt_AF.B4
5	56183305	56183306	C	T	MAP3K1	0.543	1.000" > exp
oncogemini loh \
    --patient B \
    --specific B4 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\t"  "%0.3f\n", $7,$8)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    loh.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1
6	152419919	152419920	T	C	ESR1	0.391	1.000" > exp
oncogemini loh \
    --patient C \
    --specific C1 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\t"  "%0.3f\n", $7,$8)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

###################################################################
# 5. Test expected errors
###################################################################
printf "testing expected errors...\n"
printf "    loh.specifcB5...\n"
echo "Error: Specified sample name with --specific is not found, make sure a single sample only is provided and check the sample manifest file for available sample names" > exp
echo "$(oncogemini loh \
    --patient B \
    --specific B5 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp

printf "    loh.specificB0...\n"
echo "Error: Specified sample with --specific is the first timepoint, specify a sample that has a preceding sample" > exp
echo "$(oncogemini loh \
    --patient B \
    --specific B0 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp

printf "    loh.no_normal...\n"
echo "Error: There are no normal samples; check the sample manifest file for proper format and loading" > exp
echo "$(oncogemini loh \
    --patient A \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp

printf "    loh.no_tumor...\n"
echo 'Error: There are no tumor samples; check the sample manifest file for proper format and loading' > exp
echo "$(oncogemini loh \
    --patient D \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp
