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
printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
3	10089729	10089730	C	T	FANCD2	0.000	0.006	0.008	0.268	0.146
6	152419919	152419920	T	C	ESR1	0.000	0.553	0.463	0.383	0.406" > exp
oncogemini truncal \
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

printf "    truncal.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.000	0.493	0.477
5	56183305	56183306	C	T	MAP3K1	0.000	0.416	0.342" > exp
oncogemini truncal \
    --patient C \
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
# 2. Test --maxNorm
###################################################################
printf "testing --maxNorm parameter...\n"
printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
5	56183305	56183306	C	T	MAP3K1	0.200	1.000	1.000	0.543	1.000
6	152419919	152419920	T	C	ESR1	0.000	0.553	0.463	0.383	0.406
16	68842399	68842400	G	C	CDH1	0.062	1.000	0.591	1.000	0.563" > exp
oncogemini truncal \
    --patient B \
    --maxNorm 0.2 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --increase
###################################################################
printf "testing --increase parameter...\n"
printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
6	152419919	152419920	T	C	ESR1	0.000	0.553	0.463	0.383	0.406" > exp
oncogemini truncal \
    --patient B \
    --increase 0.3 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    for(i = 7; i <= 10; i++) printf("%0.3f\t"), $i; printf("%0.3f\n", $11)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    truncal.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.000	0.493	0.477" > exp
oncogemini truncal \
    --patient C \
    --increase 0.4 \
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
# 4. Test expected errors
###################################################################
printf "testing expected errors...\n"
printf "    truncal.no_normal...\n"
echo "There are no normal samples; check the sample manifest file for proper format and loading" > exp
echo "$(oncogemini truncal \
    --patient A \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp

printf "    truncal.no_tumor...\n"
echo "There are no tumor samples; check the sample manifest file for proper format and loading" > exp
echo "$(oncogemini truncal \
    --patient D \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp
