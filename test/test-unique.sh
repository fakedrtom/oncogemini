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
printf "    unique.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A3" > exp
oncogemini unique \
    --patient A \
    --specific A3 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    unique.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B1" > exp
oncogemini unique \
    --patient B \
    --specific B1 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
6	152419921	152419922	T	A	ESR1	0.404
16	68842399	68842400	G	C	CDH1	0.588" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
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
# 2. Test --maxOthers
###################################################################
printf "testing --maxOthers parameter...\n"
printf "    unique.patientA...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.A1
16	68842399	68842400	G	C	CDH1	0.606" > exp
oncogemini unique \
    --patient A \
    --specific A1 \
    --maxOthers 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\n", $7)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    unique.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B1
6	152419919	152419920	T	C	ESR1	0.553" > exp
oncogemini unique \
    --patient B \
    --specific B1 \
    --maxOthers 0.5 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db | \
    awk '{if ($1=="chrom") print $0}; \
    {if ($1 != "chrom") \
    printf("%d\t" "%d\t" "%d\t" "%s\t" "%s\t" "%s\t", $1,$2,$3,$4,$5,$6); \
    printf("%0.3f\n", $7)}' | \
    grep -v "^0.000" > obs
check obs exp
rm obs exp

printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
16	68842399	68842400	G	C	CDH1	0.588" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
    --maxOthers 0.5 \
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
# 3. Test --increase
###################################################################
printf "testing --increase parameter...\n"
printf "    unique.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C2
16	68842399	68842400	G	C	CDH1	0.588" > exp
oncogemini unique \
    --patient C \
    --specific C2 \
    --increase 0.5 \
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
# 4. Test expected errors
###################################################################
printf "testing expected errors...\n"
printf "    unique.specificA5...\n"
echo "Sample listed with --specific, A5, is not found, check the sample manifest file for available samples" > exp
echo "$(oncogemini unique \
    --patient A \
    --specific A5 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp

printf "    unique.no_specific...\n"
echo "No sample(s) specified with --specific, please provide sample(s)" > exp
echo "$(oncogemini unique \
    --patient A \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp

printf "    unique.no_others...\n"
echo "There are no other samples to compare --specific samples to; check the sample manifest file for proper format and loading" > exp
echo "$(oncogemini unique \
    --patient D \
    --specific D0,D1 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db 2>&1 > /dev/null)" > obs
check obs exp
rm obs exp
