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

#printf "    truncal.patientA...\n"
#echo "NameError: There are no normal samples; check the ped file for proper format and loading" > exp

#oncogemini truncal \
#    --patient A \
#    --columns "chrom,start,end,ref,alt,gene" \
#    oncogemini_test.db | tail -1 > obs
#check obs exp
#rm obs exp

printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
3	10089729	10089730	C	T	FANCD2	0.0	0.00588235294118	0.00806451612903	0.267605633803	0.146496815287
6	152419919	152419920	T	C	ESR1	0.0	0.552631578947	0.463414634146	0.383333333333	0.40625" > exp
oncogemini truncal \
    --patient B \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    truncal.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.0	0.492753623188	0.476923076923
5	56183305	56183306	C	T	MAP3K1	0.0	0.415730337079	0.341772151899" > exp
oncogemini truncal \
    --patient C \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 2. Test --maxNorm
###################################################################
printf "testing --maxNorm parameter...\n"
printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
5	56183305	56183306	C	T	MAP3K1	0.2	1.0	1.0	0.542857142857	1.0
6	152419919	152419920	T	C	ESR1	0.0	0.552631578947	0.463414634146	0.383333333333	0.40625
16	68842399	68842400	G	C	CDH1	0.0617283950617	1.0	0.590909090909	1.0	0.56338028169" > exp
oncogemini truncal \
    --patient B \
    --maxNorm 0.2 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

###################################################################
# 3. Test --increase
###################################################################
printf "testing --increase parameter...\n"
printf "    truncal.patientB...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.B0	alt_AF.B1	alt_AF.B2	alt_AF.B3	alt_AF.B4
6	152419919	152419920	T	C	ESR1	0.0	0.552631578947	0.463414634146	0.383333333333	0.40625" > exp
oncogemini truncal \
    --patient B \
    --increase 0.3 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
check obs exp
rm obs exp

printf "    truncal.patientC...\n"
echo "chrom	start	end	ref	alt	gene	alt_AF.C0	alt_AF.C1	alt_AF.C2
5	56168738	56168740	AC	A	MAP3K1	0.0	0.492753623188	0.476923076923" > exp
oncogemini truncal \
    --patient C \
    --increase 0.4 \
    --columns "chrom,start,end,ref,alt,gene" \
    oncogemini_test.db > obs
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
