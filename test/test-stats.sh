check() 
{
    if diff <( sort "$1" ) <( sort "$2" ); then
        echo ok
    else
        echo fail
        exit 1
    fi
}

###################################################################
# 1. Test variants by sample
###################################################################
echo "    stat.t01...\c"
echo "sample	total
1094PC0005	1
1094PC0009	1
1094PC0012	5
1094PC0013	4" > exp
oncogemini stats --vars-by-sample test3.snpeff.db \
       > obs
check obs exp
rm obs exp

###################################################################
# 2. Test genotypes by sample
###################################################################
echo "    stat.t02...\c"
echo "sample	num_hom_ref	num_het	num_hom_alt	num_unknown	total
1094PC0005	4	0	1	3	8
1094PC0009	4	0	1	3	8
1094PC0012	3	2	3	0	8
1094PC0013	3	2	2	1	8" > exp
oncogemini stats --gts-by-sample test3.snpeff.db \
       > obs
check obs exp
rm obs exp

###################################################################
# 3. Test site freq. spectrum
###################################################################
echo "    stat.t03...\c"
echo "aaf	count
0.0	3
0.25	1
0.5	1
1.0	3" > exp
oncogemini stats --sfs test3.snpeff.db \
       > obs
check obs exp
rm obs exp

###################################################################
# 4. Test snp counts
###################################################################
echo "    stat.t04...\c"
echo "type	count
A->G	2
C->T	1
G->A	1
G->C	1
G->T	1" > exp
oncogemini stats --snp-counts test3.snpeff.db \
        > obs

check obs exp
rm obs exp

###################################################################
# 5. Test transition / transversion ratios
###################################################################
echo "    stat.t05...\c"
echo "ts	tv	ts/tv
4	5	0.8" > exp
oncogemini stats --tstv test.snpeff.vcf.db > obs
check obs exp
rm obs exp

###################################################################
# 6. Test transition / transversion ratios in coding region
###################################################################
echo "    stat.t06...\c"
echo "ts	tv	ts/tv
3	2	1.5" > exp
oncogemini stats --tstv-coding test.snpeff.vcf.db > obs
check obs exp
rm obs exp

###################################################################
# 7. Test transition / transversion ratios in noncoding region
###################################################################
echo "    stat.t07...\c"
echo "ts	tv	ts/tv
1	3	0.3333" > exp
oncogemini stats --tstv-noncoding test.snpeff.vcf.db > obs
check obs exp
rm obs exp


###################################################################
# 8. Test multi-dimensional scaling (mds)
###################################################################
echo "    stat.t08...\c"
echo "sample1	sample2	distance
1094PC0013	1094PC0013	0.0
1094PC0013	1094PC0012	0.0
1094PC0009	1094PC0009	0.0
1094PC0005	1094PC0013	0.25
1094PC0013	1094PC0005	0.25
1094PC0005	1094PC0012	0.25
1094PC0012	1094PC0013	0.0
1094PC0013	1094PC0009	0.4
1094PC0012	1094PC0012	0.0
1094PC0009	1094PC0005	0.25
1094PC0009	1094PC0013	0.4
1094PC0005	1094PC0005	0.0
1094PC0012	1094PC0009	0.4
1094PC0009	1094PC0012	0.4
1094PC0005	1094PC0009	0.25
1094PC0012	1094PC0005	0.25" > exp
oncogemini stats --mds test.snpeff.vcf.db \
      > obs
check obs exp
rm obs exp

###################################################################
# 9. Test summarize
###################################################################
echo "    stat.t09...\c"
echo "sample	total	num_het	num_hom_alt	num_hom_ref
1094PC0005	1	0	1	4
1094PC0009	1	0	1	4
1094PC0012	5	2	3	3
1094PC0013	4	2	2	3" > exp
oncogemini stats --summarize "select * from variants" test3.snpeff.db > obs
check obs exp
rm obs exp
