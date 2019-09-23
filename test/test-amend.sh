check() 
{
    if diff <( sort "$1" ) <( sort "$2" ); then
        echo ok
    else
        echo fail
	exit 1
    fi
}

####################################################################
# 1. Test amending the samples table
####################################################################
echo "    amend.t01...\c"
echo "1	1	M10475	0	0	1	1	None	brown
2	1	M10478	M10475	M10500	2	2	None	red
3	1	M10500	0	0	2	2	None	purple
4	1	M128215	M10475	M10500	1	1	None	green" > exp
oncogemini amend --sample test_amend_sample.ped test.amend.db
oncogemini query -q "select * from samples" test.amend.db > obs
check obs exp
rm obs

echo "    amend.t02...\c"
# unchanged
oncogemini amend --sample test_amend2_sample.ped test.amend.db
oncogemini query -q "select * from samples" test.amend.db > obs
check obs exp

echo "    amend.t03...\c"
echo "1	1	M10475	0	0	1	1	None	brown
2	1	M10478	M10475	M10500	2	2	None	red
3	1	M10500	0	0	2	2	None	None
4	1	M128215	M10475	M10500	1	1	None	green" > exp
# unless we use --clear
oncogemini amend --clear --sample test_amend2_sample.ped test.amend.db
oncogemini query -q "select * from samples" test.amend.db > obs
check obs exp
