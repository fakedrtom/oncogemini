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
# 1. Test the TFAM dump
####################################################################
echo "    dump.t01...\c"
echo "1 M10475 None None 1 1
1 M10478 M10475 M10500 2 2
1 M10500 None None 2 2
1 M128215 M10475 M10500 1 1" > exp
oncogemini dump --tfam test4.snpeff.ped.db > obs
check obs exp
rm obs exp
