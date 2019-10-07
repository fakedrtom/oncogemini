check()
{
	if diff $1 $2; then
    	echo ok
	else
    	echo fail
	fi
}
export -f check

SCRIPT_PATH=$(which oncogemini 2> /dev/null)

if [ $? -eq 1 ]
then
  SCRIPT_PATH=$(dirname "$(readlink -f "$0")")
  export PATH=$PATH:"${SCRIPT_PATH}/../../bin"
fi

echo "Using oncogemini found at: $SCRIPT_PATH" 1>&2
printf "\n"
echo "Moving to test dir"
cd test
echo "Removing existing DBs"
rm -f ./*.db
printf "\n"

# setup the testing databases from the testing VCF files
set -e
echo "##############################"
echo "##### Setting up the DBs #####"
echo "##############################"
printf "\n"

bash data-setup.sh

## bash test-vep-extra.sh || exit

## bash test-vcf-output.sh || exit

printf "\n"
echo "##############################"
echo "##### Begin test scripts #####"
echo "##############################"
printf "\n"

# Test region
echo "Running test-region.sh..."
bash test-region.sh
printf "\n"

# Test amending the database
echo "Running test-amend.sh..." 
bash test-amend.sh
printf "\n"

# Test query tool
echo "Running test-query.sh..."
bash test-query.sh
printf "\n"

# Test database dumping
echo "Running test-dump.sh..."
bash test-dump.sh
printf "\n"

# Test basic functionality
echo "Running test-columns.sh..."
bash test-columns.sh
printf "\n"

# Test genotype BLOB functionality
echo "Running test-genotypes.sh..."
bash test-genotypes.sh
printf "\n"

# Test EFF string derived elements in INFO column
## bash test-effstring.sh

# Test annotate functionality
echo "Running test-annotate-tool.sh..."
bash test-annotate-tool.sh
printf "\n"
## echo "Running test-annotate-tool-vcf.sh..."
## bash test-annotate-tool-vcf.sh
## printf "\n"

# Test stats tool
echo "Running test-stats.sh..."
bash test-stats.sh
printf "\n"

# Test windower
## bash test-windower.sh

# Test wildcards
echo "Running test-wildcards.sh..."
bash test-wildcards.sh
printf "\n"

# Test ROH
## bash test-roh.sh

# Test fusions
## bash test-fusions.sh

## bash test-multiple-alts.sh

# Test muliple columns
echo "Running test-multi-col.sh..."
bash test-multi-col.sh
printf "\n"

# Test EFF string
echo "Running test-eff.sh..."
bash test-eff.sh
printf "\n"

# Test common oncogemini tool parameters
echo "Running test-common-parameters.sh..."
bash test-common-parameters.sh
printf "\n"

# Test bottleneck
echo "Running test-bottleneck.sh..."
bash test-bottleneck.sh
printf "\n"

# Test loh
echo "Running test-loh.sh..."
bash test-loh.sh
printf "\n"

# Test truncal
echo "Running test-truncal.sh..."
bash test-truncal.sh
printf "\n"

# Test unique
echo "Running test-unique.sh..."
bash test-unique.sh
printf "\n"

# Test somatic tools
echo "Running test-somatic-tools.sh..."
bash test-somatic-tools.sh
printf "\n"

# cleanup
echo "Cleaning up DBs created and other files"
rm ./*.db
rm -Rf *.gts
rm _err
