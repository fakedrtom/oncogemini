set -eo pipefail

echo "##### creating oncogemini_test.db #####"
vcf2db.py test.oncogemini.vcf oncogemini_test.manifest oncogemini_test.db 
echo "##### creating somatic_tests_tool.db #####"
cp oncogemini_test.db somatic_tests_tools.db
echo "##### creating test.vep.extra.db #####"
vcf2db.py test-vep-extra.vcf test-vep-extra.ped test.vep.extra.db 
echo "##### creating test.eff.db #####"
vcf2db.py test.eff.vcf test.eff.ped test.eff.db 
echo "##### creating test.mad.db #####"
vcf2db.py test.multiple-alts.decomp.snpeff.vcf test.mad.ped test.mad.db 
echo "##### creating test.roh.vcf.db #####"
vcf2db.py test.roh.vcf test.roh.ped test.roh.vcf.db 
echo "##### creating test.snpeff.vcf.db #####"
vcf2db.py test.snpeff.vcf test.snpeff.ped test.snpeff.vcf.db
echo "##### creating test1.snpeff.db #####"
vcf2db.py test1.snpeff.vcf test1.snpeff.ped test1.snpeff.db
echo "##### creating test1.vep.db #####"
vcf2db.py test1.vep.vcf test1.snpeff.ped test1.vep.db
echo "##### creating test2.snpeff.db.db #####"
vcf2db.py test2.snpeff.vcf test2.snpeff.ped test2.snpeff.db
echo "##### creating test3.snpeff.db.db #####"
vcf2db.py test3.snpeff.vcf test3.snpeff.ped test3.snpeff.db
echo "##### creating test.query.db #####"
vcf2db.py test.query.no_vep.vcf test.query.ped test.query.db
echo "##### creating test.region.db #####"
vcf2db.py test.region.vep.vcf test_extended_ped.ped test.region.db 
echo "##### creating test4.snpeff.ped.db #####"
vcf2db.py test4.snpeff.vcf test4.snpeff.ped test4.snpeff.ped.db 
echo "##### creating test.vcf_id.snpeff.vcf.db #####"
vcf2db.py test.vcf_id.snpeff.vcf test.vcf_id.snpeff.ped test.vcf_id.snpeff.vcf.db 
echo "##### creating test.family.db #####"
vcf2db.py test.family.vcf test.de_novo.ped test.family.db 
echo "##### creating extended_ped.db #####"
vcf2db.py test4.snpeff.vcf test_extended_ped.ped extended_ped.db 
echo "##### creating test.amend.db #####"
cp extended_ped.db test.amend.db
echo "##### creating test.fusions.db #####"
vcf2db.py test.fusions.vcf test.fusions.ped test.fusions.db 
