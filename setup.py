import os
from setuptools import setup
import subprocess as sp

version_py = os.path.join(os.path.dirname(__file__), 'oncogemini', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '').strip()
long_description = """
``oncogemini`` is a database framework for exploring genetic variation in tumor samples'
"""

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines() if not
                        x.startswith(("oncogemini", "http", "git", "vcfanno", "vcf2db"))]
#with open("requirements.txt", "r") as f:
#    conda_install_requires = [x.strip() for x in f.readlines() if 
#                        x.startswith(("vcfanno", "vcf2db"))]

setup(
        name="oncogemini",
        version=version,
        install_requires=install_requires,
        requires=['python (>=2.5, <3.0)'],
        packages=['oncogemini',
                  'oncogemini.scripts',
                  'oncogemini.tests',
                  'oncogemini.data'],
        author="Thomas Nicholas",
        description='A database framework for exploring genetic variation in tumor samples',
        long_description=long_description,
        url="http://oncogemini.readthedocs.org",
        package_dir={'oncogemini': "oncogemini"},
        package_data={'oncogemini': [
            'static/css/gemini.css',
            'static/img/gemini.png',
            'static/third_party/bootstrap/css/*',
            'static/third_party/bootstrap/img/*',
            'static/third_party/bootstrap/js/*',
            'static/third_party/jquery/jquery-1.7.1.js',
            'static/third_party/jquery-ui/jquery-ui.min.js',
            'views/*'
            ]},
        zip_safe=False,
        include_package_data=True,
        scripts=['oncogemini/scripts/oncogemini'],
        author_email="thomas.nicholas@utah.edu",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )

#sp.check_call(["conda","install","-y","-c","bioconda", "vcf2db", "vcfanno"])
