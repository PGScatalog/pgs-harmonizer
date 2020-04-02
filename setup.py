from setuptools import setup

setup(name='pgs-harmonizer',
      version='0.1',
      description='A pipeline to format and harmonize Polygenic Score (PGS) Catalog Scoring Files within and between different genome builds.',
      url='https://github.com/PGScatalog/pgs-harmonizer',
      author='Sam Lambert',
      author_email='pgs-info@ebi.ac.uk',
      license='Apache 2.0',
      packages=['pgs_harmonizer'],
      zip_safe=False)