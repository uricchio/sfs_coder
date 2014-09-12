from distutils.core import setup

setup(name='sfscoder',
      version='0.1',
      description='Front end to SFS_CODE',
      author='Lawrence Uricchio',
      author_email='uricchil@gmail.com',
      url='http://uricchio.github.io/sfs_coder',
      packages=['sfscoder'],
      package_dir={'sfscoder': 'src/sfscoder'},
      package_data={'sfscoder': ['req/*.pl','req/hg19/*.gz','req/hg19/geneticMap/*.gz']},
     )
