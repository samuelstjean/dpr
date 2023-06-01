from setuptools import setup, find_packages

scripts=['scripts/dpr',
         'scripts/dpr_make_fancy_graph']

setup(scripts=scripts,
      packages=find_packages())
