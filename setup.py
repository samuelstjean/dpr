from setuptools import setup, find_packages

setup(
    name='dpr',
    version='0.1',
    author='Samuel St-Jean',
    author_email='samuel@isi.uu.nl',
    packages=find_packages(),
    # scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    url='https://github.com/samuelstjean/dpr',
    license='LICENSE',
    description='Implementation of "Reducing variability in along-tract analysis with diffusion profile realignment".',
    long_description=open('README.md').read(),
    install_requires=['numpy>=1.10',
                      'scipy>=0.19',
                      'matplotlib>=2.0'],
)
