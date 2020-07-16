from setuptools import setup, find_packages

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='dpr',
    version='0.2',
    author='Samuel St-Jean',
    author_email='samuel@isi.uu.nl',
    packages=find_packages(),
    scripts=['scripts/dpr', 'scripts/dpr_make_fancy_graph'],
    url='https://github.com/samuelstjean/dpr',
    license='MIT',
    description='Implementation of "Reducing variability in along-tract analysis with diffusion profile realignment".',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=['numpy>=1.10',
                      'scipy>=0.19',
                      'matplotlib>=2.0'],
)
