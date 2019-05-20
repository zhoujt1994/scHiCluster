from setuptools import setup, find_packages

schicluster_version = '1.1'

setup(
    name='schicluster',
    version=schicluster_version,
    author='Jingtian Zhou',
    author_email='jiz509@eng.ucsd.edu',
    packages=find_packages(),
    description='A package for single-cell Hi-C data clustering and downstream analysis.',
    long_description=open('README.md').read(),
    include_package_data=True,
    install_requires=['numpy', 'scipy', 'scikit-learn']
)
