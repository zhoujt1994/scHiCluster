from setuptools import setup, find_packages

setup(
    name='schicluster',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author='Jingtian Zhou',
    author_email='jiz509@eng.ucsd.edu',
    packages=find_packages(),
    description='A package for single-cell Hi-C data clustering and downstream analysis.',
    long_description=open('README.md').read(),
    license='MIT',
    long_description_content_type='text/markdown',
    url='https://github.com/zhoujt1994/scHiCluster',
    include_package_data=True,
    install_requires=['numpy', 'scipy', 'scikit-learn', 'h5py', 'joblib', 'clodius',
                      'tables', 'cooler', 'pandas', 'statsmodels', 'rpy2', 'anndata', 'xarray'],
    package_data={
        '': ['*.txt', '*.tsv', '*.csv', '*.fa', '*Snakefile', '*ipynb', '*R']
    },
    entry_points={
        'console_scripts': ['hicluster=schicluster.__main__:main',
                            'hic-internal=schicluster._hicluster_internal:internal_main'],
    }
)
