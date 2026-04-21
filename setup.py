from setuptools import setup, find_packages

setup(
    name='cemba_data',
    version='2.0.0',
    author='Hanqing Liu, Wubin Ding',
    author_email='hanliu@salk.edu, wding@salk.edu',
    maintainer='Mykhailo Batiuk',
    maintainer_email='mykhailo.batiuk@protonmail.ch',
    description='Pipelines for single nucleus methylome and multi-omic dataset',
    url='https://github.com/mbatiuk/cemba_data',
    license='MIT',
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
    ],
    packages=find_packages(exclude=('doc',)),
    include_package_data=True,
    package_data={
        '': ['*.txt', '*.tsv', '*.csv', '*.fa', '*Snakefile', '*ipynb', '*.smk']
    },
    install_requires=[],
    entry_points={
        'console_scripts': [
            'yap=cemba_data.__main__:main',
            'yap-internal=cemba_data._yap_internal_cli_:internal_main',
            'yap-hisat3n=cemba_data.hisat3n.cli:main',
            'yap-gcp=cemba_data.gcp:main',
        ]
    }
)