from setuptools import setup, find_packages

setup(
    name='snATAC_utilities',
    version='1.0.0',
    author='Leighton Ditchburn',
    author_email='leighton.ditchburn@outlook.com',
    description='Collection of functions for preprocessing, clustering and plotting of snATAC-seq analysis with snapATAC2',
    packages=find_packages(),
    classifiers=['Programming Language :: Python :: 3',
                 'License :: OSI Approved :: MIT License',
                 'Operating System :: OS Independent',
    ],
    python_requires='>=3.8'
)