from setuptools import setup, find_packages

setup(
    name='snATAC_utilities',
    version='0.1',
    author='Leighton Ditchburn',
    author_email='leighton.ditchburn@outlook.com',
    description='Collection of functions for preprocessing, clustering and plotting of snATAC-seq analysis with snapATAC2',
    packages=find_packages(),
    install_requires=['matplotlib',
                      'seaborn',
                      'pandas',
                      'scanpy',
                      'snapatac2',
                      'pybedtools',
                      'tqdm',
                      'numpy'],
    classifiers=['Programming Language :: Python :: 3',
                 'License :: OSI Approved :: MIT License',
                 'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)