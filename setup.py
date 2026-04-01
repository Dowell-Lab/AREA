from setuptools import setup, find_packages

setup(
    name='area',
    version='0.1.0',
    description='AREA — Association by Rank Enrichment Analysis',
    url='https://github.com/Dowell-Lab/psea',
    author='Mary Ann Allen, Rutendo Sigauke',
    author_email='Mary.A.allen@colorado.edu',
    license='GNU AFFERO GENERAL PUBLIC LICENSE',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    python_requires='>=3.9',
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'statsmodels',
        'pyyaml',
    ],
    extras_require={
        'gpu': ['cupy'],
    },
    entry_points={
        'console_scripts': [
            'area=area.cli:main',
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)
