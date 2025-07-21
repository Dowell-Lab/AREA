from setuptools import setup

setup(
    name='PSEA',
    version='0.1.0',    
    description='A Python package for participant set enrichment anaylsis',
    url='https://github.com/Dowell-Lab/psea',
    author='Mary Ann Allen, Rutendo Siguke',
    author_email='Mary.A.allen@colorado.edu',
    license='GNU AFFERO GENERAL PUBLIC LICENSE',
    packages=['psea'],
    install_requires=['pandas>1.5',                     
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)

