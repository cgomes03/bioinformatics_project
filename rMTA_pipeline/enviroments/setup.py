from setuptools import setup, find_packages

setup(
    name='MTA',
    version='0.1',
    package_dir= {'': 'src'},
    packages=find_packages(where='src', exclude=['tests']),
    install_requires=[ 'scipy',
                    'cobra',
                    'numpy',
                    'pandas<=2.2.3',
                    'mewpy',
                    'mat4py',
                    'pybiomart',
                    'gurobipy<=11.0.3',
                    'mat4py',
                    'deap',
                    'seaborn'

                       ],
     python_requires='>3.7, <3.12',
    #

)
