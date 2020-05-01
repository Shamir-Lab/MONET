from setuptools import setup

setup(name='monet',
      version=0.1,
      author='Nimrod Rappoport',
      url='https://github.com/Shamir-Lab/MONET',
      description='python implementation of the MONET algorithm for multi-omic module detection',
      license='GPL-3.0',
      python_requires='>3.0.0',
      install_requires=['numpy', 'pandas', 'networkx']
)
