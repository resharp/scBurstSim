from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

setup(
    name='scBurstSim',
    version='0.1.0',
    description='Burst simulator',
    long_description=readme,
    author='Richard Gremmen',
    author_email='gremmen@resharp.nl',
    url='https://github.com/resharp/scBurstSim',
    packages=find_packages(exclude=('tests', 'docs'))
)