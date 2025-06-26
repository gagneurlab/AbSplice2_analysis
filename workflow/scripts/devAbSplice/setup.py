from setuptools import find_packages, setup


requirements = []


setup(
    name='devAbSplice',
    packages=find_packages(),
    version='0.1.0',
    description='Scripts related with devAbSplice ',
    author='Nils Wagner',
    license='MIT',
    test_suite='tests',
    tests_require=['pytest'],
    install_requires=requirements
)
