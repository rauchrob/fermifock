from setuptools import setup

setup(name='fermifock',
      version='0.1',
      description='library for symbolic computations with ferminonic fock space',
      long_description=open('README.md').read(),
      author='Robert Rauch',
      author_email='mail@robertrauch.de',
      license='MIT',
      install_requires=['sympy'],
      extras_require={'dev': ['nose']},
      test_suite='nose.collector',
      tests_require=['nose'],
      packages=['fermifock'],
      zip_safe=False)
