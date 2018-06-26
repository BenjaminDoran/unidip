from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='unidip',
      version='0.1.1',
      description='Python port of the UniDip clustering algorithm',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/BenjaminDoran/unidip',
      author='Benjamin Doran',
      author_email='benadoran@gmail.com',
      license='GPL',
      packages=['unidip'],
      install_requires=[
          'numpy',
          'matplotlib'
      ],
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'],
     )