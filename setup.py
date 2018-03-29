from setuptools import setup

setup(name='unidip',
      version='0.1',
      description='Python port of the UniDip clustering algorithm',
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