from jme.jupy_tools.__version__ import VERSION
from setuptools import setup, find_namespace_packages

DESCRIPTION = "Thngs I've found Useful in Jupyter"
NAME = "jupy_tools"
AUTHOR = "John Eppley"
AUTHOR_EMAIL = "jmeppley@gmail.com"
MAINTAINER = "John Eppley"
MAINTAINER_EMAIL = "jmeppley@gmail.com"
URL = 'http://github.com/jmeppley/jupy_tools'
LICENSE = 'MIT'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      license=LICENSE,
      packages=find_namespace_packages(include=['jme.*']),
      install_requires=['pandas', 'matplotlib', 'ipykernel', 'polars'],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT',
          'Natural Language :: English',
          'Programming Language :: Python :: 3.7'],
      )
