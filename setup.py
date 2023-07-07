from setuptools import setup, find_packages

import cusum

with open("README.md", 'r') as fh:
    long_description = fh.read()

with open("requirements.txt") as req:
    install_req = req.read().splitlines()

setup(name='cusum',
      version=cusum.__version__,
      description='Full cumulative sum processing for change detection',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='',
      author='Antoine Pfefer',
      author_email='',
      install_requires=install_req,
      python_requires='>=3.9',
      license='MIT',
      packages=find_packages(),
      zip_safe=False)
