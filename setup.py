from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='pyspectrafuse',
      version='0.0.1',
      description='Python tools for handle cluster pipeline',
      url='https://github.com/bigbio/pyspectrafuse',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='PgAtk Team',
      author_email='ypriverol@gmail.com',
      license='LICENSE.txt',
      include_package_data=True,
      install_requires=[
          'Click',
          'numpy',
          'pandas',
          'pathlib',
          'pyarrow',
          'spectrum_utils[iplot]',
      ],
      python_requires=">=3.8",
      scripts=['pyspectrafuse/pyspectrafuse_cli.py'],
      packages=find_packages(),
      entry_points={
          'console_scripts': [
              'pyspectrafuse = pyspectrafuse.pyspectrafuse_cli:main'
          ]},
      package_data={'pypgatk': ['config/*.yaml', 'config/*.json']}, zip_safe=False)

import codecs
import os.path

from setuptools import find_packages
from setuptools import setup

with open("README.md", encoding="UTF-8") as fh:
    long_description = fh.read()


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), "r") as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")


setup(
    name="pyspectrafuse",
    version=get_version("pyspectrafuse/__init__.py"),
    author="BigBio Team",
    author_email="ypriverol@gmail.com",
    description="command utils for cluster",
    long_description_content_type="text/markdown",
    long_description=long_description,
    license="'Apache 2.0",
    data_files=[("", ["LICENSE"])],
    package_data={
        "": ["*.xml"],
    },
    url="https://github.com/bigbio/pyspectrafuse",
    packages=find_packages(),
    install_requires=[
        'Click',
        'numpy',
        'pandas',
        'pathlib',
        'pyarrow',
        'pyteomics'
    ],
    entry_points={
        "console_scripts": ['pyspectrafuse_cli = pyspectrafuse.pyspectrafuse_cli:main']
    },
    platforms=["any"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="python multiomics proteomics quantms reanalysis",
    include_package_data=True,
    python_requires=">=3.8",
)
