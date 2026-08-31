# Copyright 2023 The GenoTools Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================


from setuptools import setup, find_packages

setup(
    name='the_real_genotools', 
    version='2.0.1', 
    packages=find_packages(exclude=["tests", "tests.*"]),
    author='Dan Vitale',
    author_email='d.vitale199@gmail.com',
    description='A collection of tools for genotype quality control and analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/dvitale199/GenoTools',
    license='Apache License 2.0',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    extras_require={
        'dev': [
            'pytest>=3.7',
            'twine>=1.11.0'
        ]
    },
    entry_points={
        'console_scripts': [
            'genotools=genotools.cli:main',
            'genotools-download=genotools.download_refs:handle_download',
        ]
    },
    install_requires = [
        'dash_bio>=1.0.2',
        'joblib>=1.3.0',
        'matplotlib>=3.6.2',
        'numba>=0.57.1',
        'numpy>=1.23.5',
        'pandas>=2.0.3',
        'Pillow>=9.3.0',
        'plotly>=5.11.0',
        'psutil>=5.9.0',
        'requests>=2.28.1',
        'scikit_learn>=1.3.0',
        'scipy>=1.9.3',
        'seaborn>=0.12.1',
        'setuptools>=65.6.3',
        'statsmodels>=0.13.5',
        'streamlit>=1.15.2',
        'umap_learn==0.5.3',
        'xgboost>=1.7.6'
    ],
    # container/ is a Docker build context kept as a historical reference, not
    # library code: its Dockerfile COPYs the build context in, so nothing reads
    # these from an installed package. The two *.pkl models it holds are
    # umap_linearsvc pipelines from the 1.x era that 2.0 cannot load at all
    # (AncestryModel.load rejects them), and they were 2.16 MB of a 2.57 MB
    # wheel. Ship the text files, not the models.
    package_data={
      'genotools': ['container/*.txt','container/Dockerfile']
   }
)
