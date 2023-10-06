from setuptools import setup, find_packages

setup(
    name='genotools', 
    version='1.0', 
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'genotools=GenoTools.genotools:main',
        ]
    },
    install_requires = [
        'dash_bio==1.0.2',
        'joblib==1.3.0',
        'matplotlib==3.6.2',
        'numba==0.57.1',
        'numpy==1.23.5',
        'pandas==2.0.3',
        'Pillow==9.3.0',
        'plotly==5.11.0',
        'requests==2.28.1',
        'scikit_learn==1.3.0',
        'scipy==1.9.3',
        'seaborn==0.12.1',
        'setuptools==65.6.3',
        'streamlit==1.15.2',
        'umap_learn==0.5.3',
        'xgboost==1.7.6'
    ]
)