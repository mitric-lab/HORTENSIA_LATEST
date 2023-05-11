from setuptools import setup
from setuptools.extension import Extension

setup(name='hortensia_latest',
    version='1.0.2305.1',
    description='Autoionization dynamics',
    url='https://github.com/mitric-lab/HORTENSIA_latest',
    author='Kevin Issler',
    author_email='kevin.issler@uni-wuerzburg.de',
    license='MIT',
    packages=['hortensia_latest','hortensia_latest.wigner','hortensia_latest.gui'],
    package_data={'': ['*.pyx','*.so','*.c','*.png','basis','*.in']},
    python_requires='>=3.8',
    install_requires=[
        "scipy","joblib","cython","pyscf","matplotlib"],
    scripts=['bin/hortensia'],
    include_package_data=True,
    zip_safe=False)
