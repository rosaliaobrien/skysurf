from setuptools import setup, find_packages

setup(
    name="zodisurf",
    version="0.3.0",  
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19.0",
        "scipy>=1.7.0",
    ],
    # Extras allow users to run: pip install .[fast] or pip install .[tutorial]
    extras_require={
        "fast": ["numba>=0.55.0"],
        "tutorial": ["healpy>=1.15.0", "matplotlib>=3.0.0", "jupyterlab>=3.0.0", "numba>=0.55.0"],
        "all": ["healpy>=1.15.0", "matplotlib>=3.0.0", "jupyterlab>=3.0.0", "numba>=0.55.0", "pytest"],
    },
    author="Rosalia O'Brien, Tejovrash Acharya",
    author_email="rosalia.d.obrien@gmail.com", #    description="Python implementation of Kelsall+1998 zodiacal light model with SKYSURF enhancements",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/rosaliaobrien/skysurf/tree/main/Zodiacal_Light_Model_Python",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.7",
    package_data={
        'zodisurf': ['*.json', '*.tab'], 
    },
)