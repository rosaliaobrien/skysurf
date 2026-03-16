from setuptools import setup, find_packages

setup(
    name="zodisurf",
    version="2.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.19.0",
        "scipy>=1.7.0",
    ],
    author="Rosalia O'Brien, Tejovrash Acharya",
    author_email="your.email@domain.com",
    description="Python implementation of Kelsall+1998 zodiacal light model with SKYSURF enhancements",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/ZodiModel",
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
        'zodi_model': ['*.json', '*.tab'],
    },
)
