import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="msbrainpy",
    version="0.0.1",
    author="Abigail S. McGovern",
    author_email="Abigail.McGovern@monash.edu",
    description="Whole mouse brain image processing and analysis - integrated with AMBA",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AbigailMcGovern/msbrainpy",
    packages=setuptools.find_packages('msbrainpy'),
    package_dir={'': 'msbrainpy'},
    py_modules=['chain', 'io', '__init__', 'map', 'teraPrep', 'base'],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU License",
        "Operating System :: OS Independent",
    ],
)