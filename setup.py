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
    url="https://github.com/AbigailMcGovern/sampleproject",
    packages=['msbrainpy', 'msbrainpy.amba', 'msbrainpy.map', 'msbrainpy.quantify', 'msbrainpy.raw'],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU License",
        "Operating System :: OS Independent",
    ],
)