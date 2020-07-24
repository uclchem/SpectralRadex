import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="example-pkg-YOUR-USERNAME-HERE", # Replace with your own username
    version="0.0.1",
    author="Jonathan Holdship",
    author_email="jonholdship@gmail.com",
    description="A package for RADEX and spectral modelling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jonholdship/spectralradex",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)