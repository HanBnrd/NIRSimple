import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nirsimple",
    version="0.1.6",
    author="Johann Benerradi",
    author_email="johann.benerradi@gmail.com",
    description="fNIRS signal processing simplified",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/HanBnrd/NIRSimple",
    license='MIT',
    packages=setuptools.find_packages(),
    package_data={"nirsimple": ["tables/*.csv"]},
    install_requires=[
        "numpy",
        "pandas",
        "scipy"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
