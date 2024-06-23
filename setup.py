from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="Genomic Insights",
    version="0.1.0",
    author="yehonatanke",
    description="A DNA sequence analysis tool providing insights into sequence characteristics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yehonatanke/Genomic_Insights",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.7",
    install_requires=[
        "biopython==1.79",
        "matplotlib==3.4.3",
        "numpy==1.21.2",
        "scikit-learn==0.24.2",
        "seaborn==0.11.2",
    ],
    entry_points={
        "console_scripts": [
            "genomic_insights=scripts.run_analysis:main",
        ],
    },
)
