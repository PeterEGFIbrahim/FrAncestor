from setuptools import setup, find_packages

setup(
    name="FrAncestor",
    version="0.1",
    author="Peter E.G.F. Ibrahim",
    author_email="pibrahim001@dundee.ac.uk",
    description="FrAncestor for Fragment library Design.",
    url="https://github.com/PeterEGFIbrahim/FrAncestor",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "": ["*.so"],
    },
    install_requires=[
        "numpy",
        "tqdm",
        "timeout-decorator",
        "pandas",
        "seaborn",
        "matplotlib",
        "rdkit"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "FrAncestor=FrAncestor.FrAncestor_run:main",
        ],
    },
)
