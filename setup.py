from setuptools import setup, find_packages

setup(
    name="GeneCover",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",  # Add your external packages here
        "gurobipy",
        "scipy"
    ],
    author="An Wang",
    description="GeneCover: A Combinatorial Approach for Pre-labeling Marker Gene Selection",
    url="https://github.com/ANWANGJHU/GeneCover.git",
)