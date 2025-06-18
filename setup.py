from setuptools import setup, find_packages

setup(
    name="GeneCover",
    version="2.0",
    packages=find_packages(),
    install_requires=[
        "numpy",  # Add your external packages here
        "gurobipy",
        "scipy",
        "pyscipopt"
    ],
    author="An Wang",
    description="GeneCover: A Combinatorial Approach for Pre-labeling Marker Gene Selection",
    url="https://github.com/ANWANGJHU/GeneCover.git",
    license="MIT"
)