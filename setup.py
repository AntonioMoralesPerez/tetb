import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Transversality-Enforced Tight-Binding model",  # This is the name of the package
    version="0.1",  # The initial release version
    author="Antonio Morales-Pérez, Chiara Devescovi, Yoonseok Hwang, Mikel García-Díez, Barry Bradlyn, Juan L. Mañes, Maia G. Vergniory, and Aitzol García-Etxarri",  # Full name of the author
    author_email="antonio.morales@dipc.org",  # Email of the author at the time of first release
    keywords="tight binding, topological quantum chemistry, topological photonics, topological materials",
    description="Python package for obtaining Transversality-Enforced Tight-Binding models in photonic crystals",
    long_description=long_description,  # Long description read from the the readme file
    long_description_content_type="text/markdown",
    py_modules=[
        "tetb",
    ],  # Name of the python package
    package_dir={"": "src"},  # Directory of the source code of the package
    install_requires=[
        "numpy",
        "sympy==1.12",
        "ortools==9.10.4067",
    ],  # Install other dependencies if any
)
