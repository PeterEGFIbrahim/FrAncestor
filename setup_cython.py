from setuptools import setup, Extension
from Cython.Build import cythonize
import os

# List of Python modules to compile
modules = [
    "FrAncestor/FrAncestor_run.py",
    "FrAncestor/modules/D0_FrAncestor_PAINS.py",
    "FrAncestor/modules/D0_FrAncestor_Standardizer.py",
    "FrAncestor/modules/D1_FrAncestor_descs_filter.py",
    "FrAncestor/modules/D1_FrAncestor_descs.py",
    "FrAncestor/modules/D2_FrAncestor_scaffolds_frags.py",
    "FrAncestor/modules/D2_FrAncestor_diversity_picker.py",
    "FrAncestor/modules/D2_FrAncestor_Ph4s.py",
    "FrAncestor/modules/D2_FrAncestor_sss.py",
    "FrAncestor/modules/D2_FrAncestor_tree.py",
    "FrAncestor/modules/D3_FrAncestor_conformers.py",
    "FrAncestor/modules/D3_FrAncestor_Ph4.py",
    "FrAncestor/modules/D4_FrAncestor_ligs.py",
    "FrAncestor/modules/D5_FrAncestor_Ph4_targets.py",
    "FrAncestor/modules/D6_FrAncestor_Rx.py",
    "FrAncestor/modules/FrAncestor_MGT.py",
    "FrAncestor/modules/FrAncestor_utility.py"
    ]

# Create a list of Extension objects
extensions = [
    Extension(
        module.replace("/", ".").replace(".py", ""),
        [module],
    )
    for module in modules
]

setup(
    name="FrAncestor",
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"}),
)