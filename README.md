# FrAncestor V.0.1 Fragment Screening Library Design

---

# ![Screenshot 2025-05-26 at 10 11 04](https://github.com/user-attachments/assets/0dea040b-d8e8-4ed5-9786-af1e9999b8f0)


---

# FrAncestor V.0.1 Fragment Screening Library Design

FrAncestor is a Python package for designing fragment screening libraries through a comprehensive, six-dimensional analysis framework. 
Each dimension captures a distinct layer of molecular insight: 
## Features
- **1st-D**: Calculates sixteen molecular descriptors for property profiling.
- **2nd-D**: Identifies molecular scaffolds, frameworks, and evaluates scaffold sociability via a Molecular Genealogy Tree (MGT).
- **3rd-D**: Generates 3D conformations for structure-based virtual screening.
- **4th-D**: Performs ligand-pharmacophore virtual screening.
- **5th-D**: Enables ligand–target virtual screening using target hotspot information.
- **6th-D**: Applies reaction vector-based enumeration for synthesizable molecule generation.
- **Molecular Genalogy Tree (MGT)**: FrAncestor also supports **Near-Neighbour** Search for rapid chemical space navigation—enabling fast identification of structurally similar compounds within large Chemical Spaces.
FrAncestor offers a modular, multi-layered approach to accelerate early-stage drug discovery and fragment library design.
---

## System Requirements

- **Dependencies**  

- This package requires the following to be installed before use:

  - RDKit (for molecule parsing and manipulation) - (conda install -c rdkit rdkit)

  - pmapper (for pharmacophore modeling) (Polishchuk P, Kutlushina A, Bashirova D, Mokshyna O, Madzhidov T. Virtual Screening Using Pharmacophore Models Retrieved from Molecular Dynamic Simulations. Int J Mol Sci. 2019;20(23):5834. Published 2019 Nov 20. doi:10.3390/ijms20235834)

- **Supported File Formats**  

  - .smi — SMILES format

  - .sdf — Structure Data File

  - .pdb — Protein Data Bank format 
---

## Installation

### Prerequisites

- Python >= 3.6

### Set up Conda Environment

To ensure all dependencies are correctly installed and avoid conflicts follow these steps:

1. **Install Conda**: If Conda is not already installed, follow the installation guide [here](https://docs.conda.io/en/latest/miniconda.html).

2. **Create the Environment**:
   ```bash
   conda env create -f FrAncestor_env.yml
   conda activate FrAncestor_env
   ```

### Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/PeterEGFIbrahim/FrAncestor.git

2. Navigate to the directory:
   ```bash
   cd FrAncestor
   
4. Install the package:
   ```bash
   pip install .

## Usage

### Mandatory Parameters:
   ```bash
    --------------------------------------------------------------------------------------------------------------------------------------
         -dir,       --DIRECTORY             : Process all molecules files in a directory
         -smi,       --SMILES                : library in smiles format .smi
         -sdf,       --SDF                   : library in Structure Data File format .sdf
         -pdb,       --PDB                   : library in Protein Data Bank format .pdb
         -c          --cpus                  : Number of processors to run parallel jobs; if not defined default cpu count.
    --------------------------------------------------------------------------------------------------------------------------------------  
   ```
### Optional Parameters:
   ```bash
#     Dimensions parameters:
     --------------------------------------------------------------------------------------------------------------------------------------
     -0D         --chemspacestandardiser    : Chemical Space Standarisation Processor
         -pains      --PAINS_filter          : Filter out PAINS
     --------------------------------------------------------------------------------------------------------------------------------------
     -1D         --First_Dimension       : Library Descriptors - Generate 1D
         -r          --Rule_RO3_RO5          : 
         -MWT        --MWT                   : Rule of Three (RO3) or Lipinski's Rule of Five (RO5) e.g. 3 or 5
         -nHeavy     --nHeavy                : 
         -nHBA       --nHBA                  : ------------------------ Descriptor customization -------------------------
         -nHBD       --nHBD                  : 
         -nRB        --nRB                   : Any of these descriptors can be specified separately:
         -Fsp3       --Fsp3                  : 
         -nChiC      --nChiC                 : It must be provided in a range between two integers separated by , 
         -nRings     --nRings                :         
         -nArRings   --nArRings              :       e.g.  -MWT 250,500 
         -TPSA       --TPSA                  : 
         -Halogen    --Halogen               : Selects only molecules with MWT between 250 and 500
         -logP       --logP                  :
         -logD       --logD                  : Same for other flagged descriptors 
         -QED        --QED                   :
         -SFI        --SFI                   : 
         -SAScore    --SAScore               :  
     --------------------------------------------------------------------------------------------------------------------------------------
     -2D         --Second_Dimension      : Scaffold and Rgroup analysis - Generate 2D
         -fp         --similarity_fp_tree    : Fingerprint names followed by thresholds; pharm2d 0.5 maccs 0.2 morgan 0.1 torsion atompairs or avalon
                                             : e.g. -fp maccs 0.5 pharm2d 0.5 morgan 0.8 morgan 0.3 ... 
         -dp         --diversity_picker      : Randomly pick molecules based on Diversity Picker modul. Number of desired molecules. e.g., 1000
     --------------------------------------------------------------------------------------------------------------------------------------
     -3D         --Third_Dimension       : Pharmacophoric analysis - Generate 3D
                                             : Number of conformations to be generated per molecule. (defaults: 3)
     --------------------------------------------------------------------------------------------------------------------------------------
     -4D         --Fourth_Dimension      : Pharmacophore Virtual Screening - Generate 4D
                                             : Threshold of RMSD in angstrom (A), in decimals. (Default = 1.1)
     --------------------------------------------------------------------------------------------------------------------------------------
     -5D         --Fifth_Dimension       : Target Virtual Screening - Generate 5D
                                             : Target molecule file to compare against.
         -min        --min_features          : Number of minimum matching pharmacophoric features (defaults: 3)
         -T          --Threshold             : Threshold of RMSD in angstrom (A), in decimals e.g., 1.1
     --------------------------------------------------------------------------------------------------------------------------------------
     -6D         --Sixth_Dimension       : Reaction Vectors and R-group enumeration - Generate 6D 
                                             : Rgroup decomposition in binding site after Hit identification.
         -R_lib,     --Rx_lib                : Reagents library (Default Enamine Building blocks .smi file)
         -RxS        --Rx_Steps              : Number of reaction cycles (Reaction-steps).
     --------------------------------------------------------------------------------------------------------------------------------------
     -cs         --chemicalspace         : Path to chemical space molecules file (finding Near-Neighbours) 
                                             : One or more files used as chemical space reference (must differ from input files)
         -fp         --similarity_fp_tree    : Fingerprint names followed by thresholds; pharm2d 0.5 maccs 0.2 morgan 0.1 torsion atompairs or avalon
                                             : e.g. -fp maccs 0.5 pharm2d 0.5 morgan 0.8 morgan 0.3 ... 
         -sss        --sociability_threshold : Scaffold Sociability Threshold (default = 0.1)
         -s          --chunk_size            : Chunk size for the number of molecules in the Chemical space, (default = 1000)
     --------------------------------------------------------------------------------------------------------------------------------------
     Peter E.G.F. Ibrahim
   ```

### Example Command
   ```
FrAncestor -c 10 -dir ./ -0D -1D -r 3 -Halogen 0,3 -2D -3D 10 -4D -5D targets/target_hotspots_FMOPhore.txt -min 3 -T 1.1 -6D -R_lib Reactions/Reactions_lib_std.smi -RxS 1 -cs ChemicalSpace/Francestor_filtered_797_std.smi -sss 0.1 -s 200 -fp maccs 0.5 pharm2d 0.5
   ```
### Help
To view the full list of options and their usage:
   ```
FrAncestor --help
   ```
### Outputs

After running FrAncestor expected output can be found in the examples folder.

### Requirements

The following Python libraries are required (automatically installed with the package):

- `pmaaper`
- `RDKit`
- `argparse`

### Developer Information

Author: Peter E.G.F. Ibrahim  
Email: pibrahim001@dundee.ac.uk | 2448959@dundee.ac.uk | peteregfi@gmail.com  
GitHub: [PeterEGFIbrahim](https://github.com/PeterEGFIbrahim)  
