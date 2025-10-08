#!/usr/bin/env bash
#$ -pe smp 40
#$ -jc short
#$ -cwd
#$ -p -100

conda activate FrAncestor_env

# FrAncestor -c 10 -dir ./ -0D -1D -r 3 -Halogen 0,3 -2D -3D 10 -4D -5D targets/target_hotspots_FMOPhore.txt -min 3 -T 1.1 -6D -R_lib Reactions/Reactions_lib_std.smi -RxS 1 -cs ChemicalSpace/Francestor_filtered_797_std.smi -sss 0.1 -s 200 -fp maccs 0.5 pharm2d 0.5

# FrAncestor -h
# FrAncestor Tool-kit

# options:
#   -h, --help            show this help message and exit

#     Mandatory parameters:
#    --------------------------------------------------------------------------------------------------------------------------------------
#         -dir,       --DIRECTORY             : Process all molecules files in a directory
#         -smi,       --SMILES                : library in smiles format .smi
#         -sdf,       --SDF                   : library in Structure Data File format .sdf
#         -pdb,       --PDB                   : library in Protein Data Bank format .pdb
#         -c          --cpus                  : Number of processors to run parallel jobs; if not defined default cpu count.
#    --------------------------------------------------------------------------------------------------------------------------------------
# :
#   -These are required parameters THESE ARE REQUIRED PARAMETERS, -! THESE ARE REQUIRED PARAMETERS

#     Dimensions parameters:
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -0D         --chemspacestandardiser    : Chemical Space Standarisation Processor
#         -pains      --PAINS_filter          : Filter out PAINS
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -1D         --First_Dimension       : Library Descriptors - Generate 1D
#         -r          --Rule_RO3_RO5          : Rule of Three (RO3) or Lipinski's Rule of Five (RO5) e.g. 3 or 5
#         -MWT        --MWT                   : 
#         -nHeavy     --nHeavy                : ------------------------ Descriptor customization -------------------------
#         -nHBA       --nHBA                  : 
#         -nHBD       --nHBD                  : Any of these descriptors can be specified separately:
#         -nRB        --nRB                   : 
#         -Fsp3       --Fsp3                  : It must be provided in a range between two integers separated by , 
#         -nChiC      --nChiC                 : 
#         -nRings     --nRings                :         e.g.  -MWT 250,500 
#         -nArRings   --nArRings              :       
#         -TPSA       --TPSA                  : Selects only molecules with MWT between 250 and 500
#         -Halogen    --Halogen               : 
#         -logP       --logP                  : Same for all other flagged descriptors 
#         -logD       --logD                  : 
#         -QED        --QED                   :
#         -SFI        --SFI                   : 
#         -SAScore    --SAScore               :  
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -2D         --Second_Dimension      : Scaffold and Rgroup analysis - Generate 2D
#      --------------------------------------------------------------------------------------------------------------------------------------
#     -3D         --Third_Dimension       : Pharmacophoric analysis - Generate 3D
#                                             : Number of conformations to be generated per molecule. (defaults: 3)
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -4D         --Fourth_Dimension      : Pharmacophore Virtual Screening - Generate 4D
#                                             : Threshold of RMSD in angstrom (A), in decimals. (Default = 1.1)
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -5D         --Fifth_Dimension       : Target Virtual Screening - Generate 5D
#                                             : Target molecule file to compare against.
#         -min        --min_features          : Number of minimum matching pharmacophoric features (defaults: 3)
#         -T          --Threshold             : Threshold of RMSD in angstrom (A), in decimals e.g., 1.1
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -6D         --Sixth_Dimension       : Reaction Vectors and R-group enumeration - Generate 6D 
#                                             : Rgroup decomposition in binding site after Hit identification.
#                                             : Reagents library (Default Enamine Building blocks .smi file)
#         -RxS        --Rx_Steps              : Number of reaction cycles (Reaction-steps).
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -MGT        --chemicalspace         : Path to chemical space molecules file (finding Near-Neighbours) 
#                                             : One or more files used as chemical space reference (must differ from input files)
#         -fp         --similarity_fp_tree    : Fingerprint names followed by thresholds; pharm2d 0.5 maccs 0.2 morgan 0.1 torsion 0.1 atompairs 0.1 or avalon 0.1
#                                             : e.g. -fp maccs 0.5 pharm2d 0.5 morgan 0.1 torsion 0.1 atompairs 0.1 avalon 0.1 ... 
#         -sss        --sociability_threshold : Scaffold Sociability Threshold (default = 0.1)
#         -s          --chunk_size            : Chunk size for the number of molecules in the Chemical space, (default = 1000)
#     --------------------------------------------------------------------------------------------------------------------------------------
#     Peter E.G.F. Ibrahim

conda deactivate
conda deactivate
