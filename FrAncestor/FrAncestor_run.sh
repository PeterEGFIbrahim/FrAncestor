#!/usr/bin/env bash
#$ -pe smp 40
#$ -jc short
#$ -cwd
#$ -p -100
#$ -q ddu.q@hpc-ddu-gpu-001.compute.dundee.ac.uk

source /cluster/ddu/2448959/miniconda3/bin/activate
conda activate FMO_JCPE

# -----------------------
# CONFIGURATION VARIABLES
# -----------------------

num_cpus=40                  							# -c
input_smiles="."         	 							# -smi
run_0D=true                  							# -0D (flag only)
run_1D=true                  							# -1D (flag only)
Ro3_Ro5=3             							        # -r
halogen="0,3"                							# -Halogen
run_2D=true                  							# -2D
run_3D=10                    							# -3D <value>
run_4D=true                  							# -4D
run_5D=true                  							# -5D
target_hotspots="../targets_FMOPhore_3D/target_hotspots_FMOPhore.txt"  # -5D file
min_matches=3               							# -min
Threshold=1.1             							    # -T
run_6D=true                 							# -6D
reactions_lib="../Reactions_lib/Reactions_lib_std.smi"   		# -R_lib
rxs=1                       							# -RxS
chemicalspace_file="chemialspaces/"  # -cs
sss_threshold=0.1           							# -sss
chunk_size=5000              							# -s
fingerprints=(maccs 0.5 pharm2d 0.5)   					# -fp

# -----------------------
# BUILD COMMAND
# -----------------------

cmd="FrAncestor -c ${num_cpus} -dir ${input_smiles}"

$run_0D && cmd+=" -0D"
$run_1D && cmd+=" -1D"
cmd+=" -r ${Ro3_Ro5} -Halogen ${halogen}"
$run_2D && cmd+=" -2D"
[ -n "$run_3D" ] && cmd+=" -3D ${run_3D}"
$run_4D && cmd+=" -4D"
$run_5D && cmd+=" -5D ${target_hotspots}"
cmd+=" -min ${min_matches} -T ${Threshold}"
$run_6D && cmd+=" -6D ${reactions_lib} -RxS ${rxs}"
cmd+=" -cs ${chemicalspace_file} -sss ${sss_threshold} -s ${chunk_size}"

# Append fingerprint arguments
if [ ${#fingerprints[@]} -gt 0 ]; then
    cmd+=" -fp"
    for fp in "${fingerprints[@]}"; do
        cmd+=" ${fp}"
    done
fi

# -----------------------
# RUN
# -----------------------
echo "Running command:"
echo "$cmd"
eval "$cmd "

conda deactivate
conda deactivate

# FrAncestor -h
# FrAncestor -c 10 -dir ./ -0D -1D -r 3 -Halogen 0,3 -2D -3D 10 -4D -5D targets/target_hotspots_FMOPhore.txt -min 3 -T 1.1 -6D -R_lib Reactions/Reactions_lib_std.smi -RxS 1 -cs ChemicalSpace/Francestor_filtered_797_std.smi -sss 0.1 -s 200 -fp maccs 0.5 pharm2d 0.5

# FrAncestor 6-Dimensional Fragment Screening Library Design

# optional arguments:
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
#         -r          --Rule_RO3_RO5          : 
#         -MWT        --MWT                   : Rule of Three (RO3) or Lipinski's Rule of Five (RO5) e.g. 3 or 5
#         -nHeavy     --nHeavy                : 
#         -nHBA       --nHBA                  : ------------------------ Descriptor customization -------------------------
#         -nHBD       --nHBD                  : 
#         -nRB        --nRB                   : Any of these descriptors can be specified separately:
#         -Fsp3       --Fsp3                  : 
#         -nChiC      --nChiC                 : It must be provided in a range between two integers separated by , 
#         -nRings     --nRings                :         
#         -nArRings   --nArRings              :       e.g.  -MWT 250,500 
#         -TPSA       --TPSA                  : 
#         -Halogen    --Halogen               : Selects only molecules with MWT between 250 and 500
#         -logP       --logP                  :
#         -logD       --logD                  : Same for other flagged descriptors 
#         -QED        --QED                   :
#         -SFI        --SFI                   : 
#         -SAScore    --SAScore               :  
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -2D         --Second_Dimension      : Scaffold and Rgroup analysis - Generate 2D
#         -fp         --similarity_fp_tree    : Fingerprint names followed by thresholds; pharm2d 0.5 maccs 0.2 morgan 0.1 torsion atompairs or avalon
#                                             : e.g. -fp maccs 0.5 pharm2d 0.5 morgan 0.8 morgan 0.3 ... 
#         -dp         --diversity_picker      : Randomly pick molecules based on Diversity Picker modul. Number of desired molecules. e.g., 1000
#     --------------------------------------------------------------------------------------------------------------------------------------
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
#         -R_lib,     --Rx_lib                : Reagents library (Default Enamine Building blocks .smi file)
#         -RxS        --Rx_Steps              : Number of reaction cycles (Reaction-steps).
#     --------------------------------------------------------------------------------------------------------------------------------------
#     -cs         --chemicalspace         : Path to chemical space molecules file (finding Near-Neighbours) 
#                                             : One or more files used as chemical space reference (must differ from input files)
#         -fp         --similarity_fp_tree    : Fingerprint names followed by thresholds; pharm2d 0.5 maccs 0.2 morgan 0.1 torsion atompairs or avalon
#                                             : e.g. -fp maccs 0.5 pharm2d 0.5 morgan 0.8 morgan 0.3 ... 
#         -sss        --sociability_threshold : Scaffold Sociability Threshold (default = 0.1)
#         -s          --chunk_size            : Chunk size for the number of molecules in the Chemical space, (default = 1000)
#     --------------------------------------------------------------------------------------------------------------------------------------
#     Peter E.G.F. Ibrahim
    