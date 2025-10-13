import sys
import os
import argparse
import subprocess
import shutil, os, glob
import time
import threading
from multiprocessing import Pool
import multiprocessing
import getpass
import timeout_decorator
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import cycle
from pathlib import Path
import itertools
import socket
import warnings
import os, sys, logging, warnings
from rdkit import RDLogger
from FrAncestor import PAINSFilter
from FrAncestor import MoleculeStandardizer
from FrAncestor import DescFilters
from FrAncestor import FrAncestor1Descs
from FrAncestor import FrAncestor2Dscaffolds
from FrAncestor import ScaffoldSociabilityScore
from FrAncestor import FrAncestor2Dph4s
from FrAncestor import FrAncestor2DDiversityPicker
from FrAncestor import FrAncestor2DTree
from FrAncestor import FrAncestor3DConf
from FrAncestor import FrAncestor3Dph4s
from FrAncestor import FrAncestor4Dligs
from FrAncestor import FrAncestor5Dhotspots
from FrAncestor import FrAncestor6Drxs
from FrAncestor import FrAncestorMGT
from FrAncestor import EnvironmentGuard

###################################################################################################
# Parse command-line arguments
###################################################################################################
mandatory_params_text = """
    \033[1mMandatory parameters:\033[0m
   --------------------------------------------------------------------------------------------------------------------------------------
        -dir,       --DIRECTORY             : Process all molecules files in a directory
        -smi,       --SMILES                : library in smiles format .smi
        -sdf,       --SDF                   : library in Structure Data File format .sdf
        -pdb,       --PDB                   : library in Protein Data Bank format .pdb
        -c          --cpus                  : Number of processors to run parallel jobs; if not defined default cpu count.
   --------------------------------------------------------------------------------------------------------------------------------------
"""
optional_params_text = """
    \033[1mDimensions parameters:\033[0m
    --------------------------------------------------------------------------------------------------------------------------------------
    -0D         --chemspacestandardiser    : \033[1mChemical Space Standarisation Processor\033[0m
        -pains      --PAINS_filter          : Filter out PAINS
    --------------------------------------------------------------------------------------------------------------------------------------
    -1D         --First_Dimension       : \033[1mLibrary Descriptors - Generate 1D\033[0m
        -r          --Rule_RO3_RO5          : Rule of Three (RO3) or Lipinski's Rule of Five (RO5) e.g. 3 or 5
        -MWT        --MWT                   : 
        -nHeavy     --nHeavy                : ------------------------ \033[1mDescriptor customization\033[0m -------------------------
        -nHBA       --nHBA                  : 
        -nHBD       --nHBD                  : Any of these descriptors can be specified separately:
        -nRB        --nRB                   : 
        -Fsp3       --Fsp3                  : It must be provided in a range between two integers separated by , 
        -nChiC      --nChiC                 : 
        -nRings     --nRings                :         e.g.  -MWT 250,500 
        -nArRings   --nArRings              :       
        -TPSA       --TPSA                  : Selects only molecules with MWT between 250 and 500
        -Halogen    --Halogen               : 
        -logP       --logP                  : Same for all other flagged descriptors 
        -logD       --logD                  : 
        -QED        --QED                   :
        -SFI        --SFI                   : 
        -SAScore    --SAScore               :  
    --------------------------------------------------------------------------------------------------------------------------------------
    -2D         --Second_Dimension      : \033[1mScaffold and Rgroup analysis - Generate 2D\033[0m
     --------------------------------------------------------------------------------------------------------------------------------------
    -3D         --Third_Dimension       : \033[1mPharmacophoric analysis - Generate 3D\033[0m
                                            : Number of conformations to be generated per molecule. (defaults: 3)
    --------------------------------------------------------------------------------------------------------------------------------------
    -4D         --Fourth_Dimension      : \033[1mPharmacophore Virtual Screening - Generate 4D\033[0m
                                            : Threshold of RMSD in angstrom (A), in decimals. (Default = 1.1)
    --------------------------------------------------------------------------------------------------------------------------------------
    -5D         --Fifth_Dimension       : \033[1mTarget Virtual Screening - Generate 5D\033[0m
                                            : Target molecule file to compare against.
        -min        --min_features          : Number of minimum matching pharmacophoric features (defaults: 3)
        -T          --Threshold             : Threshold of RMSD in angstrom (A), in decimals e.g., 1.1
    --------------------------------------------------------------------------------------------------------------------------------------
    -6D         --Sixth_Dimension       : \033[1mReaction Vectors and R-group enumeration - Generate 6D\033[0m 
                                            : Rgroup decomposition in binding site after Hit identification.
                                            : Reagents library (Default Enamine Building blocks .smi file)
        -RxS        --Rx_Steps              : Number of reaction cycles (Reaction-steps).
    --------------------------------------------------------------------------------------------------------------------------------------
    -MGT        --chemicalspace         : \033[1mPath to chemical space molecules file (finding Near-Neighbours)\033[0m 
                                            : One or more files used as chemical space reference (must differ from input files)
        -fp         --similarity_fp_tree    : Fingerprint names followed by thresholds; pharm2d 0.5 maccs 0.2 morgan 0.1 torsion 0.1 atompairs 0.1 or avalon 0.1
                                            : e.g. -fp maccs 0.5 pharm2d 0.5 morgan 0.1 torsion 0.1 atompairs 0.1 avalon 0.1 ... 
        -sss        --sociability_threshold : Scaffold Sociability Threshold (default = 0.1)
        -s          --chunk_size            : Chunk size for the number of molecules in the Chemical space, (default = 1000)
    --------------------------------------------------------------------------------------------------------------------------------------
    \033[1mPeter E.G.F. Ibrahim\033[0m
"""
###################################################################################################
parser = argparse.ArgumentParser(description="""\033[1mFrAncestor 6-Dimensional Fragment Screening Library Design\033[0m""")
mandatory_group = parser.add_argument_group(mandatory_params_text)
                            ############################
mandatory_group.add_argument('-These are required parameters', '-!', type=str, help='')
mandatory_group.add_argument('-dir', '--DIRECTORY', type=str, help=argparse.SUPPRESS)
mandatory_group.add_argument('-smi', '--SMILES', type=str, help=argparse.SUPPRESS)
mandatory_group.add_argument('-sdf', '--SDF', type=str, help=argparse.SUPPRESS)
mandatory_group.add_argument('-pdb', '--PDB', type=str, help=argparse.SUPPRESS)
mandatory_group.add_argument('-c', '--cpus', type=int, help=argparse.SUPPRESS)
mandatory_group.add_argument('-split', '--split_dataset', type=int, help=argparse.SUPPRESS)
optional_group = parser.add_argument_group(optional_params_text)
                            ############################
optional_group.add_argument('-These are optinal parameters', '-...', type=str, help='')
optional_group.add_argument('-0D', '--standardise', action='store_true', help=argparse.SUPPRESS)
optional_group.add_argument('-pains', '--PAINS_filter', action='store_true', help=argparse.SUPPRESS)
                            ############################
optional_group.add_argument('-1D', '--First_Dimension', action='store_true', help=argparse.SUPPRESS)
optional_group.add_argument('-r', '--Rule_RO3_RO5', type=int, help=argparse.SUPPRESS)
optional_group.add_argument('-MWT', help=argparse.SUPPRESS)
optional_group.add_argument('-nHeavy', help=argparse.SUPPRESS)
optional_group.add_argument('-nHBA', help=argparse.SUPPRESS)
optional_group.add_argument('-nHBD', help=argparse.SUPPRESS)
optional_group.add_argument('-nRB', help=argparse.SUPPRESS)
optional_group.add_argument('-Fsp3', help=argparse.SUPPRESS)
optional_group.add_argument('-nChiC', help=argparse.SUPPRESS)
optional_group.add_argument('-nRings', help=argparse.SUPPRESS)
optional_group.add_argument('-nArRings', help=argparse.SUPPRESS)
optional_group.add_argument('-TPSA', help=argparse.SUPPRESS)
optional_group.add_argument('-Halogen', help=argparse.SUPPRESS)
optional_group.add_argument('-logP', help=argparse.SUPPRESS)
optional_group.add_argument('-logD', help=argparse.SUPPRESS)
optional_group.add_argument('-QED', help=argparse.SUPPRESS)
optional_group.add_argument('-SFI', help=argparse.SUPPRESS)
optional_group.add_argument('-SAScore', help=argparse.SUPPRESS)
                            ############################
optional_group.add_argument('-2D', '--Second_Dimension', action='store_true', help=argparse.SUPPRESS)
optional_group.add_argument('-dp', type=int, help=argparse.SUPPRESS)
                            ############################
optional_group.add_argument('-3D', '--Third_Dimension', type=int, nargs='?', const=3, help=argparse.SUPPRESS)
                            ############################
optional_group.add_argument('-4D', '--Fourth_Dimension', type=float, nargs='?', const=1.1, help=argparse.SUPPRESS)
                            ############################
optional_group.add_argument('-5D', '--Fifth_Dimension', type=str, help=argparse.SUPPRESS)
optional_group.add_argument('-min', '--min_features', type=int, default=None, help=argparse.SUPPRESS)
optional_group.add_argument('-T', '--Threshold', type=float, default=None, help=argparse.SUPPRESS)
                            ############################
optional_group.add_argument('-6D', '--Sixth_Dimension', type=str, help=argparse.SUPPRESS)
optional_group.add_argument('-RxS', '--Rx_Steps', type=int, help=argparse.SUPPRESS)
                            ############################
optional_group.add_argument("-MGT", nargs='+', metavar='CHEMSPACE', help=argparse.SUPPRESS)
optional_group.add_argument('-fp', nargs='+', help=argparse.SUPPRESS)
optional_group.add_argument("-sss", '--sociability_threshold', type=float, default=None, help=argparse.SUPPRESS)
optional_group.add_argument("-s", '--chunk_size', type=int, default=None, help=argparse.SUPPRESS)
args = parser.parse_args() 
###################################################################################################
if (args.SMILES or args.SDF or args.DIRECTORY or args.PDB) and args.cpus:
    pass
else:
    parser.print_help()
    sys.exit(1)
FrAncestor_dir = os.path.dirname(os.path.abspath(__file__))

# Silence the Boost.Python converter RuntimeWarning
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r"to-Python converter for boost::shared_ptr<RDKit::FilterHierarchyMatcher> already registered.*"
)
RDLogger.DisableLog('rdApp.*')  # or: RDLogger.logger().setLevel(RDLogger.CRITICAL)

# (Optional) If you still see Python-level DeprecationWarnings from other libs:
warnings.filterwarnings("ignore", category=DeprecationWarning)




# 1) Silence RDKit logs like: [09:49:03] DEPRECATION WARNING...
RDLogger.DisableLog('rdApp.*')                 # or: RDLogger.logger().setLevel(RDLogger.CRITICAL)

# 2) Silence Python warnings (tune patterns as needed)
warnings.filterwarnings(
    "ignore",
    category=RuntimeWarning,
    message=r"to-Python converter for boost::shared_ptr<RDKit::FilterHierarchyMatcher> already registered.*"
)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# 3) Quiet logging (you control what you emit)
logging.basicConfig(level=logging.ERROR, format='%(message)s')
###################################################################################################
#  Definition of Dynamic Six Dimensional Fragment Screening Library (D6DFSL)
###################################################################################################
def run_command(command):
    return subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
def renumber_atoms(pdb_file):
    updated_lines = []
    atom_counters = {}
    with open(pdb_file, 'r') as infile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                residue_name = line[17:20].strip()
                residue_number = line[22:26].strip()
                chain_id = line[21].strip()
                atom_name = line[12:16].strip()
                residue_key = (residue_name, residue_number, chain_id)
                if residue_key not in atom_counters:
                    atom_counters[residue_key] = {}
                if atom_name.startswith("HX"):
                    base_name = "H"
                elif atom_name[0].isdigit():
                    base_name = ''.join([char for char in atom_name[1:] if not char.isdigit()])
                else:
                    base_name = ''.join([char for char in atom_name if not char.isdigit()])
                if base_name not in atom_counters[residue_key]:
                    atom_counters[residue_key][base_name] = 0
                atom_counters[residue_key][base_name] += 1
                new_atom_name = f"{base_name}{atom_counters[residue_key][base_name]}"
                new_atom_name = new_atom_name[:4].ljust(4)
                updated_line = line[:12] + new_atom_name + line[16:]
                updated_lines.append(updated_line)
            else:
                updated_lines.append(line)
    with open(pdb_file, 'w') as outfile:
        outfile.writelines(updated_lines)
    print(f"Renamed and renumbered PDB file saved as: {pdb_file}")
####################################################################################################
#  0 Dimension
####################################################################################################
def preparelibs(molecules_file, cpus):
    folder_name = os.path.splitext(molecules_file)[0]
    os.makedirs(folder_name, exist_ok=True)
    # shutil.move(molecules_file, folder_name)
    shutil.copy(molecules_file, folder_name)
    os.chdir(folder_name)
    def libProcessor(molecules_file, cpus):
        input_file = os.path.join(folder_name, os.path.basename(molecules_file))
        base_name = os.path.splitext(os.path.basename(molecules_file))[0]
        std_smiles_file = os.path.join(folder_name, f"{base_name}_std.smi")
        d0_dir = os.path.join(folder_name, "D0")
        os.makedirs(d0_dir, exist_ok=True)
        std_smiles_d0_path = os.path.join(d0_dir, f"{base_name}_std.smi")
        if os.path.exists(std_smiles_d0_path):
            log_message(f"Standardized file already exists: {std_smiles_d0_path} — skipping standardization.")
            return 
        if any([
            args.standardise, args.MGT, args.Sixth_Dimension,
            args.Fifth_Dimension, args.Fourth_Dimension,
            args.Third_Dimension, args.Second_Dimension,
            args.First_Dimension
        ]):
            standardizer = MoleculeStandardizer()
            standardizer.standardize_file(input_file, std_smiles_file, cpus)
            if not os.path.isfile(std_smiles_file):
                raise FileNotFoundError(f"Standardized file not created: {std_smiles_file}")
        else:
            shutil.move(input_file, std_smiles_file)
        shutil.move(std_smiles_file, std_smiles_d0_path)
        if args.PAINS_filter:
            PAINS_REOS_removed_file = os.path.join(folder_name, "PAINS_REOS_removed.smi")
            pains_filter = PAINSFilter()
            pains_filter.filter_smiles_file(os.path.join(d0_dir, f"{base_name}_std.smi"), PAINS_REOS_removed_file, cpus)
            shutil.move(PAINS_REOS_removed_file, os.path.join(d0_dir, "PAINS_REOS_removed.smi"))
    try:
        libProcessor(molecules_file, cpus)
    except Exception as e:
        log_message(f"🚫 An error occurred in libProcessor: {str(e)}")
    finally:
        os.chdir("..")
    #############################################################
def preparelibrary_all(libs, cpus):
    files = [lib for lib in libs if os.path.isfile(lib)]
    if not files:
        log_message("No valid molecule files found for 0D processing.")
        return
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = {executor.submit(preparelibs, lib, cpus): lib for lib in files}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log_message(f"🚫 An error occurred in D0_Processor: {str(e)}")
####################################################################################################
#  1 Dimension
####################################################################################################
def D1(molecules_file, cpus, sascore):
    folder_name = os.path.splitext(molecules_file)[0]
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    try:
        ##############################################################
        def update_molecules_file(original_file, filtered_file):
            with open(filtered_file, 'r') as f:
                lines = f.readlines()
            with open(original_file, 'w') as f:
                for line in lines[1:]:
                    columns = line.strip().split('\t')  
                    first_three_columns = columns[:3]
                    f.write('\t'.join(first_three_columns) + '\n') 
            ##############################################################
        def D1_run(molecules_file, cpus, sascore):
            d1_dir = os.path.join(folder_name, "D1")
            os.makedirs(d1_dir, exist_ok=True)
            base_name = os.path.splitext(molecules_file)[0]
            descriptor_txt_name = f"{base_name}_descriptors.txt"
            descriptor_txt_path = os.path.join(d1_dir, descriptor_txt_name)
            std_smi_path = os.path.join(folder_name, "D0", f"{base_name}_std.smi")
            # ====== STEP 1: Generate descriptors if not already existing ======
            if not os.path.exists(descriptor_txt_path):
                log_message(f"Generating descriptors: {descriptor_txt_name}")
                D1_processor = FrAncestor1Descs(sascore_file=sascore)
                D1_processor.calculate1Descs(std_smi_path, descriptor_txt_name, cpus, sascore)
                shutil.move(descriptor_txt_name, descriptor_txt_path)
                for file in glob.glob(os.path.join(folder_name, "*descriptors.png")):
                    shutil.move(file, os.path.join(d1_dir, os.path.basename(file)))
            else:
                log_message(f"Descriptor file already exists: {descriptor_txt_path} — skipping descriptor generation.")
            # ====== STEP 2: Check if filtering is required ======
            filtering_requested = (
                args.Rule_RO3_RO5 or
                any([
                    args.MWT, args.nHeavy, args.nHBA, args.nHBD, args.nRB, args.Fsp3, args.nChiC, 
                    args.nRings, args.nArRings, args.TPSA, args.Halogen, args.logP, args.logD, 
                    args.QED, args.SFI, args.SAScore
                ])
            )
            if filtering_requested:
                log_message("Applying descriptor filtering...")
                filter_criteria = []
                if args.Rule_RO3_RO5:
                    if args.Rule_RO3_RO5 == 3:
                        filter_criteria = ['MWT', '0,300', 'nHeavy', '0,23', 'nRB', '0,3',
                                           'nHBD', '0,3', 'nHBA', '0,3', 'logP', '-3,3', 'TPSA', '0,60']
                    elif args.Rule_RO3_RO5 == 5:
                        filter_criteria = ['MWT', '0,500', 'nHeavy', '0,32', 'nRB', '0,8',
                                           'nHBD', '0,5', 'nHBA', '0,10', 'logP', '-5,5', 'TPSA', '0,140']
                    else:
                        raise ValueError("Unsupported rule specified. Use 3 for RO3 or 5 for RO5.")
                custom_filters = {
                    'MWT': args.MWT, 'nHeavy': args.nHeavy, 'nHBA': args.nHBA, 'nHBD': args.nHBD,
                    'nRB': args.nRB, 'Fsp3': args.Fsp3, 'nChiC': args.nChiC, 'nRings': args.nRings,
                    'nArRings': args.nArRings, 'TPSA': args.TPSA, 'Halogen': args.Halogen,
                    'logP': args.logP, 'logD': args.logD, 'QED': args.QED, 'SFI': args.SFI,
                    'SAScore': args.SAScore
                }
                for desc, val in custom_filters.items():
                    if val:
                        filter_criteria.append(desc)
                        filter_criteria.append(val)
                if filter_criteria:
                    D1_filter_processor = DescFilters(filter_criteria=filter_criteria)
                    D1_filter_processor.processD1filter(descriptor_txt_path, descriptor_txt_path, cpus)
                    update_molecules_file(std_smi_path, descriptor_txt_path)
                    for file in glob.glob(os.path.join(folder_name, "*descriptors_filtered.png")):
                        shutil.move(file, os.path.join(d1_dir, os.path.basename(file)))
                else:
                    log_message("No descriptor filtering criteria provided.")
        ##############################################################
        D1_run(os.path.basename(molecules_file), cpus, sascore)
    except Exception as e:
        log_message(f"🚫 An error occurred in D1: {str(e)}")
    finally:
        os.chdir("..")
    #############################################################
def D1_Processor(libs, cpus, sascore):
    files = [lib for lib in libs if os.path.isfile(lib)]
    if not files:
        log_message("No valid molecule files found for 1D processing.")
        return
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = {executor.submit(D1, lib, cpus, sascore): lib for lib in files}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log_message(f"🚫 An error occurred in D1_Processor: {str(e)}")
####################################################################################################
#  2 Dimension
####################################################################################################
def D2(molecules_file, cpus):
    folder_name = os.path.splitext(molecules_file)[0]
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    try:
        def D2_run(molecules_file, cpus):
            base_name = os.path.splitext(os.path.basename(molecules_file))[0]
            d2_dir = os.path.join(folder_name, "D2")
            os.makedirs(d2_dir, exist_ok=True)
            # Path to expected scaffold output
            scaffold_file = os.path.join(d2_dir, f"{base_name}_std_scaffold.txt")
            if os.path.isfile(scaffold_file):
                log_message(f"ChemicalSpace already prepared... D2 path: {scaffold_file}")
                return 
            # Step 1: 2D Scaffold calculation
            D2_processor = FrAncestor2Dscaffolds()
            D2_processor.calculate2Dscaffolds(f"D0/{molecules_file[:-4]}_std.smi", cpus)
            # Step 2: Pharmacophore features
            os.chdir(d2_dir)
            try:
                molecules_file_ph4 = f"{molecules_file[:-4]}_std_scaffold.txt"
                D2ph4s_processor = FrAncestor2Dph4s()
                D2ph4s_processor.calculate2Dph4s(molecules_file_ph4, cpus)
            finally:
                os.chdir("..")
            # Step 3: Scaffold Sociability Score
            if args.sociability_threshold:
                D2_processor = ScaffoldSociabilityScore()
                D2_processor.calculateSSS(scaffold_file, scaffold_file, cpus)
            # Step 4: Fingerprint-based processing
            if args.fp and not args.MGT:
                if args.dp:
                    D2_fptree_processor = FrAncestor2DDiversityPicker(fingerprints=args.fp)
                    D2_fptree_processor.calculate2DdiversityPicker(molecules_file, molecules_file, args.dp, cpus)
                else:
                    D2_fptree_processor = FrAncestor2DTree(fingerprints=args.fp)
                    D2_fptree_processor.calculate2DTree(molecules_file, cpus)
        D2_run(os.path.basename(molecules_file), cpus)
    except Exception as e:
        log_message(f"🚫 An error occurred in D2: {str(e)}")
    finally:
        os.chdir("..")
    #############################################################
def D2_Processor(libs, cpus):
    files = [lib for lib in libs if os.path.isfile(lib)]
    if not files:
        log_message("No valid molecule files found for 2D processing.")
        return
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = {executor.submit(D2, lib, cpus): lib for lib in files}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log_message(f"🚫 An error occurred in D2_Processor: {str(e)}")
####################################################################################################
#  3 Dimension
####################################################################################################
def D3(molecules_file, cpus):
    folder_name = os.path.splitext(molecules_file)[0]
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    try:
        #############################################################
        def D3_run(molecules_file, cpus):
            base_name = os.path.splitext(os.path.basename(molecules_file))[0]
            d3_dir = os.path.join(folder_name, "D3")
            os.makedirs(d3_dir, exist_ok=True)
            sdf_file = os.path.join(d3_dir, f"{base_name}_std_scaffold.sdf")
            if os.path.isfile(sdf_file):
                log_message(f"ChemicalSpace already in 3D: {sdf_file}")
                return
            nconfs = args.Third_Dimension if getattr(args, 'Third_Dimension', None) is not None else 3
            D3Dconf_processor = FrAncestor3DConf()
            D3Dconf_processor.calculate3Dconf(
                molecules_file=os.path.join("D2", f"{base_name}_std_scaffold.txt"),
                nconfs=nconfs,
                seed=42,
                cpus=cpus
            )
            while not os.path.exists(sdf_file):
                time.sleep(10)
            os.chdir(d3_dir)
            try:
                D3ph4s_processor = FrAncestor3Dph4s()
                D3ph4s_processor.calculate3Dph4s(f"{base_name}_std_scaffold.sdf", cpus)
            finally:
                os.chdir("..")
        D3_run(os.path.basename(molecules_file), cpus)
    except Exception as e:
        log_message(f"🚫 An error occurred in D3: {str(e)}")
    finally:
        os.chdir("..")
    #############################################################
def D3_Processor(libs, cpus):
    files = [lib for lib in libs if os.path.isfile(lib)]
    if not files:
        log_message("No valid molecule files found for 3D processing.")
        return
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = {executor.submit(D3, lib, cpus): lib for lib in files}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log_message(f"🚫 An error occurred in D3_Processor: {str(e)}")
####################################################################################################
#  4 Dimension
####################################################################################################
def D4(molecules_file, cpus):
    folder_name = os.path.splitext(molecules_file)[0]
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    try:
    ##############################################################
        def D4_run(molecules_file, cpus):

            base_name = os.path.splitext(os.path.basename(molecules_file))[0]
            d4_dir = os.path.join(folder_name, "D4")
            os.makedirs(d4_dir, exist_ok=True)
            # Path to expected scaffold output
            scaffold_reps_file = os.path.join(d4_dir, f"{base_name}_std_scaffold_reps.txt")
            if os.path.isfile(scaffold_reps_file):
                log_message(f"ChemicalSpace already prepared... D4 path: {scaffold_reps_file}")
                return 
            D4_processor = FrAncestor4Dligs()
            Threshold = args.Fourth_Dimension if getattr(args, 'Fourth_Dimension', None) is not None else 1.1
            D4_processor.calculate4D(
                molecules_file =f"D3/{molecules_file[:-4]}_std_scaffold_3d_conf.ph4s",
                scaf_file =f"D2/{molecules_file[:-4]}_std_scaffold.txt",
                Threshold=Threshold,
                cpus=cpus)
        D4_run(os.path.basename(molecules_file), cpus)
    except Exception as e:
        log_message(f"🚫 An error occurred in D4: {str(e)}")
    finally:
        os.chdir("..")
    #############################################################
def D4_Processor(libs, cpus):
    files = [lib for lib in libs if os.path.isfile(lib)]
    if not files:
        log_message("No valid molecule files found for 4D processing.")
        return
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = {executor.submit(D4, lib, cpus): lib for lib in files}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log_message(f"🚫 An error occurred in D4_Processor: {str(e)}")
####################################################################################################
#  5 Dimension
####################################################################################################
def D5(molecules_file, cpus, Fifth_Dimension):
    folder_name = os.path.splitext(molecules_file)[0]
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    try:
        def D5_run(molecules_file, cpus, Fifth_Dimension):
            d5_dir = os.path.join(folder_name, "D5")
            os.makedirs(d5_dir, exist_ok=True)
            Fifth_Dimension = os.path.abspath(Fifth_Dimension)
            if not os.path.isfile(Fifth_Dimension):
                raise FileNotFoundError(f"Fifth_Dimension file not found: {Fifth_Dimension}")
            D5_processor = FrAncestor5Dhotspots(Fifth_Dimension=Fifth_Dimension, min_matching_features=args.min_features)
            D5_processor.calculate5D(
                molecules_file=f"D4/{os.path.splitext(molecules_file)[0]}_std_scaffold_3d_conf_reps.ph4s",
                scaffolds=f"D4/{os.path.splitext(molecules_file)[0]}_std_scaffold_reps.txt",
                Fifth_Dimension=Fifth_Dimension,
                Threshold=args.Threshold,
                cpus=cpus)
        D5_run(os.path.basename(molecules_file), cpus, Fifth_Dimension)
    except Exception as e:
        log_message(f"🚫 An error occurred in D5: {str(e)}")
    finally:
        os.chdir("..")
    #############################################################
def D5_Processor(libs, cpus, Fifth_Dimension):
    files = [lib for lib in libs if os.path.isfile(lib)]
    if not files:
        log_message("No valid molecule files found for 5D processing.")
        return
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = {executor.submit(D5, lib, cpus, Fifth_Dimension): lib for lib in files}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log_message(f"🚫 An error occurred in D5_Processor: {str(e)}")
####################################################################################################
#  6 Dimension
####################################################################################################
def D6(molecules_file, cpus, Rx_file, Rx_lib, Rx_Steps):
    folder_name = os.path.splitext(molecules_file)[0]
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    try:
        #############################################################
        def D6_run(molecules_file, cpus, Rx_file, Rx_lib, Rx_Steps):
            d6_dir = os.path.join(folder_name, "D6")
            os.makedirs(d6_dir, exist_ok=True)
            base_name = os.path.splitext(molecules_file)[0]
            os.chdir(d6_dir)
            d4_file = os.path.join("../D4", f"{base_name}_std_scaffold_reps.txt")
            d2_file = os.path.join("../D2", f"{base_name}_std_scaffold.txt")
            if os.path.isfile(d4_file):
                selected_file = d4_file
            else:
                selected_file = d2_file
            D6_processor = FrAncestor6Drxs(Rx_lib=Rx_lib, Rx_Steps=Rx_Steps)
            D6_processor.calculate6D(
                molecules_file=selected_file,
                Rx_file=Rx_file,
                Rx_lib=Rx_lib,
                Rx_Steps=Rx_Steps,
                cpus=cpus)
            os.chdir("..")
        Rx_Steps = int(args.Rx_Steps) if args.Rx_Steps else 1
        D6_run(os.path.basename(molecules_file), cpus, Rx_file, Rx_lib, Rx_Steps)
    except Exception as e:
        log_message(f"🚫 An error occurred in D6: {str(e)}")
    finally:
        os.chdir("..")
    #############################################################
def D6_Processor(libs, cpus, Rx_file, Rx_lib, Rx_Steps):
    files = [lib for lib in libs if os.path.isfile(lib)]
    if not files:
        log_message("No valid molecule files found for 6D processing.")
        return
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = {executor.submit(D6, lib, cpus, Rx_file, Rx_lib, Rx_Steps): lib for lib in files}
        for future in as_completed(futures):
            lib = futures[future]
            try:
                future.result()
            except Exception as e:
                log_message(f"🚫 An error occurred in D6_Processor for {lib}: {str(e)}")
####################################################################################################
#  Chemical Space 
####################################################################################################
def chemicalspace_task(params):
    library_file, chemicalspace_file, sociability_threshold, chunk_size, cpus, sascore, fp = params
    if not isinstance(library_file, str) or not os.path.isfile(library_file):
        raise ValueError(f"Invalid library_file: {library_file}")
    if not isinstance(chemicalspace_file, str) or not os.path.isfile(chemicalspace_file):
        raise ValueError(f"Invalid chemicalspace_file: {chemicalspace_file}")
    lib_base = os.path.splitext(os.path.basename(library_file))[0]
    lib_dir = os.path.dirname(library_file)
    d4_file = os.path.join(lib_dir, "D4", f"{lib_base}_std_scaffold_reps.txt")
    d2_file = os.path.join(lib_dir, lib_base, "D2", f"{lib_base}_std_scaffold.txt")
    selected_file = None
    if os.path.isfile(d4_file):
        selected_file = d4_file
    elif os.path.isfile(d2_file):
        selected_file = d2_file
    else:
        raise FileNotFoundError(f"No scaffold file found for library_file: {selected_file}")

    cs_base = os.path.splitext(os.path.basename(chemicalspace_file))[0]
    cs_dir = os.path.dirname(chemicalspace_file)
    d4_cs_file = os.path.join(cs_dir, cs_base, "D4", f"{cs_base}_std_scaffold_reps.txt")
    d2_cs_file = os.path.join(cs_dir, cs_base, "D2", f"{cs_base}_std_scaffold.txt")
    selected_cs_file = None
    if os.path.isfile(d4_cs_file):
        selected_cs_file = d4_cs_file
    elif os.path.isfile(d2_cs_file):
        selected_cs_file = d2_cs_file
    else:
        raise FileNotFoundError(f"No scaffold file found for chemicalspace_file: {chemicalspace_file}")

    MGT_processor = FrAncestorMGT(fingerprints=fp)
    MGT_processor.calculateMGT(
        selected_file,
        selected_cs_file,
        sociability_threshold,
        chunk_size,
        cpus)

def run_all_combinations(input_files, chemicalspace_files, sociability_threshold, chunk_size, cpus, sascore, fp):
    if isinstance(input_files, str):
        input_files = [input_files]
    if isinstance(chemicalspace_files, str):
        chemicalspace_files = [chemicalspace_files]
    all_tasks = list(itertools.product(input_files, chemicalspace_files))
    per_cpu = max(1, cpus // len(all_tasks))
    params_list = [
        (inp, cs, sociability_threshold, chunk_size, per_cpu, sascore, fp)
        for inp, cs in all_tasks
    ]
    with Pool(min(len(params_list), cpus)) as pool:
        pool.map(chemicalspace_task, params_list)

def run_NearNeighbours_task(
                           library_file, 
                           chemicalspace_file, 
                           sociability_threshold, 
                           chunk_size, 
                           cpus, 
                           sascore,
                           fp=args.fp
                           ):

    base_name = os.path.splitext(os.path.basename(library_file))[0]
    folder_name = base_name
    os.makedirs(folder_name, exist_ok=True)
    os.chdir(folder_name)
    try:
        run_all_combinations(
            [library_file],        
            [chemicalspace_file],
            sociability_threshold,
            chunk_size,
            cpus,
            sascore, 
            fp
        )
    except Exception as e:
        log_message(f"🚫 An error occurred in run_NearNeighbours_task: {str(e)}")
    finally:
        os.chdir("..")
    #############################################################
def NearNeighbours_Processor(molecules_file, chemicalspace_processing, sociability_threshold, chunk_size, cpus, sascore):
    valid_libs = [lib for lib in molecules_file if os.path.isfile(lib)]
    if not valid_libs:
        log_message("No valid molecule files found for NearNeighbours processing.")
        return
    task_pairs = [
        (chem_proc, lib) 
        for chem_proc, lib in itertools.product(chemicalspace_processing, valid_libs)
    ]

    with ProcessPoolExecutor(max_workers=cpus) as executor:
        futures = {
            executor.submit(run_NearNeighbours_task, lib, chem_proc, sociability_threshold, chunk_size, cpus, sascore): (chem_proc, lib)
            for chem_proc, lib in task_pairs
        }
        for future in as_completed(futures):
            chem_proc, lib = futures[future]
            try:
                future.result()
            except Exception as e:
                log_message(f"Error processing {lib} with {chem_proc}: {str(e)}")
###################################################################################################
#  Multiple processes
###################################################################################################
message3 = '''######################################################
# Copyright "©", FrAncestor - Peter E.G.F. Ibrahim.  #
######################################################\n'''
with open("FrAncestor.log", "a") as f:
    f.write(message3)
###################################################################################################
def log_message(message):
    log_file_path = "FrAncestor.log"
    with open(log_file_path, "a") as log_file:
        log_file.write(f"{message}\n")
#################################################
def main_libProcessor(molecules_file, 
                      cpus):
    log_message(f"Starting Chemical Space Processor... {molecules_file}")
    if isinstance(molecules_file, str):
        molecules_file = [molecules_file]
    files = [file for file in molecules_file if os.path.isfile(file)]
    if not files:
        log_message("No valid molecule files found to process.")
        return
    preparelibrary_all(files, cpus)
    log_message(f"✅ ✊ Finished Chemical Space Processor successfully: {molecules_file}")
#################################################
def main_D1_Processor(molecules_file, 
                      cpus, 
                      sascore):
    log_message(f"1D run... {molecules_file}")
    if isinstance(molecules_file, str):
        molecules_file = [molecules_file]
    files = [file for file in molecules_file if os.path.isfile(file)]
    if not files:
        log_message("No valid molecule files found to process.")
        return
    D1_Processor(files, cpus, sascore)
    log_message(f"✅ ☝️ Finished 1D successfully: {molecules_file}")
#################################################
def main_D2_Processor(molecules_file, 
                      cpus):
    log_message(f"2D run... {molecules_file}")
    if isinstance(molecules_file, str):
        molecules_file = [molecules_file]
    files = [file for file in molecules_file if os.path.isfile(file)]
    if not files:
        log_message("No valid molecule files found to process.")
        return
    D2_Processor(files, cpus)
    log_message(f"✅ ✌️ Finished 2D successfully: {molecules_file}")
# #################################################
def main_D3_Processor(molecules_file, 
                      cpus):
    log_message(f"3D run... {molecules_file}")
    if isinstance(molecules_file, str):
        molecules_file = [molecules_file]
    files = [file for file in molecules_file if os.path.isfile(file)]
    if not files:
        log_message("No valid molecule files found to process.")
        return
    D3_Processor(files, cpus)
    log_message(f"✅ 🖖 Finished 3D successfully: {molecules_file}")
#################################################
def main_D4_Processor(molecules_file, 
                      cpus):
    log_message(f"4D run... {molecules_file}")
    if isinstance(molecules_file, str):
        molecules_file = [molecules_file]
    files = [file for file in molecules_file if os.path.isfile(file)]
    if not files:
        log_message("No valid molecule files found to process.")
        return
    D4_Processor(files, cpus)
    log_message(f"✅ 🤟 Finished 4D successfully: {molecules_file}")
#################################################
def main_D5_Processor( molecules_file,
                       cpus,
                       Fifth_Dimension):
    log_message(f"5D run... {molecules_file}")
    if isinstance(molecules_file, str):
        molecules_file = [molecules_file]

    combined_files = (molecules_file or [])
    combined_files = [os.path.abspath(f) for f in combined_files if f]  

    files = sorted(set(combined_files))
    files = [file for file in files if os.path.isfile(file)]
    if isinstance(files, str):
        files = [files]
    files = sorted(set(files))
    files = [file for file in files if os.path.isfile(file)]
    if not files:
        log_message("No valid molecules file found to process.")
        return
    to_process = []
    for chem_file in files:
        base_name = os.path.splitext(os.path.basename(chem_file))[0]
        base_dir = os.path.dirname(chem_file)

        d4_file = os.path.join(base_dir, base_name, "D4", f"{base_name}_std_scaffold_reps.txt")
        d2_file = os.path.join(base_dir, base_name, "D2", f"{base_name}_std_scaffold.txt")

        if os.path.isfile(d4_file):
            log_message(f"Molecules file already prepared... D4 path: {d4_file}")
        elif os.path.isfile(d2_file):
            log_message(f"Molecules file already prepared... D2 path: {d2_file}")
        else:
            to_process.append(chem_file)
    if to_process:
        log_message(f"Molecules file to process: {to_process}")
        preparelibrary_all(to_process, cpus)
        D2_Processor(to_process, cpus)
        if args.Fifth_Dimension and not (args.Third_Dimension or args.Fourth_Dimension):
            D3_Processor(to_process, cpus)
            D4_Processor(to_process, cpus)
            D5_Processor(to_process, cpus, Fifth_Dimension)
            log_message(f"✅ 🖐 Finished 5D successfully: {molecules_file}")
        elif args.Fifth_Dimension and (args.Third_Dimension or args.Fourth_Dimension):
            D5_Processor(to_process, cpus, Fifth_Dimension)
            log_message(f"✅ 🖐 Finished 5D successfully: {molecules_file}")
    else: 
        log_message("✅ All molecules files are already prepared — skipping processing.")
        D3_Processor(molecules_file, cpus)
        D4_Processor(molecules_file, cpus)
        D5_Processor(molecules_file, cpus, Fifth_Dimension)
        log_message(f"✅ 🖐 Finished 5D successfully: {molecules_file}")
        
#################################################
def main_D6_Processor(molecules_file, cpus, Rx_file, Rx_lib, Rx_Steps):
    def ensure_flat_list(inputs):
        if isinstance(inputs, str):
            return [inputs]
        elif isinstance(inputs, list):
            flat = []
            for item in inputs:
                if isinstance(item, list):
                    flat.extend(item)
                else:
                    flat.append(item)
            return flat
        else:
            return []

    def collect_molecule_files(input_sources):
        all_files = set()
        for source in input_sources:
            if not isinstance(source, str):
                continue
            abs_path = os.path.abspath(source)
            if os.path.isfile(abs_path) and abs_path.endswith((".smi", ".sdf")):
                all_files.add(abs_path)
            elif os.path.isdir(abs_path):
                for ext in ("*.smi", "*.sdf"):
                    for file in glob.glob(os.path.join(abs_path, ext)):
                        all_files.add(os.path.abspath(file))
        return sorted(all_files)
    log_message("🔁 6D run started...")
    molecules_file = ensure_flat_list(molecules_file)
    Rx_lib = ensure_flat_list(Rx_lib)
    mol_files = collect_molecule_files(molecules_file)
    rxlib_files = collect_molecule_files(Rx_lib)

    if not mol_files:
        log_message("🚫 No valid molecule files found.")
        return
    if not rxlib_files:
        log_message("🚫 No valid Rx_lib files found.")
        return
    to_process = []
    prepared_rxlib_d0 = []
    for file in rxlib_files:
        base_name = os.path.splitext(os.path.basename(file))[0]
        base_dir = os.path.dirname(file)
        d0_file = os.path.join(base_dir, base_name, "D0", f"{base_name}_std.smi")
        if not os.path.isfile(d0_file):
            to_process.append(file)
        prepared_rxlib_d0.append(d0_file)
    if to_process:
        log_message(f"🛠 Preparing D0 files: {to_process}")
        preparelibrary_all(to_process, cpus)
    else:
        log_message("✅ All D0 Rx_lib files already prepared.")
    valid_rxlibs = [f for f in prepared_rxlib_d0 if isinstance(f, str) and os.path.isfile(f)]
    if not valid_rxlibs:
        log_message("🚫 No valid prepared Rx_lib D0 files found.")
        return
    log_message(f"🚀Reaction libraries : {len(valid_rxlibs)} Rx_lib(s).")
    for rx_lib in valid_rxlibs:
        log_message(f"🚀 Running D6_Processor: {molecules_file} × {os.path.basename(rx_lib)}")
        D6_Processor(molecules_file, cpus, Rx_file, rx_lib, Rx_Steps)
    log_message(f"✅ 🎉 Finished full 6D run on {molecules_file} × {len(valid_rxlibs)} Rx_lib(s).")

#################################################
def main_chemicalspace_processor( 
                       molecules_file,
                       cpus, 
                       chemicalspace_files,
                       sascore):
    if isinstance(molecules_file, str):
        molecules_file = [molecules_file]

    combined_files = (molecules_file or []) + (chemicalspace_files or [])
    combined_files = [os.path.abspath(f) for f in combined_files if f]  

    files = sorted(set(combined_files))
    files = [file for file in files if os.path.isfile(file)]
    if not files:
        log_message("No valid ChemicalSpace files found to process.")
        return
    to_process = []
    for chem_file in files:
        base_name = os.path.splitext(os.path.basename(chem_file))[0]
        base_dir = os.path.dirname(chem_file)
        d4_file = os.path.join(base_dir, base_name, "D4", f"{base_name}_std_scaffold_reps.txt")
        d2_file = os.path.join(base_dir, base_name, "D2", f"{base_name}_std_scaffold.txt")
        if os.path.isfile(d4_file):
            log_message(f"ChemicalSpace already prepared... D4 path: {d4_file}")
        elif os.path.isfile(d2_file):
            log_message(f"ChemicalSpace already prepared... D2 path: {d2_file}")
        else:
            to_process.append(chem_file)
    if to_process:
        log_message(f"ChemicalSpaces to process: {to_process}")
        preparelibrary_all(to_process, cpus)
        if (args.First_Dimension 
            or args.Rule_RO3_RO5 
            or args.MWT 
            or args.nHeavy 
            or args.nHBA 
            or args.nHBD 
            or args.nRB 
            or args.Fsp3 
            or args.nChiC 
            or args.nRings 
            or args.nArRings 
            or args.TPSA
            or args.Halogen 
            or args.logP
            or args.logD 
            or args.QED 
            or args.SFI 
            or args.SAScore):
            D1_Processor(to_process, cpus, sascore)
        D2_Processor(to_process, cpus)
        log_message(f"✅ Finished ChemicalSpace processing successfully: {files}")
    else:
        log_message("✅ All ChemicalSpace files are already prepared — skipping processing.")

#################################################
def main_NearNeighbours(molecules_file, 
                        cpus, 
                        chemicalspace_files,
                        sociability_threshold, 
                        chunk_size, 
                        sascore
                        ):
    if isinstance(chemicalspace_files, str):
        chemicalspace_files = [chemicalspace_files]
    chemicalspace_files = sorted(set(chemicalspace_files))
    chemicalspace_files = [file for file in chemicalspace_files if os.path.isfile(file)]
    if not chemicalspace_files:
        log_message("No valid ChemicalSpace files found to process.")
        return
    if isinstance(molecules_file, str):
        molecules_file = [molecules_file]
    files = [file for file in molecules_file if os.path.isfile(file)]
    NearNeighbours_Processor(files, 
                             chemicalspace_files, 
                             sociability_threshold, 
                             chunk_size, 
                             cpus, 
                             sascore
                            )
    log_message(f"✅ 🙌 All Near-Neighbours have been identified gracefully for  {molecules_file}!")
#################################################
def run_FrAncstor_for_libs(molecules_file,
                           total_files,
                           threads,               # <- single int (per-task threads)
                           sascore,
                           Fifth_Dimension,
                           Rx_file,
                           Rx_lib,
                           Rx_Steps,
                           chemicalspace_files,
                           sociability_threshold,
                           chunk_size,
                           args):
    
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["MKL_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(threads)


    try:
        main_libProcessor(molecules_file, threads)
    except Exception as e:
        print(f"🚫 An error occurred in main_libProcessor: {e}")
    time.sleep(5)
    ######################
    if args.First_Dimension:
        try:
            main_D1_Processor(molecules_file, threads, sascore)
        except Exception as e:
            print(f"🚫 An error occurred in main_D1_Processor: {e}")
        time.sleep(5)
    ######################    
    if any([
        args.Second_Dimension,
        args.Third_Dimension,
        args.Fourth_Dimension,
        args.Fifth_Dimension,
        args.Sixth_Dimension,
        args.MGT
    ]):
        try:
            main_D2_Processor(molecules_file, threads)
        except Exception as e:
            print(f"🚫 An error occurred in main_D2_Processor: {e}")
        time.sleep(5)
    ######################    
    if args.Third_Dimension:
        try:
            main_D3_Processor(molecules_file, threads)
        except Exception as e:
            print(f"🚫 An error occurred in main_D3_Processor: {e}")
        time.sleep(5)
    ######################
    if args.Fourth_Dimension:
        try:
            main_D3_Processor(molecules_file, threads)
        except Exception as e:
            print(f"🚫 An error occurred in main_D3_Processor: {e}")
        time.sleep(5)
        try:
            main_D4_Processor(molecules_file, threads)
        except Exception as e:
            print(f"🚫 An error occurred in main_D4_Processor: {e}")
        time.sleep(5)
    ######################    
    if args.Fifth_Dimension:
        try:
            main_D5_Processor(molecules_file, threads, Fifth_Dimension)
        except Exception as e:
            print(f"🚫 An error occurred in main_D5_Processor: {e}")
        time.sleep(5)
    ######################    
    if args.Sixth_Dimension:
        try:
            main_D6_Processor(molecules_file, threads, Rx_file, Rx_lib, Rx_Steps)
        except Exception as e:
            print(f"🚫 An error occurred in main_D6_Processor: {e}")
        time.sleep(5)
    ######################    
    if args.MGT:
        try:
            total_cpus = args.cpus
            num_files = max(len(chemicalspace_files), 1)
            cpus_per_file = max(1, total_cpus // num_files)
            log_message(f"{cpus_per_file} cpus per Chemical space")
            main_chemicalspace_processor(
                molecules_file=chemicalspace_files,
                cpus=cpus_per_file,
                chemicalspace_files=chemicalspace_files,
                sascore=sascore
            )
        except Exception as e:
            print(f"🚫 An error occurred in main_chemicalspace_processor: {e}")
        try:
            main_NearNeighbours(
                molecules_file, threads,
                chemicalspace_files,
                sociability_threshold,
                chunk_size,
                sascore
            )
        except Exception as e:
            print(f"🚫 An error occurred in main_NearNeighbours: {e}")
        time.sleep(5)
##################################################################################################
def main():
    ####################################
#     EnvironmentGuard().enforce()
    ####################################
    # Step 1: Preload chemicalspace paths
    chemicalspace_files = []
    chemicalspace_set = set()
    chemicalspace_vendors = []
    if args.MGT:
        default_base = '/cluster/ddu/pibrahim001/FrAncestor/chemical_space'
        for entry in args.MGT:
            abs_path = os.path.abspath(entry)
            if entry.lower() == 'chemical_space' and os.path.isdir(default_base):
                for subfolder in os.listdir(default_base):
                    subfolder_path = os.path.join(default_base, subfolder)
                    if os.path.isdir(subfolder_path):
                        chemicalspace_vendors.append(subfolder)
                        for ext in ("*.smi", "*.sdf"):
                            for file in glob.glob(os.path.join(subfolder_path, ext)):
                                abs_file = os.path.abspath(file)
                                if abs_file not in chemicalspace_set:
                                    chemicalspace_files.append(abs_file)
                                    chemicalspace_set.add(abs_file)
            elif os.path.isdir(os.path.join(default_base, entry)):
                subfolder_path = os.path.join(default_base, entry)
                chemicalspace_vendors.append(entry)
                for ext in ("*.smi", "*.sdf"):
                    for file in glob.glob(os.path.join(subfolder_path, ext)):
                        abs_file = os.path.abspath(file)
                        if abs_file not in chemicalspace_set:
                            chemicalspace_files.append(abs_file)
                            chemicalspace_set.add(abs_file)
            elif os.path.isdir(abs_path):
                chemicalspace_vendors.append(os.path.basename(abs_path))
                for ext in ("*.smi", "*.sdf"):
                    for file in glob.glob(os.path.join(abs_path, ext)):
                        abs_file = os.path.abspath(file)
                        if abs_file not in chemicalspace_set:
                            chemicalspace_files.append(abs_file)
                            chemicalspace_set.add(abs_file)
            elif os.path.isfile(abs_path):
                if abs_path.endswith((".smi", ".sdf")) and abs_path not in chemicalspace_set:
                    chemicalspace_files.append(abs_path)
                    chemicalspace_set.add(abs_path)
    molecules_file = []
    def add_if_not_duplicate(path):
        abs_path = os.path.abspath(path)
        if abs_path not in molecules_file:
            molecules_file.append(abs_path)
    if args.SMILES:
        add_if_not_duplicate(args.SMILES)
    if args.SDF:
        add_if_not_duplicate(args.SDF)
    if args.PDB:
        abs_path = os.path.abspath(args.PDB)
        renumber_atoms(abs_path)
        add_if_not_duplicate(abs_path)
    if args.DIRECTORY:
        for ext in ("*.smi", "*.sdf", "*.pdb"):
            for f in glob.glob(os.path.join(args.DIRECTORY, ext)):
                abs_path = os.path.abspath(f)
                if f.endswith(".pdb"):
                    renumber_atoms(abs_path)
                add_if_not_duplicate(abs_path)
    log_message(f"👍 Total input molecule files: {len(molecules_file)}")
    if chemicalspace_vendors:
        vendors_str = ', '.join(set(chemicalspace_vendors))
        log_message(f"🌌 Vendor's Chemical Space ({vendors_str}): {len(chemicalspace_files)} files")
    else:
        log_message("🌌 No chemical space files loaded.")
    molecules_file = list(set(molecules_file))
    valid_chemspace_files = [f for f in chemicalspace_files if os.path.basename(f).lower() != 'done']
    
    # --- CPU planning ---
    requested_cpus = max(1, int(args.cpus))
    requested_cpus = min(requested_cpus, os.cpu_count() or requested_cpus)

    num_files = max(len(chemicalspace_files), len(molecules_file))
    files = chemicalspace_files if len(chemicalspace_files) >= len(molecules_file) else molecules_file
    if num_files == 0:
        raise ValueError("No input files provided")

    # how many jobs in parallel:
    concurrency = min(num_files, requested_cpus)

    # split CPUs across those concurrent jobs (distribute remainder first)
    base = requested_cpus // concurrency
    rem  = requested_cpus %  concurrency
    per_file_cpus = [base + 1] * rem + [base] * (concurrency - rem)   # len == concurrency

    # logging
    log_message(
        f"num_files={num_files}, requested_cpus={requested_cpus}, "
        f"concurrency={concurrency}, base={base}, rem={rem}, "
        f"counts={{2:{per_file_cpus.count(2)}, 1:{per_file_cpus.count(1)}}}"
    )
    print("Per-file CPU allocations:", per_file_cpus)
    
    def get_absolute_path(file_path):
        return os.path.abspath(file_path) if file_path else None
    Fifth_Dimension = get_absolute_path(args.Fifth_Dimension) if args.Fifth_Dimension else None
    FrAncestor_dir = os.path.dirname(os.path.abspath(__file__))
    sascore = os.path.join(FrAncestor_dir, "modules", ".D1_fpscores.pkl.gz") 
    Rx_file = os.path.join(FrAncestor_dir, "modules", ".D6_Reactions_corrected.csv") 
    Rx_file = get_absolute_path(Rx_file) if args.Sixth_Dimension else None
    Rx_lib = get_absolute_path(args.Sixth_Dimension) if args.Sixth_Dimension else None
    Rx_Steps = int(args.Rx_Steps) if args.Rx_Steps is not None else (1 if args.Sixth_Dimension else None)

    def get_absolute_path(path):
        return os.path.abspath(path)
    def collect_chemicalspace_files(cs_args):
        default_base = '/cluster/ddu/pibrahim001/FrAncestor/chemical_space'
        collected_files = []
        if not cs_args:
            return collected_files
        for arg in cs_args:
            if arg.lower() == 'chemical_space' and os.path.isdir(default_base):
                for subdir in os.listdir(default_base):
                    subdir_path = os.path.join(default_base, subdir)
                    if os.path.isdir(subdir_path):
                        for f in os.listdir(subdir_path):
                            if f.endswith('.smi') or f.endswith('.sdf'):
                                collected_files.append(get_absolute_path(os.path.join(subdir_path, f)))
            elif os.path.isdir(os.path.join(default_base, arg)):
                subdir_path = os.path.join(default_base, arg)
                for f in os.listdir(subdir_path):
                    if f.endswith('.smi') or f.endswith('.sdf'):
                        collected_files.append(get_absolute_path(os.path.join(subdir_path, f)))
            elif os.path.exists(arg):
                if os.path.isdir(arg):
                    for f in os.listdir(arg):
                        if f.endswith('.smi') or f.endswith('.sdf'):
                            collected_files.append(get_absolute_path(os.path.join(arg, f)))
                elif os.path.isfile(arg) and (arg.endswith('.smi') or arg.endswith('.sdf')):
                    collected_files.append(get_absolute_path(arg))
            else:
                for path in glob.glob(arg):
                    if os.path.isfile(path) and (path.endswith('.smi') or path.endswith('.sdf')):
                        collected_files.append(get_absolute_path(path))
        return collected_files
    chemicalspace_files = collect_chemicalspace_files(args.MGT) or []
    sociability_threshold = float(args.sociability_threshold) if args.sociability_threshold is not None else 0.0
    chunk_size = int(args.chunk_size) if args.chunk_size is not None else 1000
    ###########################################
    with ProcessPoolExecutor(max_workers=concurrency) as executor:
        cpu_cycle = cycle(per_file_cpus)  # assign threads per job
        futures = []
        for mol_file in molecules_file:
            threads = next(cpu_cycle)
            futures.append(
                executor.submit(
                    run_FrAncstor_for_libs,
                    mol_file,
                    len(molecules_file),
                    threads,                 # <- pass single int
                    sascore,
                    Fifth_Dimension,
                    Rx_file,
                    Rx_lib,
                    Rx_Steps,
                    chemicalspace_files,
                    sociability_threshold,
                    chunk_size,
                    args,
                )
            )
        for future in as_completed(futures):
            try:
                _ = future.result()
            except Exception as e:
                import traceback
                log_message("🚫 Full traceback:")
                log_message(traceback.format_exc())

    log_message("✅ FrAncestor's job is Completed🫡! - You are awesome 🫵😎, keep up the great work!🙏")
##################################################################################################
if __name__ == '__main__':
    main()
