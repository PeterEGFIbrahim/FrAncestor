import os
from rdkit import Chem
from rdkit.Chem import QED
from rdkit.Chem import rdMolDescriptors
import rdkit.Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Crippen
import os.path as op
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import multiprocessing
from multiprocessing import Pool, cpu_count
import math
import gzip
import pickle 
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pathlib import Path
import argparse
from .FrAncestor_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

class FrAncestor1Descs:
    def __init__(self, sascore_file=None):
        self._fscores = None
        if sascore_file:
            self.readFragmentScores(sascore_file)

    def readFragmentScores(self, sascore_file):
        try:
            with gzip.open(sascore_file, 'rb') as file:
                raw_scores = pickle.load(file)
                outDict = {entry[j]: float(entry[0]) for entry in raw_scores for j in range(1, len(entry))}
                self._fscores = outDict
        except FileNotFoundError:
            raise ValueError(f"File {sascore_file} not found.")
        except Exception as e:
            raise ValueError(f"Error reading fragment scores: {e}")

    def numBridgeheadsAndSpiro(self, mol):
        nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        return nBridgehead, nSpiro

    def calculateScore(self, m):
        if self._fscores is None:
            self.readFragmentScores()
        fp = rdMolDescriptors.GetMorganFingerprint(m, 2)
        fps = fp.GetNonzeroElements()
        score1 = 0.0
        nf = 0
        for bitId, v in fps.items():
            nf += v
            sfp = bitId
            score1 += self._fscores.get(sfp, -4) * v
        score1 /= nf if nf > 0 else 1 
        nAtoms = m.GetNumAtoms()
        nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
        nBridgeheads, nSpiro = self.numBridgeheadsAndSpiro(m)
        ri = m.GetRingInfo()
        nMacrocycles = sum(1 for ring in ri.AtomRings() if len(ring) > 8)
        sizePenalty = nAtoms**1.005 - nAtoms
        stereoPenalty = math.log10(nChiralCenters + 1)
        spiroPenalty = math.log10(nSpiro + 1)
        bridgePenalty = math.log10(nBridgeheads + 1)
        macrocyclePenalty = math.log10(2) if nMacrocycles > 0 else 0.0
        score2 = -sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty
        score3 = 0.0
        if nAtoms > len(fps):
            score3 = math.log(float(nAtoms) / len(fps)) * 0.5
        sascore = score1 + score2 + score3
        min_score, max_score = -4.0, 2.5
        sascore = 11.0 - (sascore - min_score + 1.0) / (max_score - min_score) * 9.0
        if sascore > 8.0:
            sascore = 8.0 + math.log(sascore + 1.0 - 9.0)
        if sascore > 10.0:
            sascore = 10.0
        elif sascore < 1.0:
            sascore = 1.0
        return sascore

    def calculate1Descs(self, molecules_file, output_file, num_cpus, sascore):
        self.descriptors(molecules_file, output_file, num_cpus)
        # self.plot_descriptors(output_file)
        self.plot_descriptors_better(output_file)

    def process_chunk(self, smiles_chunk):
        return [self.process_molecule(smiles_data) for smiles_data in smiles_chunk]

    def molecule_rings(self, mol):
        try:
            ring_info = mol.GetRingInfo()
            if not ring_info or not ring_info.AtomRings():
                return []
            rings = [set(ring_atoms) for ring_atoms in ring_info.AtomRings()]
            ring_systems = {}
            for idx, ring in enumerate(rings):
                if idx not in ring_systems:
                    ring_systems[idx] = idx
                for other_idx, other_ring in enumerate(rings[idx + 1:], start=idx + 1):
                    if other_idx not in ring_systems:
                        common_atoms = ring & other_ring
                        if common_atoms:
                            ring_systems[other_idx] = ring_systems[idx]
            grouped_rings = {}
            for idx, ring in enumerate(rings):
                ring_idx = ring_systems[idx]
                grouped_rings.setdefault(ring_idx, set()).update(ring)
            return list(grouped_rings.values())
        except Exception as e:
            print(f"Error in molecule_rings: {e}")
            return [] 
       
    def molecular_features(self, mol):
        feature_patterns = {"Halogen": "[F,Cl,Br,I]"}
        counts = {name: len(mol.GetSubstructMatches(Chem.MolFromSmarts(smart)))
                  for name, smart in feature_patterns.items()}
        return counts["Halogen"]

    def calculate_molecular_descriptors(self, mol):
        if mol is None:
            raise ValueError("Invalid molecule: `mol` is None or cannot be processed.")
        MWT = rdMolDescriptors.CalcExactMolWt(mol)
        nHeavy = mol.GetNumHeavyAtoms()
        nHBA = rdMolDescriptors.CalcNumHBA(mol)
        nHBD = rdMolDescriptors.CalcNumHBD(mol)
        nRB = rdMolDescriptors.CalcNumRotatableBonds(mol)
        Fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
        nChiC = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        nRg = rdMolDescriptors.CalcNumRings(mol)
        nAr = rdMolDescriptors.CalcNumAromaticRings(mol)
        TPSA = rdMolDescriptors.CalcTPSA(mol)
        Halogen = self.molecular_features(mol)
        logP = Crippen.MolLogP(mol)
        logd_ph = 0.74
        logD = logP - logd_ph
        qed = QED.qed(mol)
        SFI = logD + nAr 
        SAScore = self.calculateScore(mol)
        return [MWT, nHeavy, nHBA, nHBD, nRB,  Fsp3, nChiC, 
                nRg, nAr, 
                TPSA, Halogen, logP, logD, 
                qed, SFI, SAScore
                ]

    def process_molecule(self, smiles_data):
        smiles, mol_id, inchikey = smiles_data
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            descriptors = self.calculate_molecular_descriptors(mol)
            formatted_descriptors = ["{:.2f}".format(value) if isinstance(value, (float, int)) else value for value in descriptors]
            return [smiles, mol_id, inchikey] + formatted_descriptors
        return None

    def plot_descriptors(self, descriptors_file):
        data = pd.read_csv(descriptors_file, sep='\t')
        descriptors = [
            'MWT', 'nHeavy', 'nHBA', 'nHBD', 'nRB', 'Fsp3', 'nChiC',
            'nRings', 'nArRings', 
            'TPSA', 'Halogen', 'logP',
            'logD', 'QED', 'SFI', 'SAScore'
        ]
        colors = sns.color_palette('coolwarm', n_colors=len(descriptors))
        fig, axes = plt.subplots(nrows=6, ncols=3, figsize=(18, 12))
        axes = axes.flatten()
        for idx, (descriptor, color) in enumerate(zip(descriptors, colors)):
            sns.histplot(data[descriptor], bins=20, kde=True, color=color, ax=axes[idx])
            axes[idx].set_xlabel(descriptor, fontsize=12, fontweight='bold')
            axes[idx].set_ylabel('Count', fontsize=10)
            axes[idx].set_title(f'{descriptor} Distribution', fontsize=14, fontweight='bold')
            axes[idx].grid(True, linestyle='--', alpha=0.7)
        for ax in axes[len(descriptors):]:
            fig.delaxes(ax)
        plt.tight_layout(pad=2.0)
        save_path = f'{descriptors_file[:-4]}_enhanced.png'
        plt.savefig(save_path, dpi=300)

    def plot_descriptors_better(self, descriptors_file):
        combined_df = pd.read_csv(descriptors_file, sep="\t")
        column_rename_map = {
            "Molecular_Weight": "MWT",
            "Num_AromaticRings": "Aromatic_Rings",
            "Num_H_Acceptors": "HBA",
            "Num_H_Donors": "HBD",
            "Molecular_PolarSurfaceArea": "TPSA",
        }
        combined_df.rename(columns=column_rename_map, inplace=True)
        columns_to_cast = [
            "MWT", "nHeavy", "nHBA", "nHBD", "nRB", "nChiC", "nRings", "nArRings", "Halogen"
        ]
        for col in columns_to_cast:
            if col in combined_df.columns:
                combined_df[col] = combined_df[col].fillna(0).astype(int)
        numeric_df = combined_df.iloc[:, 3:].select_dtypes(include=["number"])
        melted_df = numeric_df.melt(var_name="key", value_name="value")
        unique_keys = melted_df['key'].unique()
        colors = sns.color_palette("husl", len(unique_keys)) 
        color_map = dict(zip(unique_keys, colors))
        g = sns.FacetGrid(melted_df, col="key", col_wrap=4, sharex=False, sharey=False)
        
        for ax, key in zip(g.axes.flat, unique_keys):
            subset = melted_df[melted_df['key'] == key]
            sns.histplot(
                data=subset, 
                x="value", 
                ax=ax, 
                color=color_map[key], 
                kde=False, 
                alpha=0.7
            )
            median_value = subset["value"].median()
            ax.axvline(median_value, color='black', linestyle='--', linewidth=1.5)
            ax.text(
                median_value, 
                ax.get_ylim()[1] * 0.9, 
                f"{median_value:.2f}", 
                color='black', 
                fontsize=9, 
                ha='center'
            )
            if key == "MWT": 
                max_value = subset["value"].max()
                tick_spacing = max_value // 5
                ax.xaxis.set_major_locator(plt.MultipleLocator(tick_spacing))
            elif key in columns_to_cast:
                ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
            ax.set_title(key, fontweight='bold')  
        
        g.set_axis_labels("Value", "Count")
        g.fig.subplots_adjust(top=0.9)
        filename = Path(descriptors_file).stem
        g.fig.suptitle(f"{filename}", fontsize=16, fontweight='bold')
        output_file = f"1st-D_FrAncestor_Molecular_Descriptors_for_{filename}.png"
        plt.savefig(output_file)

    def descriptors(self, smiles_file, output_file, num_cpus):
        with open(smiles_file, 'r') as f:
            lines = [line.strip().split('\t') for line in f.readlines()]
            smiles_data = [(line[0], line[1], line[2]) for line in lines]
        chunk_size = len(smiles_data)
        smiles_chunks = [smiles_data[i:i + chunk_size] for i in range(0, len(smiles_data), chunk_size)]
        with Pool(processes=num_cpus) as pool:
            results = pool.map(self.process_chunk, smiles_chunks)
        with open(output_file, 'w') as output:
            header = "smiles\tMOL_ID\tINCHI_Key\tMWT\tnHeavy\tnHBA\tnHBD\tnRB\tFsp3\tnChiC\tnRings\tnArRings\tTPSA\tHalogen\tlogP\tlogD\tQED\tSFI\tSAScore\n"
            output.write(header)
            for chunk_results in results:
                for res in chunk_results:
                    if res:
                        output.write("\t".join(map(str, res)) + "\n")

if __name__ == "__main__":
    """
    FrAncestor V.0.1 - 1D Descriptors - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="D1: Descriptor calculation from SMILES input")
    parser.add_argument('-i', '--molecules_file', required=True, help='Input SMILES or molecule file')
    parser.add_argument('-o', '--output_file', required=True, help='Output file for descriptors')
    parser.add_argument('-c', '--cpus', type=int, default=1, help='Number of CPUs to use')
    parser.add_argument('-s', '--sascore', required=True, help='Path to SAScore file')
    args = parser.parse_args()
    D1_processor = FrAncestor1Descs(sascore_file=args.sascore)
    D1_processor.calculate1Descs(args.molecules_file, args.output_file, args.cpus, args.sascore)
