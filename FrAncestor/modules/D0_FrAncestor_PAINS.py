from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from multiprocessing import Pool, cpu_count
import shutil
import argparse
from .FrAncestor_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

class PAINSFilter:
    def __init__(self, filter_type="PAINS"):
        self.catalog = self._init_catalog(filter_type)

    def _init_catalog(self, filter_type):
        if filter_type == "PAINS":
            params = FilterCatalogParams()
            params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
            return FilterCatalog(params)
        else:
            raise ValueError(f"Unknown filter type: {filter_type}")

    def _passes_reos(self, mol):
        if not mol:
            return False
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Chem.Lipinski.NumHDonors(mol)
        hba = Chem.Lipinski.NumHAcceptors(mol)
        return mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10

    def _filter_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            if self.catalog.HasMatch(mol) or not self._passes_reos(mol):
                return (smiles, False)  
            else:
                return (smiles, True) 
        return None

    def _process_smiles(self, smiles_batch):
        filtered_batch = []
        removed_batch = []
        for smiles in smiles_batch:
            result = self._filter_smiles(smiles)
            if result:
                if result[1]: 
                    filtered_batch.append(result[0])  
                else:
                    removed_batch.append(result[0])  
        return filtered_batch, removed_batch

    def filter_smiles_file(self, smiles_file, output_removed_file, num_cpus=None):
        if num_cpus is None:
            num_cpus = cpu_count()
        with open(smiles_file, 'r') as input_file:
            smiles_list = [line.strip() for line in input_file]
        chunk_size = max(1, len(smiles_list) // num_cpus)
        smiles_chunks = [smiles_list[i:i + chunk_size] for i in range(0, len(smiles_list), chunk_size)]
        with Pool(processes=num_cpus) as pool:
            results = pool.map(self._process_smiles, smiles_chunks)
        filtered_smiles = []
        removed_smiles = []
        for filtered_batch, removed_batch in results:
            filtered_smiles.extend(filtered_batch)
            removed_smiles.extend(removed_batch)
        temp_filtered_file = "temp_filtered.smi"
        with open(temp_filtered_file, 'w') as temp_filtered, \
             open(output_removed_file, 'w') as output_removed:
            for smiles in filtered_smiles:
                temp_filtered.write(smiles + '\n')
            for smiles in removed_smiles:
                output_removed.write(smiles + '\n')
        shutil.move(temp_filtered_file, smiles_file)

if __name__ == "__main__":
    """
    FrAncestor V 0.1 - PAINS Filter - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="PAINS filtering of SMILES files")
    parser.add_argument('-i', '--input_file', required=True, help='Input SMILES file')
    parser.add_argument('-o', '--output_file', required=True, help='Output file for removed SMILES')
    parser.add_argument('-c', '--cpus', type=int, default=None, help='Number of CPUs to use')
    args = parser.parse_args()
    standardizer = PAINSFilter()
    standardizer.filter_smiles_file(args.input_file, args.output_file, args.cpus)
    