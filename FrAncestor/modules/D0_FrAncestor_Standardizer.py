import os
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from standardiser import break_bonds, neutralise, rules
from multiprocessing import Pool, cpu_count
import argparse
from .FrAncestor_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

class MoleculeStandardizer:
    def __init__(self):
        RDLogger.DisableLog('rdApp.*')

    def smiles_file(self, smiles_molecules):
        allowed_elements = {"C", "N", "O", "S", "P", "H", "F", "Cl", "Br", "I", "B"}
        molecule_counter = 1
        with open(smiles_molecules) as smi:
            for line in smi:
                parts = line.strip().split()
                if len(parts) < 1:
                    continue
                smiles_structure = parts[0]
                mol_ID = parts[1] if len(parts) > 1 else f"smi_{molecule_counter}"
                molecule_counter += 1
                try:
                    smiles_molecule = Chem.MolFromSmiles(smiles_structure, sanitize=False)
                    if smiles_molecule:
                        if all(atom.GetSymbol() in allowed_elements for atom in smiles_molecule.GetAtoms()):
                            smiles_molecule.SetProp('_ID', mol_ID)
                            yield smiles_molecule, mol_ID
                        else:
                            print(f"Unsupported elements detected in SMILES: {smiles_structure}. Skipping line.")
                    else:
                        print(f"Invalid SMILES structure in line: {line.strip()}. Skipping this line.")
                except Exception as e:
                    print(f"Error processing SMILES '{smiles_structure}': {e}. Skipping this line.")

    def sanitizer(self, molecule):
        periodic_table = Chem.GetPeriodicTable()
        for atom in molecule.GetAtoms():
            symbol = atom.GetSymbol()
            valence_bond_values = [3, 5] if symbol == 'P' else periodic_table.GetValenceList(symbol)
            available_states = [atom.GetExplicitValence() - valence for valence in valence_bond_values]
            if 0 in available_states:
                atom.SetFormalCharge(0)
            else:
                max_state = max(available_states)
                if max_state < 0:
                    atom.SetFormalCharge(0)
                    atom.SetNumExplicitHs(-max_state)
                else:
                    atom.SetFormalCharge(min(available_states))
        Chem.SanitizeMol(molecule)
        return Chem.RemoveHs(molecule)

    def standardize_molecule(self, molecule):
        try:
            molecule = Chem.RemoveHs(molecule, sanitize=True)
            molecule = neutralise.run(molecule)
            molecule = self.sanitizer(molecule)
            molecule = break_bonds.run(molecule)
            molecule = rules.run(molecule, verbose=True)
            
            # Keep only the largest fragment
            fragments = Chem.GetMolFrags(molecule, asMols=True, sanitizeFrags=True)
            if fragments:
                molecule = max(fragments, key=lambda frag: frag.GetNumAtoms())

            Chem.SanitizeMol(molecule)
        except Exception:
            return None, "Standardization failed"
        return molecule, "Standardization successful"

    def generate_inchikey(self, molecule):
        try:
            generated_inchi = inchi.MolToInchi(molecule)
            if not generated_inchi:
                return False, "Failed to generate InChI"
            inchi_key = inchi.InchiToInchiKey(generated_inchi)
            if not inchi_key:
                return False, "Failed to generate InChIKey"
            return inchi_key, "InChI and InChIKey generated successfully"
        except Exception as e:
            return False, f"Error during InChI generation: {str(e)}"

    def process_molecule(self, molecule, molecule_id, index):
        try:
            molecule_name = molecule_id.strip()
        except Exception as e:
            print(f"Error retrieving molecule ID for molecule {index}: {str(e)}")
            molecule_name = f"unnamed_{index:09d}"
        standardized_molecule, std_success = self.standardize_molecule(molecule)
        if std_success == "Standardization successful":
            inchikey, inchikey_msg = self.generate_inchikey(standardized_molecule)
            if inchikey:
                return Chem.MolToSmiles(standardized_molecule, isomericSmiles=True), molecule_name, inchikey
        return None

    def process_molecule_chunk(self, molecules_chunk, chunk_id):
        results = []
        for mol, mol_id in molecules_chunk:
            result = self.process_molecule(mol, mol_id, chunk_id)
            results.append(result)
        return results

    def standardize_file(self, molecules_filepath, output, num_cpus=None, batch_size=10000):
        if num_cpus is None:
            num_cpus = cpu_count()
        file_name, file_extension = os.path.splitext(molecules_filepath)
        file_extension = file_extension.lower()[1:]
        if file_extension not in {"smi", "txt", "csv", "sdf", "pdb"}:
            raise ValueError(f"File format '{file_extension}' is not supported!")
        if file_extension == "pdb":
            sdf_path = f"{file_name}.sdf"
            print(f"Converting {molecules_filepath} to {sdf_path}...")
            os.system(f'obabel -i pdb {molecules_filepath} -o sdf -O {sdf_path}')
            molecules_filepath = sdf_path
            file_extension = "sdf"

        def sdf_batch_generator(path, batch_size):
            suppl = Chem.SDMolSupplier(path, sanitize=False)
            batch = []
            for i, mol in enumerate(suppl):
                if mol is not None:
                    # Try common ID fields in priority order
                    if mol.HasProp("IDNUMBER"):
                        mol_id = mol.GetProp("IDNUMBER")
                    elif mol.HasProp("ID"):
                        mol_id = mol.GetProp("ID")
                    elif mol.HasProp("_Name") and mol.GetProp("_Name").strip():
                        mol_id = mol.GetProp("_Name")
                    elif mol.HasProp("Catalog ID") and mol.GetProp("Catalog ID").strip():
                        mol_id = mol.GetProp("Catalog ID")
                    else:
                        mol_id = f"sdf_{i+1}"  
                    batch.append((mol, mol_id))
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
            if batch:
                yield batch

        def smiles_batch_generator(path, batch_size):
            batch = []
            molecule_counter = 1
            with open(path) as smi:
                for line in smi:
                    parts = line.strip().split()
                    if len(parts) < 1:
                        continue
                    smiles_structure = parts[0]
                    mol_ID = parts[1] if len(parts) > 1 else f"smi_{molecule_counter}"
                    molecule_counter += 1
                    try:
                        mol = Chem.MolFromSmiles(smiles_structure, sanitize=False)
                        if mol:
                            mol.SetProp('_ID', mol_ID)
                            batch.append((mol, mol_ID))
                    except:
                        continue
                    if len(batch) >= batch_size:
                        yield batch
                        batch = []
            if batch:
                yield batch
        batch_generator = sdf_batch_generator if file_extension == "sdf" else smiles_batch_generator
        with open(output, "w") as output_file:
            seen_inchikeys = set()
            with Pool(processes=num_cpus) as pool:
                for batch_index, mega_batch in enumerate(batch_generator(molecules_filepath, batch_size * num_cpus)):
                    chunks = [mega_batch[i:i + batch_size] for i in range(0, len(mega_batch), batch_size)]
                    chunk_results = pool.starmap(self.process_molecule_chunk, [(chunk, idx) for idx, chunk in enumerate(chunks)])
                    for chunk in chunk_results:
                        for result in chunk:
                            if result:
                                smiles, mol_ID, inchikey = result
                                if inchikey not in seen_inchikeys:
                                    seen_inchikeys.add(inchikey)
                                    output_file.write(f"{smiles}\t{mol_ID}\t{inchikey}\n")

if __name__ == "__main__":
    """
    FrAncestor V.0.1 - MoleculeStandardizer - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="PAINS filtering of SMILES files")
    parser.add_argument('-i', '--input_file', required=True, help='Input SMILES file')
    parser.add_argument('-o', '--output_file', required=True, help='Output file for standardised')
    parser.add_argument('-c', '--cpus', type=int, default=None, help='Number of CPUs to use')
    args = parser.parse_args()
    standardizer = MoleculeStandardizer()
    standardizer.standardize_file(args.input_file, args.output_file, args.cpus)
