import os
from multiprocessing import Pool, cpu_count
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from .FrAncestor_utility import EnvironmentGuard
# EnvironmentGuard().enforce()

class DescFilters:
    def __init__(self, filter_criteria=None):
        self.filter_type_thresholds = {}
        if filter_criteria:
            self.set_filter_criteria(filter_criteria)

    def set_filter_criteria(self, filter_criteria):
        for i in range(0, len(filter_criteria), 2):
            filter_type = filter_criteria[i]
            thresholds = tuple(map(float, filter_criteria[i + 1].split(','))) if i + 1 < len(filter_criteria) else None
            self.filter_type_thresholds[filter_type] = thresholds

    def processD1filter(self, smi_file, output_file, num_cpus):
        with open(smi_file, 'r') as file:
            lines = file.readlines()
        if not lines:
            print(f"Input file {smi_file} is empty. Exiting...")
            return
        header = lines[0]
        data = lines[1:]
        chunk_size = max(len(data) // (num_cpus * 10), 1)
        data_chunks = list(self.generate_chunks(data, chunk_size))
        filtered_lines = []
        with Pool(processes=num_cpus) as pool:
            results = pool.map(self.filter_data_chunk, data_chunks)
        for chunk_results in results:
            filtered_lines.extend(chunk_results)
        self.write_filtered_data(output_file, header, filtered_lines)
        # self.plot_descriptors(output_file)
        self.plot_descriptors_better(output_file)

    def filter_data_chunk(self, chunk):
        filtered_lines = []
        for line in chunk:
            columns = line.strip().split('\t')
            if len(columns) < 15:
                print(f"Skipping malformed row: {line.strip()}")
                continue
            pass_filter = True
            for filter_type, threshold in self.filter_type_thresholds.items():
                try:
                    value = self.extract_descriptor_value(columns, filter_type)
                    if filter_type == 'logP':
                        if threshold == (-3, 3) and not (-3 <= value <= 3):
                            pass_filter = False
                            break
                        elif threshold == (-5, 5) and not (-5 <= value <= 5):
                            pass_filter = False
                            break
                    else:
                        if threshold and not (threshold[0] <= value <= threshold[1]):
                            pass_filter = False
                            break
                except (ValueError, IndexError) as e:
                    print(f"Error processing row: {line.strip()} ({e})")
                    pass_filter = False
                    break
            if pass_filter:
                filtered_lines.append(line)
        return filtered_lines

    @staticmethod
    def extract_descriptor_value(columns, filter_type):
        descriptor_index_map = {
            'MWT': 3, 'nHeavy': 4, 'nHBA': 5, 'nHBD': 6, 'nRB': 7, 'Fsp3': 8,
            'nChiC': 9, 'nRings': 10, 'nArRings': 11, 'TPSA': 12, 'Halogen': 13, 'logP': 14,
            'logD': 15, 'QED': 16, 'SFI': 17, 'SAScore': 18
        }
        return float(columns[descriptor_index_map[filter_type]])

    def generate_chunks(self, data, chunk_size):
        for i in range(0, len(data), chunk_size):
            yield data[i:i + chunk_size]

    @staticmethod
    def write_filtered_data(output_filename, header, filtered_lines):
        with open(output_filename, 'w') as output_file:
            output_file.write(header)
            output_file.writelines(filtered_lines)

    @staticmethod
    def plot_descriptors(output_filename):
        data = pd.read_csv(output_filename, sep='\t')
        descriptors = [
            'MWT', 'nHeavy', 'nHBA', 'nHBD', 'nRB', 'Fsp3',
            'nChiC', 'nRings', 'nArRings', 'TPSA', 'Halogen', 'logP',
            'logD', 'QED', 'SFI', 'SAScore'
        ]
        colors = sns.color_palette('coolwarm', n_colors=len(descriptors))
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))
        axes = axes.flatten()
        for idx, (descriptor, color) in enumerate(zip(descriptors, colors)):
            sns.histplot(data[descriptor], bins=20, kde=True, color=color, ax=axes[idx])
            axes[idx].set_xlabel(descriptor, fontsize=12, fontweight='bold')
            axes[idx].set_ylabel('Count', fontsize=10)
            axes[idx].set_title(f'{descriptor}', fontsize=14, fontweight='bold')
        for ax in axes[len(descriptors):]:
            fig.delaxes(ax)
        plt.tight_layout(pad=2.0)
        plt.savefig(f'{output_filename[:-4]}_filtered.png', dpi=300)

    @staticmethod
    def plot_descriptors_better(output_filename):
        combined_df = pd.read_csv(output_filename, sep="\t")
        column_rename_map = {
            "Molecular_Weight": "MWT",
            "Num_AromaticRings": "Aromatic_Rings",
            "Num_H_Acceptors": "HBA",
            "Num_H_Donors": "HBD",
            "Molecular_PolarSurfaceArea": "TPSA",
        }
        combined_df.rename(columns=column_rename_map, inplace=True)
        columns_to_cast = [
            'MWT', 'nHeavy', 'nHBA', 'nHBD', 'nRB', 'Fsp3',
            'nChiC', 'nRings', 'nArRings', 'TPSA', 'Halogen', 'logP',
            'logD', 'QED', 'SFI', 'SAScore'
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
            # Add a dotted line for the median
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
        filename = Path(output_filename).stem
        g.fig.suptitle(f"{filename}", fontsize=16, fontweight='bold')
        output_file = f"1st-D_FrAncestor_Molecular_Descriptors_for_{filename}_filtered.png"
        plt.savefig(output_file)

if __name__ == "__main__":
    """
    FrAncestor V.0.1 - 1D Descriptors filter - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="D1: Descriptor calculation from SMILES input")
    parser.add_argument('-i', '--desc_file', required=True, help='Input SMILES or molecule file')
    parser.add_argument('-o', '--output_file', required=True, help='Output file for descriptors')
    parser.add_argument('-c', '--cpus', type=int, default=1, help='Number of CPUs to use')
    parser.add_argument('-1D', '--First_Dimension', required=True, help='Path to SAScore file')
    args = parser.parse_args()
    D1_filter_processor = DescFilters(filter_criteria=args.First_Dimension)
    D1_filter_processor.processD1filter(args.desc_file, args.output_file, args.cpus)
    