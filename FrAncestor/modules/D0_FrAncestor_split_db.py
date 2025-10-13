import argparse
import os
import shutil
# EnvironmentGuard().enforce()

class DBSpliter:
    def __init__(self):
        RDLogger.DisableLog('rdApp.*')

    def DBSpliter_file(file, splits):
        if not os.path.exists(file):
            return
        self.recursive_split(file, splits)

    def split_file(input_file, lines_per_file, output_prefix):
        with open(input_file, 'r') as file:
            file_number = 1
            output_file = open(f"{output_prefix}_{file_number}.smi", 'w')
            for i, line in enumerate(file, start=1):
                if i % lines_per_file == 0:
                    output_file.close()
                    file_number += 1
                    output_file = open(f"{output_prefix}_{file_number}.smi", 'w')
                output_file.write(line)
            output_file.close()

    def recursive_split(input_file, splits, level=0):
        if level >= len(splits):
            return 
        num_files = splits[level]
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        with open(input_file, 'r') as file:
            total_lines = sum(1 for _ in file)
        lines_per_file = max(1, total_lines // num_files)
        self.split_file(input_file, lines_per_file, base_name)
        if level + 1 < len(splits):
            for i in range(1, num_files + 1):
                split_file_name = f"{base_name}_{i}.smi"
                if os.path.exists(split_file_name):
                    folder_name = os.path.splitext(split_file_name)[0]
                    os.makedirs(folder_name, exist_ok=True)
                    new_file_path = os.path.join(folder_name, split_file_name)
                    shutil.move(split_file_name, new_file_path)
                    os.chdir(folder_name)
                    self.recursive_split(split_file_name, splits, level + 1)
                    if os.path.exists(split_file_name):
                        os.remove(split_file_name)
                    os.chdir('..')
        else:
            print(f"Final split at level {level + 1} completed for {input_file}.")

if __name__ == "__main__":
    """
    FrAncestor V.0.1 - data-spliter - Copyright "©" 2024, Peter E.G.F. Ibrahim.
    """
    parser = argparse.ArgumentParser(description="Split a file into multiple smaller files and organize them.")
    parser.add_argument('-f', '--file', help='File name to be split into smaller files', required=True)
    parser.add_argument('-n', '--splits', type=int, nargs='+', help='Sequence of splits at each level', required=True)
    args = parser.parse_args()
    Spliter = DBSpliter()
    Spliter.DBSpliter_file(args.file, args.splits)
