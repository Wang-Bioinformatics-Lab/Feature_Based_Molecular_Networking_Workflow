import sys
import argparse
import pandas as pd
from .pfa_parser import (
    ReadPFAFile,
    GetSampleGroupingFromPFA,
    GetCompoundResultsFromPFA,
    GetSCompositeMSMSFromPFA,
    ExportCompoundResults
)
import logging
from pathlib import Path
from typing import List, Dict


def write_to_mgf(output_path: Path, ms2_scans: List[Dict]):
    with open(output_path, 'w', encoding='utf-8') as f:
        for scan in ms2_scans:
            _ion_mode = str(scan['ion_mode']).upper().strip()
            if _ion_mode == 'N':
                _ion_mode = '-'
            elif _ion_mode == 'P':
                _ion_mode = '+'
            else:
                raise ValueError(f"Unknown ion mode: {_ion_mode}")
            f.write("BEGIN IONS\n")
            f.write(f"SCANS={scan['num']}\n")
            f.write(f"PEPMASS={scan['precursormz']}\n")
            f.write(f"CHARGE={scan['charge']}{_ion_mode}\n")
            f.write(f"COLLISION_ENERGY={scan.get('collision_energy', 0.0)}\n")
            for peak in scan['peaks']:
                f.write(f"{peak[0]}\t{peak[1]}\n")
            f.write("END IONS\n\n")


def write_metadata(output_file: Path, metadata: List[Dict]):
    metadata = pd.DataFrame(metadata)
    metadata.rename(columns={metadata.columns[0]: 'filename'}, inplace=True)
    metadata.columns = [col.replace(',', '').replace(';', '') for col in metadata.columns]
    metadata.columns = [f"ATTRIBUTE_{'_'.join(col.split())}" if col != 'filename' else 'filename' for col in metadata.columns]
    metadata.to_csv(output_file, sep='\t', index=False, header=True)


def reshape_compound_results(compound_results: List[Dict]) -> pd.DataFrame:
    df = pd.DataFrame(compound_results)
    df = df.rename(columns={'num': 'row ID', 'Mass': 'row m/z', 'RT': 'row retention time'})
    df = df.pivot_table(
        index=['row ID', 'row m/z', 'row retention time'],
        columns='File_Name',
        values='Response'
    ).reset_index()
    df.columns = [
        f"{col} Peak area" if col not in ['row ID', 'row m/z', 'row retention time'] else col
        for col in df.columns
    ]

    df.fillna(0.0, inplace=True)
    return df

def convert_to_feature_csv(df: pd.DataFrame, output_path: Path):
    """
    Convert the DataFrame to a CSV file with specific formatting.
    Args:
        df (pd.DataFrame): DataFrame to convert.
        output_path (Path): Path to save the CSV file.
    """
    df.to_csv(output_path, sep=',', index=False, header=True)

def process_pfa_file(input_file: Path, response_type: str, unzipped_files_dir: Path):
    """
    Process the PFA file and extract compound results, metadata, and MS2 scans to a directory.
    Args:
        input_file (Path): Path to the input PFA file.
        response_type (str): Type of response to extract ('Height' or 'Algo').
        unzipped_files_dir (Path): Directory to store unzipped files.

    Returns:
        pd.DataFrame: Reshaped compound results.
        List[Dict]: Metadata from the PFA file.
        List[Dict]: MS2 scans from the PFA file.
    """

    if not isinstance(input_file, Path):
        input_file = Path(input_file)
    if not isinstance(unzipped_files_dir, Path):
        unzipped_files_dir = Path(unzipped_files_dir)

    if not unzipped_files_dir.exists():
        unzipped_files_dir.mkdir(parents=True, exist_ok=True)
    else:
        raise ValueError(f"Unzipped files directory {unzipped_files_dir} already exists. Please remove it.")

    pfa = ReadPFAFile(input_file)
    metadata = GetSampleGroupingFromPFA(pfa)
    compound_results = GetCompoundResultsFromPFA(pfa, unzipped_files_dir)
    ms2_scans = GetSCompositeMSMSFromPFA(pfa)

    compound_table = ExportCompoundResults(compound_results, ResponseType=response_type)
    reshaped_results = reshape_compound_results(compound_table)

    return reshaped_results, metadata, ms2_scans

def main():
    parser = argparse.ArgumentParser(description='Extract and write data from PFA file.')
    parser.add_argument('--input_file', type=str, required=True, help='Input PFA file')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory')
    parser.add_argument('--response_type', choices=['Height', 'Algo'], default='Height')
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG)
    logging.info("Script Arguments:")
    for arg in vars(args):
        logging.info("%s: %s", arg, getattr(args, arg))

    input_path = Path(args.input_file)
    output_dir = Path(args.output_dir)
    unzipped_path = Path('./unzipped_files')

    if not input_path.exists():
        print(f"Input file {input_path} does not exist.")
        sys.exit(1)
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    results_df, metadata, ms2 = process_pfa_file(input_path, args.response_type, unzipped_path)

    results_df.to_csv(output_dir / '_featuretable_reformated.csv')
    write_metadata(output_dir / '_merged_metadata.tsv', metadata)
    write_to_mgf(output_dir / '_specs_ms.mgf', ms2)


if __name__ == "__main__":
    main()
