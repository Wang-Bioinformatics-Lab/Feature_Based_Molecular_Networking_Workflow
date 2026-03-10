#!/usr/bin/env python
# coding: utf-8

"""
Export GNPS Feature-Based Molecular Networking data for MetaboAnalyst analysis.
Supports two export formats:
1. Single file format (one factor) - metadata and features combined
2. Two file format (metadata table) - separate metadata and quantification tables
"""

import argparse
import os
import sys
import pandas as pd


def export_one_factor(ftable_path, metadata_path, output_path):
    """
    Export combined metadata and feature table for MetaboAnalyst one-factor analysis.
    
    Args:
        ftable_path: Path to input feature table CSV
        metadata_path: Path to input metadata TSV
        output_path: Path for output CSV file
    """
    # Read feature table and metadata
    ftable = pd.read_csv(ftable_path, sep=',')
    metadata = pd.read_csv(metadata_path, sep='\t', index_col='filename')


    
    # Round the m/z values to 4 decimal places and the rt values to 2 decimal places
    ftable['row m/z'] = ftable['row m/z'].round(4)
    ftable['row retention time'] = ftable['row retention time'].round(2)
    
    # Concatenate row id and row m/z to create a unique feature id
    ftable['feature_id'] = (ftable['row ID'].astype(str) + '_' + 
                           ftable['row m/z'].astype(str) + '_' + 
                           ftable['row retention time'].astype(str))
    
    # Keep only columns relative to samples and feature ids
    ftable = ftable[[col for col in ftable.columns if col.lower().endswith('Peak area') or col == 'feature_id']]
    
    # Move feature_id to first column
    cols = ftable.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    ftable = ftable[cols]
    
    # Remove 'Peak area' from column names
    ftable.columns = [col.replace(' Peak area', '') for col in ftable.columns]
    
    # Get the order of columns in ftable and reorganize metadata accordingly
    sample_order = ftable.columns.tolist()
    sample_order.remove('feature_id')
    
    # Reorder metadata rows to match
    metadata_reordered = metadata.reindex(sample_order)
    metadata_reordered.fillna('NA', inplace=True)  # Replace NaN with NA for MetaboAnalyst
    
    # Transpose metadata
    meta_T = metadata_reordered.T.reset_index()
    
    # Rename the first column to "Samples" in both tables
    meta_T.rename(columns={'index': 'Samples'}, inplace=True)
    ftable.rename(columns={'feature_id': 'Samples'}, inplace=True)
    
    # Stack metadata on top of features
    merged = pd.concat([meta_T, ftable])
    merged.index.name = "Samples"
    
    # Export merged table
    merged.to_csv(output_path, sep=',', index=False)
    print(f"One-factor format exported to: {output_path}")


def export_metadata_table(ftable_path, metadata_path, output_quant, output_metadata):
    """
    Export separate quantification and metadata tables for MetaboAnalyst metadata table analysis.
    
    Args:
        ftable_path: Path to input feature table CSV
        metadata_path: Path to input metadata TSV
        output_quant: Path for output quantification table CSV
        output_metadata: Path for output metadata CSV
    """
    # Read feature table and metadata
    ftable = pd.read_csv(ftable_path, sep=',')
    metadata = pd.read_csv(metadata_path, sep='\t', index_col='filename')
    
    #### Reformat quant table for MetaboAnalyst
    # Round the m/z values to 4 decimal places and the rt values to 2 decimal places
    ftable['row m/z'] = ftable['row m/z'].round(4)
    ftable['row retention time'] = ftable['row retention time'].round(2)
    
    # Concatenate row id and row m/z to create a unique feature id
    ftable['Sample'] = (ftable['row ID'].astype(str) + '_' + 
                       ftable['row m/z'].astype(str) + '_' + 
                       ftable['row retention time'].astype(str))
    
    # Keep only columns relative to samples and feature ids
    ftable = ftable[[col for col in ftable.columns if col.endswith('Peak area') or col == 'Sample']]
    
    # Move Sample to first column
    cols = ftable.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    ftable = ftable[cols]
    
    # Remove 'Peak area' from column names
    ftable.columns = [col.replace(' Peak area', '') for col in ftable.columns]
    
    # Remove .mzML or .mzXML from column names
    ftable.columns = [col.replace('.mzML', '').replace('.mzXML', '') for col in ftable.columns]
    
    #### Reformat metadata for MetaboAnalyst
    # Rename filename column to Sample
    metadata = metadata.reset_index()
    metadata = metadata.rename(columns={'filename': 'Sample'})
    
    # Remove .mzML or .mzXML from the Sample column
    metadata['Sample'] = metadata['Sample'].str.replace('.mzML', '').str.replace('.mzXML', '')
    
    # Replace NaN with NA for MetaboAnalyst
    metadata.fillna('NA', inplace=True)
    
    # Keep only the columns relative to samples with metadata
    valid_samples = set(metadata['Sample'])
    cols_to_keep = ['Sample'] + [col for col in ftable.columns if col in valid_samples]
    ftable_filtered = ftable[cols_to_keep]
    
    # Keep only metadata for samples in the feature table
    valid_samples = set(ftable_filtered.columns) - {'Sample'}
    metadata_filtered = metadata[metadata['Sample'].isin(valid_samples)]
    
    ### Export
    ftable_filtered.to_csv(output_quant, index=False)
    metadata_filtered.to_csv(output_metadata, index=False)
    print(f"Quantification table exported to: {output_quant}")
    print(f"Metadata table exported to: {output_metadata}")


def generate_readme(output_path):
    """
    Generate a README file with instructions for using the exported files.
    
    Args:
        output_path: Path for the output readme.txt file
    """
    readme_content = """MetaboAnalyst Export Instructions
===================================

These exports can be used as input in MetaboAnalyst for further statistical analysis.

ONE-FACTOR FORMAT:
- Should be used as input in Statistical Analysis [one factor]
- Make sure to specify which metadata row will be used prior to upload to MetaboAnalyst

METADATA TABLE FORMAT (quant + metadata files):
- Should be used as input in Statistical Analysis [metadata table]
- Imputing missing values in MetaboAnalyst might be necessary

For more information, visit: https://www.metaboanalyst.ca/
"""
    
    with open(output_path, 'w') as f:
        f.write(readme_content)
    print(f"README file generated: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Export GNPS FBMN data for MetaboAnalyst analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Export one-factor format
  %(prog)s -f feature_table.csv -m metadata.tsv --metaboanalyst-onefactor output.csv
  
  # Export metadata table format
  %(prog)s -f feature_table.csv -m metadata.tsv \\
    --metaboanalyst-quant quant.csv --metaboanalyst-metadata meta.csv
  
  # Generate README file
  %(prog)s -f feature_table.csv -m metadata.tsv \\
    --metaboanalyst-onefactor output.csv --readme readme.txt
        """
    )
    
    # Input arguments
    parser.add_argument('-f', '--feature-table', 
                       required=True,
                       help='Input feature table CSV file from GNPS')
    parser.add_argument('-m', '--metadata',
                       required=True,
                       help='Input metadata TSV file')
    
    # Output arguments - one-factor format
    parser.add_argument('--metaboanalyst-onefactor',
                       help='Output file for MetaboAnalyst one-factor format (single CSV)')
    
    # Output arguments - metadata table format
    parser.add_argument('--metaboanalyst-quant',
                       help='Output quantification table for MetaboAnalyst metadata table format')
    parser.add_argument('--metaboanalyst-metadata',
                       help='Output metadata table for MetaboAnalyst metadata table format')
    
    # README generation
    parser.add_argument('--readme',
                       help='Generate a README file with usage instructions')
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not os.path.exists(args.feature_table):
        print(f"Error: Feature table file not found: {args.feature_table}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.metadata):
        print(f"Error: Metadata file not found: {args.metadata}", file=sys.stderr)
        sys.exit(1)
    
    # Check that at least one output format is specified
    if not args.metaboanalyst_onefactor and not (args.metaboanalyst_quant and args.metaboanalyst_metadata):
        print("Error: Must specify at least one output format:", file=sys.stderr)
        print("  --metaboanalyst-onefactor OR", file=sys.stderr)
        print("  --metaboanalyst-quant + --metaboanalyst-metadata", file=sys.stderr)
        sys.exit(1)
    
    # Validate metadata table format requires both outputs
    if (args.metaboanalyst_quant and not args.metaboanalyst_metadata) or \
       (args.metaboanalyst_metadata and not args.metaboanalyst_quant):
        print("Error: Both --metaboanalyst-quant and --metaboanalyst-metadata must be specified together",
              file=sys.stderr)
        sys.exit(1)
    
    # Export one-factor format
    if args.metaboanalyst_onefactor:
        try:
            export_one_factor(args.feature_table, args.metadata, args.metaboanalyst_onefactor)
        except Exception as e:
            print(f"Error exporting one-factor format: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Export metadata table format
    if args.metaboanalyst_quant and args.metaboanalyst_metadata:
        try:
            export_metadata_table(args.feature_table, args.metadata, 
                                args.metaboanalyst_quant, args.metaboanalyst_metadata)
        except Exception as e:
            print(f"Error exporting metadata table format: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Generate README if requested
    if args.readme:
        try:
            generate_readme(args.readme)
        except Exception as e:
            print(f"Error generating README: {e}", file=sys.stderr)
            sys.exit(1)
    
    print("\nExport completed successfully!")


if __name__ == '__main__':
    main()




