#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the input channels with default values
params.featurefindingtool = "MZMINE"

params.inputfeatures = "./data/MZmine-GNPS_AG_test_featuretable.csv"
params.inputspectra = "data/inputSpectra"

// Define the paths to the required YAML files
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

// Set the path to the tool folder
TOOL_FOLDER = "$baseDir/bin"

// Define the process that will reformat the quantification table


// This process will reformat the input 
process quantification_table_reformatted {

    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_features 
    file input_spectra 

    output:
    path "featuretable_reformated.csv"
    path "specs_ms.mgf"

    script:
    """
    python $TOOL_FOLDER/reformat_quantification.py \
    $params.featurefindingtool \
    $input_features \
    featuretable_reformated.csv \
    $input_spectra \
    specs_ms.mgf
    """
}


process filter_spectra{

    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path _reformatted_features_ch 
    path _reformatted_spectra_ch

    output:
    path 'spectra_reformatted.mgf'
    path 'spectra_filtered.mgf'

    script:
    """
    python $TOOL_FOLDER/filter_spectra.py \
    --FILTER_PRECURSOR_WINDOW 0 \
    --WINDOW_FILTER 0 \
    $_reformatted_spectra_ch \
    spectra_reformatted.mgf \
    spectra_filtered.mgf 
    """
}


process prep_molecular_networking_parameters{

    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path spectra_filtered_ch 

    output:
    path 'networking_parameters_folder/*'

    script:
    """
    mkdir networking_parameters_folder

    python $TOOL_FOLDER/prep_molecular_networking_parameters.py \
    $spectra_filtered_ch \
    networking_parameters_folder \
    --parallelism 5 \
    --MIN_MATCHED_PEAKS 2 \
    --ion_tolerance 0.1 \
    --pm_tolerance 0.1 \
    --PAIRS_MIN_COSINE 0.5 \
    --MAX_SHIFT 0.1
    """
}


process molecular_networking_parallel_step {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    each path(networking_parameters)
    path(spectra_mgf)

    output:
    path '*.aligns'

    script:
    """
    $TOOL_FOLDER/main_execmodule \
    ExecMolecularParallelPairs \
    $networking_parameters \
    -ccms_output_aligns ${networking_parameters}.aligns \
    -ccms_INPUT_SPECTRA_MS2 $spectra_mgf
    """
}

process merge_networking_results {
    cache false 

    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file "aligns/*"

    output:
    path "merged_alignments.tsv"

    script:
    """
    python $TOOL_FOLDER/merge_tsv_files_efficient.py \
    aligns \
    merged_alignments.tsv
    """
}


process filter_networking_edges{

    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"


    input:
    path networking_pairs_results_file_ch

    output:
    path 'networking_pairs_results_file_filtered.tsv' 
    path 'networkedges_legacy_output.tsv'


    script:
    """
    python $TOOL_FOLDER/filter_networking_edges.py \
    $networking_pairs_results_file_ch \
    networking_pairs_results_file_filtered.tsv \
    networkedges_legacy_output.tsv
    """
}


// process merge_tsv_efficient{

//     publishDir "./nf_output", mode: 'copy'

//     conda "$TOOL_FOLDER/conda_env.yml"


//     input:
//     file input from Channel.fromPath(params.input)
//     val x from _val_ch4

//     output:
//     val 5 into _val_ch5


//     script:
//     """
//     python $TOOL_FOLDER/bin/merge_tsv_files_efficient.py \
//     $TOOL_FOLDER/test/ \
//     $TOOL_FOLDER/test/reference_stats/librarysearch.tsv
//     """
// }





workflow {
  
    def input_features = Channel.fromPath(params.inputfeatures)
    def input_spectra = Channel.fromPath(params.inputspectra)
    (_features_reformatted_ch, _spectra_reformatted_ch)  =  quantification_table_reformatted(input_features, input_spectra)

    // Filter the spectra
    (_spectra_reformatted_ch2, _spectra_filtered_ch) = filter_spectra(_features_reformatted_ch, _spectra_reformatted_ch)

    // Prepping the networking parameters
    // _networking_parameters_ch = prep_molecular_networking_parameters(_spectra_filtered_ch)

    // // 
    // _networking_pairs_results_ch = molecular_networking_parallel_step(_networking_parameters_ch.flatten(), _spectra_filtered_ch)

    // // Merging all the individual results into one file
    // _merge_results_ch = merge_networking_results(_networking_pairs_results_ch.collect())

    // // Filtering network edges
    // (_fitered_edges_ch, _legacy_edges_ch) = filter_networking_edges(_merge_results_ch)
}