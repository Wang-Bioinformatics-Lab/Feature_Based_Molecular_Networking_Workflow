#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the input channels with default values
params.featurefindingtool = "MZMINE"

params.inputfeatures = "data/mzmine2/gnps_featurefinding/features_quant.csv"
params.inputspectra = "data/mzmine2/gnps_featurefinding/spectra"

// Additional IIN Edges
params.input_supplemental_edges = ""

// Raw Data Including
params.input_raw_spectra = ""

// Normalization
params.normalization = "None" // Can also be RowSum, None

// Metadata
params.metadata_filename = "data/mzmine2/gnps_featurefinding/metadata.tsv"

// Libraries
params.input_libraries = "data/library"

params.pm_tolerance = "2.0"
params.fragment_tolerance = "0.5"

// Filtering Parameters
params.min_peak_intensity = "0.0" // This is a no op now
params.precursor_filter = "1"
params.window_filter = "1"

// Molecular Networking Options
params.similarity = "gnps"

params.parallelism = 24
params.networking_min_matched_peaks = 6
params.networking_min_cosine = 0.7
params.networking_max_shift = 1000

// Topology Filtering
params.topology_topk = 10
params.topology_maxcomponent = 100

// Library Search Parameters
params.library_topk = 1

params.library_min_cosine = 0.7
params.library_min_matched_peaks = 6

//TODO: Implement This
params.library_filter_precursor = 1
params.library_filter_window = 1

//TODO: Implement This
params.library_analog_search = "0"
params.library_analog_max_shift = 1999

// Define the paths to the required YAML files
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

// Set the path to the tool folder

params.publishdir = "$baseDir"
TOOL_FOLDER = "$baseDir/bin"


process filesummary {
    publishDir "$params.publishdir/nf_output/filesummary", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    errorStrategy 'ignore'

    input:
    file inputSpectra
    val ready

    output:
    file 'summaryresult.tsv'

    """
    python $TOOL_FOLDER/scripts/filesummary.py $inputSpectra summaryresult.tsv $TOOL_FOLDER/binaries/msaccess
    """
}

// Define the process that will reformat the quantification table
process quantification_table_reformatted {

    publishDir "$params.publishdir/nf_output/clustering", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_features 
    file input_spectra 

    output:
    path "featuretable_reformated.csv"
    path "specs_ms.mgf"

    script:
    """
    python $TOOL_FOLDER/scripts/reformat_quantification.py \
    $params.featurefindingtool \
    $input_features \
    featuretable_reformated.csv \
    $input_spectra \
    specs_ms.mgf \
    --QUANT_FILE_NORM $params.normalization
    """
}


process filter_spectra{

    publishDir "$params.publishdir/nf_output/clustering", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path _reformatted_features_ch 
    path _reformatted_spectra_ch

    output:
    path 'spectra_reformatted.mgf'
    path 'spectra_filtered.mgf'

    script:
    """
    python $TOOL_FOLDER/scripts/filter_spectra.py \
    --FILTER_PRECURSOR_WINDOW ${params.precursor_filter} \
    --WINDOW_FILTER ${params.window_filter} \
    $_reformatted_spectra_ch \
    spectra_reformatted.mgf \
    spectra_filtered.mgf 
    """
}


// Molecular Networking
process networkingGNPSPrepParams {
    publishDir "$params.publishdir/nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file spectrum_file

    output:
    file "params/*"

    """
    mkdir params
    python $TOOL_FOLDER/scripts/prep_molecular_networking_parameters.py \
        "$spectrum_file" \
        "params" \
        --parallelism "$params.parallelism" \
        --min_matched_peaks "$params.networking_min_matched_peaks" \
        --ms2_tolerance "$params.fragment_tolerance" \
        --pm_tolerance "$params.pm_tolerance" \
        --min_cosine "$params.networking_min_cosine" \
        --max_shift "$params.networking_max_shift"
    """
}

process calculatePairs {
    publishDir "$params.publishdir/nf_output/temp_pairs", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file spectrum_file
    each file(params_file)

    output:
    file "*_aligns.tsv" optional true

    """
    $TOOL_FOLDER/binaries/main_execmodule \
        ExecMolecularParallelPairs \
        "$params_file" \
        -ccms_INPUT_SPECTRA_MS2 $spectrum_file \
        -ccms_output_aligns ${params_file}_aligns.tsv
    """
}

// Filtering the network
process filterNetwork {
    publishDir "$params.publishdir/nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_pairs

    output:
    file "filtered_pairs.tsv"

    """
    python $TOOL_FOLDER/scripts/filter_networking_edges.py \
    $input_pairs \
    filtered_pairs.tsv \
    filtered_pairs_old_format.tsv  \
    --top_k_val $params.topology_topk \
    --max_component_size $params.topology_maxcomponent
    """
}

// Creating the metadata file
// TODO: Finish this
process createMetadataFile {
    publishDir "$params.publishdir/nf_output/metadata", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_metadata

    output:
    file "merged_metadata.tsv"

    //script in case its NO_FILE
    """
    python $TOOL_FOLDER/scripts/merge_metadata.py \
    $input_metadata \
    merged_metadata.tsv
    """
}

// Calculating the groupings
process calculateGroupings {
    publishDir "$params.publishdir/nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_metadata
    file input_featuretable

    output:
    file "clustersummary_with_groups.tsv"

    """
    python $TOOL_FOLDER/scripts/group_abundances.py \
    $input_featuretable \
    $input_metadata \
    clustersummary_with_groups.tsv
    """
}

// Library Search
process librarySearchData {
    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    each file(input_library)
    each file(input_spectrum)

    output:
    file 'search_results/*' optional true

    """
    mkdir search_results
    python $TOOL_FOLDER/scripts/library_search_wrapper.py \
    $input_spectrum $input_library search_results \
    $TOOL_FOLDER/binaries/convert \
    $TOOL_FOLDER/binaries/main_execmodule.allcandidates \
    --pm_tolerance $params.pm_tolerance \
    --fragment_tolerance $params.fragment_tolerance \
    --topk $params.library_topk \
    --library_min_cosine $params.library_min_cosine \
    --library_min_matched_peaks $params.library_min_matched_peaks \
    --analog_search $params.library_analog_search
    """
}

process librarymergeResults {
    publishDir "$params.publishdir/nf_output/library_intermediate", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path "results/*"

    output:
    path 'merged_results.tsv'

    """
    python $TOOL_FOLDER/scripts/tsv_merger.py \
    results merged_results.tsv \
    --topk $params.library_topk \
    --key_column "#Scan#" \
    --sort_column MQScore
    """
}

process librarygetGNPSAnnotations {
    publishDir "$params.publishdir/nf_output/library", mode: 'copy'

    //cache 'lenient'
    cache 'false'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path "merged_results.tsv"
    path "library_summary.tsv"

    output:
    path 'merged_results_with_gnps.tsv'

    """
    python $TOOL_FOLDER/scripts/getGNPS_library_annotations.py \
    merged_results.tsv \
    merged_results_with_gnps.tsv \
    --librarysummary library_summary.tsv
    """
}


// Enriching the Cluster Summary
process enrichClusterSummary {
    publishDir "$params.publishdir/nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_clustersummary
    file input_filtered_pairs
    file input_library_matches

    output:
    file "clustersummary_with_network.tsv"

    """
    python $TOOL_FOLDER/scripts/enrich_cluster_summary.py \
    $input_clustersummary \
    $input_filtered_pairs \
    $input_library_matches \
    clustersummary_with_network.tsv
    """
}

process createNetworkGraphML {
    publishDir "$params.publishdir/nf_output/networking", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false

    input:
    file input_clustersummary
    file input_filtered_pairs
    file input_library_matches
    file input_supplemental_edges

    output:
    file "network.graphml"
    file "network_singletons.graphml"

    """
    python $TOOL_FOLDER/scripts/create_network_graphml.py \
    $input_clustersummary \
    $input_filtered_pairs \
    $input_library_matches \
    $input_supplemental_edges \
    network.graphml \
    network_singletons.graphml
    """
}

process summaryLibrary {
    publishDir "$params.publishdir/nf_output", mode: 'copy'

    cache 'lenient'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path library_file

    output:
    path '*.tsv' optional true

    """
    python $TOOL_FOLDER/scripts/library_summary.py \
    $library_file \
    ${library_file}.tsv
    """
}

process createTallRawData {
    publishDir "$params.publishdir/nf_output/clustering", mode: 'copy'

    errorStrategy 'ignore'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path input_reformatted_quant
    path input_raw_spectra


    output:
    path 'tall_raw_data.tsv'

    """
    python $TOOL_FOLDER/scripts/create_tall_quant.py \
    $input_reformatted_quant \
    $input_raw_spectra \
    tall_raw_data.tsv
    """
}

workflow {
    // File Summary
    input_spectra_ch = Channel.fromPath(params.input_raw_spectra)
    filesummary(input_spectra_ch, 1)

    // Converting the quantification table
    input_features = Channel.fromPath(params.inputfeatures)
    input_spectra = Channel.fromPath(params.inputspectra)
    (_features_reformatted_ch, _spectra_reformatted_ch)  =  quantification_table_reformatted(input_features, input_spectra)

    // Filter the spectra
    (_spectra_reformatted_ch2, _spectra_filtered_ch) = filter_spectra(_features_reformatted_ch, _spectra_reformatted_ch)

    // Library Search
    libraries_ch = Channel.fromPath(params.input_libraries + "/*.mgf" )
    search_results_ch = librarySearchData(libraries_ch, _spectra_filtered_ch)
    merged_results_ch = librarymergeResults(search_results_ch.collect())
    merged_results_ch = merged_results_ch.ifEmpty(file("NO_FILE"))

    // Lets create a summary for the library files
    library_summary_ch = summaryLibrary(libraries_ch)

    // Merging all these tsv files from library_summary_ch within nextflow
    library_summary_merged_ch = library_summary_ch.collectFile(name: "library_summary.tsv", keepHeader: true)
    library_summary_merged_ch = library_summary_merged_ch.ifEmpty(file("NO_FILE"))

    gnps_library_results_ch = librarygetGNPSAnnotations(merged_results_ch, library_summary_merged_ch)
    gnps_library_results_ch = gnps_library_results_ch.ifEmpty(file("NO_FILE"))

    // Networking
    params_ch = networkingGNPSPrepParams(_spectra_filtered_ch)
    networking_results_temp_ch = calculatePairs(_spectra_filtered_ch, params_ch.collect())

    merged_networking_pairs_ch = networking_results_temp_ch.collectFile(name: "merged_pairs.tsv", storeDir: "./nf_output/networking", keepHeader: true)

    // Filtering the network
    filtered_networking_pairs_ch = filterNetwork(merged_networking_pairs_ch)

    // Handling Metadata, if we don't have one, we'll set it to be empty
    if(params.metadata_filename.length() > 0){
        if(params.metadata_filename == "NO_FILE"){
            input_metadata_ch = Channel.of(file("NO_FILE"))
        }
        else{
            input_metadata_ch = Channel.fromPath(params.metadata_filename).first()
        }
    }
    else{
        input_metadata_ch = Channel.of(file("NO_FILE"))
    }

    merged_metadata_ch = createMetadataFile(input_metadata_ch)

    // Enriching the network with group mappings
    clustersummary_with_groups_ch = calculateGroupings(merged_metadata_ch, _features_reformatted_ch)

    // // Adding component and library informaiton
    clustersummary_with_network_ch = enrichClusterSummary(clustersummary_with_groups_ch, filtered_networking_pairs_ch, gnps_library_results_ch)

    // Supplemental Edges
    supplemental_edges_ch = Channel.fromPath(params.input_supplemental_edges)

    // // Creating the graphml Network
    createNetworkGraphML(clustersummary_with_network_ch, filtered_networking_pairs_ch, gnps_library_results_ch, supplemental_edges_ch)

    // Creating the tall quant table
    createTallRawData(_features_reformatted_ch, input_spectra_ch)

}