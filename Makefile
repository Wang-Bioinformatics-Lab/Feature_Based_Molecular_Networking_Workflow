run:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config \
	--input_raw_spectra ./data/raw_data \
	--input_supplemental_edges ./data/raw_data

run_docker_test:
	nextflow run ./nf_workflow.nf -resume -c nextflow_dockertest.config \
	--input_raw_spectra ./data/raw_data \
	--input_supplemental_edges ./data/raw_data

run_no_metadata:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config --metadata_filename=""  \
	--input_raw_spectra ./data/raw_data \
	--input_supplemental_edges ./data/raw_data

run_agilent2:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config \
		--inputfeatures=data/agilent/new_version/Compound_Groups.tsv \
		--inputspectra=data/agilent/new_version/spectra/ \
		--featurefindingtool=AGILENT2  \
		--input_raw_spectra ./data/raw_data \
		--input_supplemental_edges ./data/raw_data
