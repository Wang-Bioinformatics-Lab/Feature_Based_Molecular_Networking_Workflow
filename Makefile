run:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config

run_no_metadata:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config --metadata_filename=""

run_agilent2:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config \
		--inputfeatures=data/agilent/new_version/Compound_Groups.tsv \
		--inputspectra=data/agilent/new_version/spectra/ \
		--featurefindingtool=AGILENT2
