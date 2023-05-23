run:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config

run_no_metadata:
	nextflow run ./nf_workflow.nf -resume -c nextflow.config --metadata_filename=""

run_hpcc:
	nextflow run ./nf_workflow.nf -resume -c nextflow_hpcc.config

run_docker:
	nextflow run ./nf_workflow.nf -resume -with-docker <CONTAINER NAME>