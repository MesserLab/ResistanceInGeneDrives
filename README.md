# Reducing resistance allele formation in CRISPR/Cas9 gene drives

Files and scripts used to simulate and plot the evolutionary dinameics of resistance
allele formation in CRISPR gene drives.

**Citation:** Jackson Champer, Jingxian Liu, Suh Yeon Oh, Riona Reeves, Anisha
Luthra, Nathan Oakes, Andrew G. Clark, Philipp W. Messer. (2017). Reducing
resistance allele formation in CRISPR gene drives
**doi** https://doi.org/10.1101/150276


# Requirements
	python 3
		matplotlib
		pandas
		numpy

	SLiM (https://messerlab.org/slim/) in your PATH
		to append to a PATH using the console:
			export PATH=$PATH:/path/to/your/slim/directory

---------------------
# Examples

	Execute a single simulation run:
		slim Scripts/example_SLiM_inputs/autosome_2gRNA_high.slim

	Generate data for 1000 simulations:
		Each simulation has 5 rows: allele frequencies of Driver allele, Wildtype, R01, R10, R11
			Each row has 40 generations of allele frequencies
		./Scripts/generate_autosome_2gRNA_data.py

	Generate Figures from example data:
		./Scripts/generate_figures.py -example

	Generate Tiled Figure from example data:
		./Scripts/generate_tiled_figure.py -example
