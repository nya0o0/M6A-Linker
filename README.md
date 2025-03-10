# m6alinker

M6A-Linker is a Python tool designed to bridge the gap between m6Anet and M6ADD，creating a script that takes m6Anet output and formats it for M6ADD input, ensuring compatibility between the two systems.

The tool processes m6Anet output files, which typically contain information such as transcript IDs, positions, and modification probabilities. It then transforms this data into a format compatible with M6ADD, including chromosome information, genomic positions, and confidence scores. The resulting output is a CSV file that can be directly used as input for M6ADD, facilitating seamless integration of m6A site predictions with disease association data.

### Installation

To install and set up the m6anet-to-m6add tool for development, follow these steps:

```bash
git clone https://github.com/nya0o0/m6alinker
cd  m6alinker
pip install -e .
```

### Usage

After installation, you can use the tool from the command line to run a data annotation for testing:

```bash
m6alinker -I tests/test_data_m6anet_mod_sites.csv -G gencode.v32.annotation.gtf -O output_perfix
```
Options:
1. Required
	`-I` or `--input_file` The path to the input file (i.e., the output csv file from m6anet)
	`-G` or `--gtf_file` The path to the gtf files for reference.
2. Optional
	`-O` or `--output_prefix` The Prefix for the output file. (Defalt: annotated_output)

