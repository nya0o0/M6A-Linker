# m6anet-to-m6add

m6anet-to-m6add is a Python tool designed to bridge the gap between m6Anet and M6ADD，creating a script that takes m6Anet output and formats it for M6ADD input, ensuring compatibility between the two systems.

The tool processes m6Anet output files, which typically contain information such as transcript IDs, positions, and modification probabilities. It then transforms this data into a format compatible with M6ADD, including chromosome information, genomic positions, and confidence scores. The resulting output is a CSV file that can be directly used as input for M6ADD, facilitating seamless integration of m6A site predictions with disease association data.

### Installation

To install and set up the m6anet-to-m6add tool for development, follow these steps:

```bash
git clone https://github.com/yourusername/m6anet-to-m6add.git
cd m6anet-to-m6add
pip install -e .
```

### Usage

After installation, you can use the tool from the command line:

```bash
m6anet-to-m6add input_m6anet_file.csv output_m6add_file.csv
```
