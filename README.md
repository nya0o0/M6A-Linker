# m6anet-to-m6add

m6anet-to-m6add is a Python tool designed to bridge the gap between m6Anet and M6ADD in m6A modification research. This mini-project focuses on creating a script that takes m6Anet output and formats it for M6ADD input, ensuring compatibility between the two systems.

The tool processes m6Anet output files, which typically contain information such as transcript IDs, positions, and modification probabilities. It then transforms this data into a format compatible with M6ADD, including chromosome information, genomic positions, and confidence scores. The resulting output is a CSV file that can be directly used as input for M6ADD, facilitating seamless integration of m6A site predictions with disease association data.

### Installation

To install and set up the m6anet-to-m6add tool for development, follow these steps:

1. First, ensure you have Conda installed. Then, create a new environment with the necessary dependencies:

```bash
conda create -n m6anet-to-m6add python=3.8 pandas -c conda-forge
conda activate m6anet-to-m6add
```

2. Clone the repository and install the package in editable mode:

```bash
git clone https://github.com/yourusername/m6anet-to-m6add.git
cd m6anet-to-m6add
pip install -e .
```

This will install the package in editable mode, allowing you to make changes to the code and immediately see the effects without reinstalling.

### Usage

After installation, you can use the tool from the command line:

```bash
m6anet-to-m6add input_m6anet_file.csv output_m6add_file.csv
```

For developers, you can also import and use the main function in your Python scripts:

```python
from m6anet_to_m6add import convert_m6anet_to_m6add

convert_m6anet_to_m6add('input_file.csv', 'output_file.csv')
```
