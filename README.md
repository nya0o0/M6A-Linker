# m6alinker

M6A-Linker is a Python tool designed to bridge the gap between [**m6Anet**](https://m6anet.readthedocs.io/en/latest/) and [**M6ADD**](http://m6add.edbc.org/)，creating a script that takes m6Anet output and formats it for queries in the M6ADD database, ensuring compatibility between the two systems.

### Features

- Annotate **m6Anet** output into a format compatible with **M6ADD**.
- Extracts **transcript IDs, positions, and modification probabilities**.
- Maps transcript positions to **genomic coordinates**, including **chromosome information**.
- Outputs a **CSV file** that can inclues detailed information required by input for **M6ADD**.
- Efficient processing with **multi-threaded progress tracking**.

### Installation

To install and set up the m6anet-to-m6add tool for development, follow these steps:

```bash
git clone https://github.com/nya0o0/m6alinker
cd  m6alinker
pip install -e .
```

### Dependencies
- `pandas`
- `numpy`
- `gffutils`
- `threading`
- `warnings`

These dependencies should be installed automatically when you install the package.

### Usage

After installation, you can first download the refenrence gtf file form https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz to the `tests/` directory, decompress it 

```
gzip -d gencode.v47.annotationgtf.gz
```

and use the tool from the command line to run a data annotation for testing:

```bash
m6alinker -I tests/test_data_m6anet_mod_sites.csv -G tests/gencode.v47.annotation.gtf -O output_perfix
```

This will generate an output file named:`my_annotated_results_m6anet_results.csv` to the current pathway.

#### Required Arguments:
| Option | Description |
|---------|-------------|
| `-I` or `--input_file` | Path to the input file (CSV output from **m6Anet**). |
| `-G` or `--gtf_file` | Path to the **GTF** annotation file (reference genome). |

#### Optional Arguments:
| Option | Description | Default |
|---------|-------------|---------|
| `-O` or `--output_prefix` | Prefix for the output file. | `annotated_output` |

