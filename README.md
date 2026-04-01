# AREA

Attribute Rank Enrichment Analysis

## Summary

The goal of AREA is to link boolean attributes to rankable attributes via samples. In this way, we may be able to find rankable attributes that cause the boolean attributes or vice versa. Similar to correlation analysis, causality is not clear from this analysis alone. Instead, these linkages should be considered potentially causal, and downstream experiments should be used to determine causality.

In our example, we analyze gene expression data from individuals. These individuals have Down syndrome, and the boolean attributes are their associated medical conditions. Our objective is to identify genes that, when highly expressed, may influence the likelihood of specific medical conditions.

![results example](https://github.com/Dowell-Lab/psea/blob/main/src/images/results_example_NES.png "results example")

NES stands for Normalized Enrichment Score. A negative NES indicates that low levels of the rankable attribute (in our case gene expression) are associated with the true condition of the boolean attribute. Conversely, a positive NES suggests that high levels of the rankable attribute are associated with the true condition of the boolean attribute.

## Installation

AREA requires Python 3.9 or later.

### Setting up a virtual environment (SLURM supercomputer)

If you are on a shared supercomputer with SLURM, create a dedicated virtual environment for AREA:

```bash
module load python/3.9.15
mkdir -p ~/VENVs
cd ~/VENVs
python -m venv areavenv39
source areavenv39/bin/activate
pip install -r ~/AREA/requirements.txt
```

Activate this environment before each AREA run:

```bash
source ~/VENVs/areavenv39/bin/activate
```

### Install with pip

```bash
cd AREA
pip install -e .
```

For GPU acceleration via CuPy:

```bash
pip install -e ".[gpu]"
```

## AREA inputs and outputs

AREA takes two CSV files as input. Both CSVs share the same samples (in our example, individuals), but one has boolean attributes for each sample and the other has rankable attributes. The output is a CSV with the statistical significance of the linkages between the boolean columns and the rank columns.

### Input files

Both input CSV files must have a common sample column. In our case, that common sample name column is "Participant".

### Binary attribute file

One of the CSV files needs to contain boolean attributes. In our case, the boolean attributes are the disease or disorder associated with each patient, which we call a comorbidity.

![boolean attributes csv](https://github.com/Dowell-Lab/psea/blob/main/src/images/boolean_attributes_df.png "boolean attributes csv")

### Value file

The other file must have columns that can be ranked by the values within the column. In our example, the rank file has genes and the expression level of those genes in each patient.

![Value csv](https://github.com/Dowell-Lab/psea/blob/main/src/images/rank_df.png "Value csv")

### Output file

The output file lists each boolean attribute and rank column in pairs. The rest of the row has the scores for that pair, including NES, raw p-value, and four adjusted p-value columns (Bonferroni, Holm, Benjamini-Hochberg, and Benjamini-Yekutieli).

![results example](https://github.com/Dowell-Lab/psea/blob/main/src/images/results_example_NES.png "results example")

### Filtering the output

The items in the output file have NOT been filtered for significance. To filter for significance, pick an adjusted p-value column (four are provided) and apply your chosen cutoff.

## Running AREA

There are three ways to run AREA after installation.

### Option 1: Direct script

No install needed — run from inside the `AREA/` directory:

```bash
python run_area.py --config area_config.yaml
```

### Option 2: Installed command

After `pip install -e .`, the `area` command is available anywhere:

```bash
area --config area_config.yaml
```

### Option 3: Module syntax

```bash
python -m area --config area_config.yaml
```

### Using a config file

Copy `area_config.yaml` and fill in your paths:

```yaml
boolean_file: /path/to/bools.csv
rank_file: /path/to/ranks.csv
join_column: Participant
out_dir: /path/to/results/
threads: 4
gpu: false
verbose: false
```

Then run:

```bash
python run_area.py --config my_config.yaml
```

### Using command-line flags

All parameters can also be passed directly on the command line:

```bash
python run_area.py \
  -bf /path/to/bools.csv \
  -rf /path/to/ranks.csv \
  -jc Participant \
  -od /path/to/results/
```

### Mixing config file and CLI flags

Command-line flags always override values from the config file. This lets you set standard parameters in your config and override specific ones per run:

```bash
python run_area.py --config my_config.yaml --verbose --threads 8
```

### Full list of parameters

| Flag | Config key | Required | Description |
|------|-----------|----------|-------------|
| `-bf` / `--boolean-file` | `boolean_file` | Yes | CSV of binary attributes (samples x attributes) |
| `-rf` / `--rank-file` | `rank_file` | Yes | CSV of continuous values (samples x features) |
| `-jc` / `--join-column` | `join_column` | Yes | Column name shared between both input files |
| `-od` / `--out-dir` | `out_dir` | Yes | Directory for all output files |
| `-t` / `--threads` | `threads` | No | Number of parallel threads (default: 4) |
| `--gpu` | `gpu` | No | Use GPU acceleration via CuPy (default: false) |
| `--keep-rank-columns` | `keep_rank_columns` | No | Text file listing rank columns to include |
| `--keep-bool-columns` | `keep_bool_columns` | No | Text file listing bool columns to include |
| `--keep-samples` | `keep_samples` | No | Text file listing sample IDs to include |
| `--verbose` | `verbose` | No | Print diagnostic output (default: false) |
| `-c` / `--config` | — | No | Path to a YAML config file |

### Filtering inputs with keep files

Often a user will not want to run all data through AREA. If a boolean attribute is true for all samples or no samples, it will not produce meaningful results. Similarly, genes with no expression variance across samples are uninformative. AREA automatically excludes constant columns, but you can further limit which columns are analyzed using keep files:

- `--keep-bool-columns` — a text file with one boolean column name per line
- `--keep-rank-columns` — a text file with one rank column name per line
- `--keep-samples` — a text file with one sample ID per line

## Testing

AREA includes a test suite in `unittest/test_area.py` that covers both unit tests and integration tests.

### Test data

Place your test files in the `testdata/` directory:

- `testdata/genes.csv` — a rank file (samples x genes)
- `testdata/comorbid_file.csv` — a boolean file (samples x comorbidities)

Unit tests for the backend, config loader, enrichment math, and CLI parsing run without test data. The integration tests (planning, runner, correction, and the full end-to-end pipeline) require the test data files and will skip gracefully if they are not present.

### Running the tests

From the `AREA/` directory:

```bash
python run_tests.py
python run_tests.py -v   # verbose output
```

### What the tests cover

| Test class | Module tested | Requires test data |
|------------|--------------|--------------------|
| `TestBackend` | `backend.py` — CPU backend, `to_numpy`, `trapz` helper | No |
| `TestConfig` | `config.py` — YAML loading, unknown-key rejection | No |
| `TestEnrichment` | `enrichment.py` — enrichment score sign, permutation count, NES p-value range | No |
| `TestCLI` | `cli.py` — argument parsing, missing-required errors, config/CLI override precedence | No |
| `TestPlanning` | `planning.py` — plan file creation, row structure, constant-column exclusion | Yes |
| `TestRunner` | `runner.py` — raw p-values output, column validation, p-values in [0, 1] | Yes |
| `TestCorrection` | `correction.py` — adjusted file creation, Bonferroni/Holm/BH/BY columns, sort order | Yes |
| `TestEndToEnd` | Full pipeline via `cli.main()` with `--verbose` and 2 threads | Yes |

## Project structure

```
AREA/
├── pyproject.toml           # Package metadata and dependencies
├── run_area.py              # Top-level entry script
├── run_tests.py             # Top-level test runner
├── area_config.yaml         # Documented config template
├── examples/                # Example SLURM shell scripts
├── testdata/                # Test CSV files (not in repo)
│   ├── genes.csv
│   └── comorbid_file.csv
├── unittest/                # Test suite
│   └── test_area.py
└── src/
    └── area/
        ├── __init__.py      # Package docstring
        ├── __main__.py      # Enables "python -m area"
        ├── backend.py       # GPU/CPU array backend selection
        ├── enrichment.py    # Core math (enrichment score, permutations, NES)
        ├── planning.py      # Run plan construction and column filtering
        ├── runner.py        # Thread-pool orchestration
        ├── correction.py    # Multiple-testing p-value correction
        ├── config.py        # YAML config file loader
        └── cli.py           # Argument parsing and main() entry point
```

## Simulated data

To discover what patterns AREA is best at finding, we created simulated data similar to our real data in `notebook_examples/simulateddata-bothdirs.ipynb`. This code creates a certain number of simulated genes based on the expression of real genes. We then create comorbidities that are biased in a sub-group of people with higher or lower expression.

## Special thank you to the Human Trisomy Project and the INCLUDE Data Coordinating Center

The data used in our example code can be sourced from the INCLUDE Data Hub at https://portal.includedcc.org/ and anyone who would like to use the original data can register for an account on the data hub. Data used from the Human Trisomy Project was converted from kallisto to DESeq2 normalized counts using tximport.
