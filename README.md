# TEAnnotate

TEAnnotate is a Nextflow workflow for annotating transposable elements (TEs) in animal, plant and fungal genome assemblies. It builds a consensus library, records the evidence and provenance for each consensus, and uses the library to produce a soft-masked assembly and repeat annotations.

## Workflow

| Input lineage | Discovery tools |
| --- | --- |
| Animal | Earl Grey and HiTE |
| Plant | EDTA, HELIANO and HiTE |
| Fungus | Earl Grey and HiTE |

TEAnnotate combines the discovery libraries and refines clusters that mix TE classes. Libraries with at least 20,000 consensi use a both-strand MMseqs2 search with greedy assignment by default. Other libraries use CD-HIT-EST. The evidence gate uses protein domains, support from another tool, and genomic copy counts. TEtrimmer curation is optional. Domain-based relabelling and a two-step classifier vote resolve some of the remaining unknown consensi. By default the vote can assign DNA transposon and Helitron classes. RepeatMasker reports interspersed and simple repeats together. Fungal runs can also report Starfish captains and repeat-induced point mutation (RIP) signatures.

## Requirements

- Linux with Nextflow 24.04.0 or later and a Java version supported by that Nextflow release.
- Conda or Mamba for tools packaged as environments.
- Apptainer or Singularity for the recommended container profile. Docker and Podman profiles are also available. A Conda-only run must use `--skip_hite`.
- Network access for initial downloads, or pre-populated container, Conda and asset caches. Allow space for the Nextflow work directory and genome masking outputs.

## Run an assembly

Save a CSV such as `samplesheet.csv` and replace the paths with your genome assemblies. Each sample ID must be unique. The organism types are `animal`, `plant`, and `fungus`. An empty type uses the animal route.

```csv
sample,genome,organism_type
plant1,/data/plant1.fna,plant
fungus1,/data/fungus1.fna,fungus
```

Launch from the TEAnnotate directory. If the copy was transferred from Windows, restore executable permissions on Linux with `chmod +x bin/*`.

```bash
nextflow run . -profile singularity --input samplesheet.csv --outdir results
```

On Slurm, use `-profile singularity,slurm`. Configure queue, resources and shared cache paths in your site configuration. Use Nextflow's `-resume` option to reuse successful tasks when their inputs and settings are unchanged.

| Option | Purpose |
| --- | --- |
| `--skip_tetrimmer` | Skip automated consensus curation |
| `--skip_reclassification` | Skip the classifier ensemble |
| `--user_te_library <fasta>` | Include a curated consensus library with tracked provenance |
| `--skip_masking` | Produce the library without genome masking |
| `--skip_starfish`, `--skip_rip` | Skip either fungal report |
| `--dedup_mmseqs_min_seqs <n>` | Change the large-library threshold. Set it to zero to use CD-HIT-EST throughout |

Parameters and defaults are defined in `nextflow.config` and `nextflow_schema.json`. Use `run_reclassify.nf` to start from an existing raw library. The entry point accepts the raw FASTA via `--input`, the assembly via `--genome`, the sample ID via `--id`, and the lineage via `--organism_type`. Supply `--dedup` to cluster a raw multi-tool library and `--provenance` to pass its tool provenance table.

## Read the outputs

- `01_genome_stats/` and `02_denovo_te/` contain assembly statistics, discovery libraries and their provenance.
- `03_library/` and `04_curation/` contain refined and gated libraries, evidence tables and curation records.
- `05_classification/<sample>.te_library.fa` is the final consensus library used for masking. This stage also records classifier decisions.
- `06_masking/` and `07_visualization/` contain the masked assembly, repeat coordinates and summary plots when enabled.
- `pipeline_info/` contains software versions and run metadata.

## Test the installation

The test profile uses a small bundled fungal assembly. A stub run checks workflow wiring and may still pull environments or images. It produces placeholders for annotation outputs.

```bash
nextflow run . -profile test,singularity -stub --outdir results_stub
nextflow run . -profile test,singularity --outdir results_test
python3 tests/check_versions.py results_test/pipeline_info/software_versions.yml
python -m pytest tests/unit
```

The full run exercises the fungal route on a small input. The smoke checks
require Python 3 with PyYAML. To validate the all-route stub, run the workflow
with the bundled samplesheet and user library, then run the output checker.

~~~bash
nextflow run . -profile test,singularity -stub \
    --input tests/stub/samplesheet_all_routes.csv \
    --user_te_library tests/stub/user_library.fa --outdir results_all_routes
bash tests/stub/check_outputs.sh results_all_routes \
    tests/stub/samplesheet_all_routes.csv . tests/stub/user_library.fa
~~~

Set TEANNOTATE_TEST_PYTHON to the Python executable with PyYAML if it differs
from python3.

For expected test outputs, see [the test dataset notes](tests/test_data/README.md).

### Tested setup

The `test,singularity` smoke test was run on Linux with Nextflow 24.10.6,
Java 11, and Apptainer 1.2.5. See [Troubleshooting](#troubleshooting) for the
Zhores setup used to complete the test with fresh caches.

## Troubleshooting

### Zhores fresh-cache test

Apptainer 1.2.5 returned HTTP 403 while pulling the Terrier image from GitHub
Container Registry. The test passed after the public Terrier, CREATE, and
NeuralTE images were downloaded from the registry and verified against their
published SHA-256 checksums. Place the verified images in the cache selected
by `--singularity_cache_dir` using the filenames shown in
[the container prefetch script](bin/prepull_containers.sh).

Use node-local scratch for `APPTAINER_TMPDIR` and set `PYTHONNOUSERSITE=1`
before launching Nextflow. The test used micromamba 2.3.3 with
`conda.useMicromamba = true` in its Nextflow site configuration.

## License and benchmark reproducibility

TEAnnotate source is distributed under the [MIT License](LICENSE). Third-party tools, databases, model weights, containers and test data retain their own terms. The [benchmark repository](https://github.com/lesidon/TEAnnotate_benchmark) contains the workflow snapshots and maize-specific scripts used to produce the manuscript results. Use those snapshots when reproducing the reported benchmark.
