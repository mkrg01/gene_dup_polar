# gene_dup_polar

Code repository to run analyses for Nishiguchi, T. and Ishikawa, A. "Convergent gene duplication in Arctic and Antarctic teleost fish."

[Snakemake](https://snakemake.github.io/) is used to manage the workflow.

## Workflow execution

To execute the workflow, Singularity must be installed on your system.

[Clone this repository](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository).

Navigate to the cloned repository, substituting `/path/to/repo` with the actual path on your system:

```
cd /path/to/repo
```

Pull [a snakemake image](https://hub.docker.com/r/snakemake/snakemake):

```
singularity pull docker://snakemake/snakemake:v7.32.4
```

Run the workflow:

```
singularity run snakemake_v7.32.4.sif snakemake --use-singularity
```

The output files will be generated in the `data` folder.

## License

The code in this repository is licensed under the [MIT license](LICENSE).