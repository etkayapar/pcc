<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [PCC (Phylogenetic dataset Compiler Collection)](#pcc-phylogenetic-dataset-compiler-collection)
  - [Psyche analysis instructions for the CSC puhti cluster](#psyche-analysis-instructions-for-the-csc-puhti-cluster)
    - [Creating a fake `conda` command](#creating-a-fake-conda-command)
    - [Double-checking the profiles before attempting to run the pipeline](#double-checking-the-profiles-before-attempting-to-run-the-pipeline)
    - [Start the pipeline](#start-the-pipeline)
  - [Troubleshooting](#troubleshooting)
    - [I get a WorkflowError, but no sign of slurm log paths in the screen log](#i-get-a-workflowerror-but-no-sign-of-slurm-log-paths-in-the-screen-log)
    - [Conda is not found](#conda-is-not-found)

<!-- markdown-toc end -->

# PCC (Phylogenetic dataset Compiler Collection)

## Psyche analysis instructions for the CSC puhti cluster

This branch includes some changes and a helper script to get the pipeline running with the snakemake SLURM executors for maximum parallel efficiency possible.

First of all you should pull all the newest changes from upstream as this will give you the new branch specific for CSC:

```bash
git pull
```

If you see an error about conflicting or unsaved changes then it must be because you have changed some files that are tracked by git. Could be the Snakefile, and the rule definitions under `rules/` or configs. If you made changes that you wish not to lose, you can copy the changed files to a different directory and then `git pull`, or you can safely stash away those changes with git:

```bash
git stash push
```

this should make all the changes go away but make them accessible later on if needed (with the `git stash pop` command). After removing the changes and `git pull` ing, you need to switch to this new branch:

```bash
git switch psyche-csc
```

Since conda is not directly available on the cluster you can use the provided container image that provides the conda environments and conda itself. Due to a bug in Snakemake, the container image needs to be pulled first before attempting it to run the workflow. To do so, you need to run the following command after loading the snakemake module provided by the cluster ( `ml load snakemake`)

```bash
snakemake --sdm conda apptainer --conda-create-envs-only
```

I tried running this step inside a SLURM job via a batch script but it failed to due to not having enough space on disk in the temporary directory, but running it inside the login node was successful.

After pulling the container successfully, you should see a singularity/apptainer container image (`*.simg` file) somewhere in the `.snakemake/singularity/` directory under your working directory.

### Creating a fake `conda` command

We need to create a fake `conda` command that can provide information about the conda installation that exists inside the container image. I put a helper script to that in this branch:

If you run the below command while still being in the top-level `pcc` directory,

```bash
utils/create_dummy_conda.sh auto
```

It should create such a fake `conda` script that outputs a json string for snakemake to parse. There should be adequate information in the output of this above script to tell if it worked, but to explicity test if it worked try running:

```bash
type conda
```

this should output a path that looks like `/users/USERNAME/.local/bin/conda`.


### Double-checking the profiles before attempting to run the pipeline

It is a bit tedious to set up the profiles so that you are not over- or under-requesting resources for the partition of your choice unfortunately...

But please do check the `.yaml` files under both `profiles/default` and `profiles/psyche-slurm` (the default and the slurm profiles hereafter) to make sure you have reasonable times and memory there. I tried to set them up so that it is at least somewhat reasonable for the `small` partition on puhti, so they may be plug and play for you except the missing SLURM account information you need to fill in in the `.yaml`  file for the slurm profile.

Since the resource specifications in the default profile seems to override everything else, please make sure that you are not asking for more time than your partition allows. Below is the excerpt form the default profile I set up for the `small` partition.

```yaml
set-resources:
  align_aa:
    runtime: "3d"
  realign_outliers_before_trimal:
    runtime: "3d"
  realign_outliers_after_trimal:
    runtime: "3d"
  infer_gene_trees_before_trimal:
    runtime: "3d"
  infer_gene_trees_after_trimal:
    runtime: "3d"
  final_gene_trees:
    runtime: "3d"
```

Since this partition has a maximum 3-day runtime I gave all these rules that many days.

### Start the pipeline

For this I recommend starting a `tmux` or `screen` session on the login node you are now and try to remember which login front-end you are connected since you need to login back to this specific login node to be able to check on your running Snakemake process.

Inside the tmux session load the snakemake module provided by your cluster, as before

```bash
ml load snakemake
```

and make sure that you have the latest available snakemake module loaded.

now we can start the pipeline:

```bash
snakemake --sdm conda apptainer --profile profiles/psyche-slurm  --apptainer-args='--bind="/users,/projappl,/scratch,$TMPDIR,$LOCAL_SCRATCH"' --local-storage-prefix='"$LOCAL_SCRATCH"' --remote-job-local-storage-prefix='"$LOCAL_SCRATCH"' collect_gene_trees_before_trimal
```
Note: For my last successful runs, I actually typed out the value I have for `$LOCAL_SCRATCH` instead of passing it by the variable (as I shown above) for the two `storage-prefix=` arguments in the above command. So you may try both if one does not work because of quotation or some other reason.

## Troubleshooting

### I get a WorkflowError, but no sign of slurm log paths in the screen log

This is usually because the executor was not even able to submit jobs because of invalid resource specifications. Please check that the time you are requesting for a job (in the default and the slurm profiles) does not exceed the maximum time allowed for the partition you are submitting to, which can be changed from the file #2.

### Conda is not found

If your Snakemake process complains that it could not find conda then make sure that you successfully completed the fake conda [section](#creating-a-fake-conda-command)

