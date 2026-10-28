# Running the lab notebooks

Each chapter ends with a lab notebook (`lab-*.ipynb`). You can run them in any of three places.

## 1. Georgia Tech instructional server (recommended)

Use the **launch** (rocket) button at the top of a lab page and choose *JupyterHub*.
The notebook and its data are copied into your account on the course server, where Python, R, PLINK and BWA
are already installed. Re-launching later does not overwrite your edits.

## 2. Binder (no account needed)

Choose *Binder* from the same launch button. Binder builds a temporary environment from the
repository's `environment.yml`. It can take a few minutes to start, and **your work is not saved**,
so download your notebook before closing the tab.

## 3. Your own computer

Install [Miniforge](https://github.com/conda-forge/miniforge), then:

```bash
git clone https://github.com/ikingjordan3/gab-book.git
cd gab-book
conda env create -f environment.yml
conda activate gab
jupyter lab
```

## Languages in the labs

Labs run on the Python kernel. R code runs in cells that start with `%%R`, and command-line tools
(PLINK, BWA, samtools) run in cells that start with `%%bash`, so one notebook can mix all three.
