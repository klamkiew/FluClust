# Influenza

## Environment recreation

### Automatic Install
- conda env create -f clustering.yml
- conda activate clustering
- jupyter lab

### Manual Install
- mamba create -n clustering jupyterlab jupyter-resource-usage jupyterlab-git jupyterlab-fasta jupyterlab-lsp pandas numpy scipy biopython hdbscan scikit-learn nodejs tqdm kneed rpy2 python-lsp-server r-languageserver raxml mafft r-base r-irkernel r-tidyverse r-plotly r-ggstar bioconductor-ggtree bioconductor-ggtreeextra r-phangorn bioconductor-tanggle (r-svglite)
- conda activate clustering
- jupyter lab

## Pipeline execution

Nearly fully automated as of 2022/01/19 in two Juyter Notebooks.

1. Clusterer.ipynb (Python): FASTA + Settings -> CSV tables
2. Graphics.ipynb (R):  CSV tables + Settings -> Grafics

## Problems

- Browser tab crashes unexpectedly when starting the jupyter server
  - fix: jupyter lab --notebook-dir=/home/user/Influenza/
  - cached content makes the notebook start in the wrong root folder
  - notebook-dir lets you reset the root location
