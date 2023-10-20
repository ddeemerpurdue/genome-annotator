# Genome Annotator

This is a simple example package. You can use
[Github-flavored Markdown](https://guides.github.com/features/mastering-markdown/)
to write your content.

# Installation

```bash
conda create -y --name gen-annotator python=3.10
source activate /path/to/environment/gen-annotator
conda install -y -c conda-forge mamba
# Download all of the dependencies
mamba install -y -c conda-forge -c bioconda python=3.10 \
        sqlite prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript
mamba install -y -c bioconda fastani
# Grab the latest anvio version
curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz --output anvio-8.tar.gz
pip install anvio-8.tar.gz
pip install fastani
```

```bash
- anvi-setup-pfams
  https://github.com/merenlab/anvio/blob/master/anvio/pfam.py#L55
- OR, I can download it directly for them into the lib

- anvi-setup-ncbi-cogs (Run_COGs)

- anvi-setup-kegg-kofams (Run_Kegg)
```

## Download NCBI's datasets

```bash
conda install -y -c conda-forge ncbi-datasets-cli
```

- pip install rast

#### Issue with CAZYME:

If the genome is too small and nothing was annotated, it gives an error!

- Probably have some sort of:
  $ genome-annotator --setup-dbs  
  Line

### Workflow for users:

1. Install the conda package ($ conda install -c bioconda genome-annotator)
2. Run the database configuration script ($ annotator --init-dbs {all,cazyme,cogs,kegg,pfam,rast,tigr
