# Genome Annotator

This is a simple example package. You can use
[Github-flavored Markdown](https://guides.github.com/features/mastering-markdown/)
to write your content.

## Start Here :-)

So installation here is a little robust due to the nature of how many programs are used. In order for everything to work smoothly, there are 3 main considerations:  
<code style="color : red">First</code> is the virtual environment you will be in while executing the pipeline. You want to make sure you are in a virtual environment so all the dependencies can play nice.  
<code style="color : red">Second</code> is the programs you install. You will need access to a lot of different programs, which should be housed within your virtual environment.  
<code style="color : red">Third</code> is the databases needed to annotate your data. Some programs have built in data permutations, but many require outside databases in order to compare your data to external data.

The next thing to consider is how to run the program. There are <code style="color : green">3 snakemake</code> files: Download.smk, Filter.smk, and Annotate.smk. Each can be assess with the command: <code>\$ genotator [download | filter | annotate] ...params</code>. So you choose either the download, filter, or annotate command. Note that everything is chained together such as <code>download --> filter --> annotate</code>, so if you run filter, you're also running download. If you run annotate, you're also running download and filter. See below for detailed instructions on how to run the program.

## Installation

<code style="color : red">Part 1: Create an isolated virtual environment using conda.</code> See [conda documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for more details on conda.

```bash
conda create -y --name gen-annotator python=3.10
source activate /path/to/environment/gen-annotator
conda install -y -c conda-forge mamba
```

### Download all of the dependencies

<code style="color : red">Part 2: Install all the dependencies</code>. Sorry there are so many :anguished:

Start by downloading mamba, which is how I prefer to install packages when I can (over conda). Almost always mamba is a lot faster at resolving dependencies and the whole installation process. See [mamba documentation](https://github.com/mamba-org/mamba) for more information.

```bash
conda install -y -c conda-forge mamba

mamba install -y -c conda-forge -c bioconda python=3.10 \
        sqlite prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript
mamba install -y -c bioconda fastani
conda install -c bioconda dbcan
pip install checkm-genome

# Grab the latest anvio version
curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz --output anvio-8.tar.gz
pip install anvio-8.tar.gz
# pip install fastani (deprecated)
```

Below is where we setup some databases through the Anvi'o suite.

```bash
- anvi-setup-pfams
  https://github.com/merenlab/anvio/blob/master/anvio/pfam.py#L55
- OR, I can download it directly for them into the lib

- anvi-setup-ncbi-cogs (Run_COGs)

- anvi-setup-kegg-kofams (Run_Kegg)
```

### Download NCBI's datasets

```bash
conda install -y -c conda-forge ncbi-datasets-cli
```

### Now for some miscellaneous installs

<code style="color : red">Note that RAST may be finnicky</code> :bowtie:

- pip install rast

##<span style="color:red">Common Issues</span>

#### CazYME:

If the genome is too small and nothing was annotated, it gives an error!

- Probably have some sort of:
  $ genome-annotator --setup-dbs  
  Line

### Workflow for users:

1. Install the conda package ($ conda install -c bioconda genome-annotator)
2. Run the database configuration script ($ annotator --init-dbs {all,cazyme,cogs,kegg,pfam,rast,tigr

<hr style="border:5px solid">

#<span style="color:green">User Guide</span>
Note: Any command starting with <code style="color : blue">\$</code> is used to denote being in a terminal (on the command line), and the <code style="color : blue">\$</code> marks the prompt.

### Part 1. Setting up the directories

1. Clone the repository: <code style="color : blue">$ git clonehttps://github.com/ddeemerpurdue/genome-annotator.git</code>
2. Enter the new directory: <code style="col`or : blue">$ cd genome-annotator</code>
3. Create the fresh conda environment: <code style="color : green">SEE INSTALLATION</code>
4. Download and initialize databases: <code style="color : green">SEE INSTALLATION</code> [Link to Installation](#installation)
5. <code style="color : red">Add genotator.py to path</code>: export \$PATH=\$PATH:/path/to/genome-annotator/bin/
6. Pass
