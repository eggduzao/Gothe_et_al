# **Gothe et al. Spatial chromosome folding and active transcription drive DNA fragility and formation of oncogenic MLL translocations.**

This repository contains the *in house* computational scripts in python and R that were used in the study by Gothe et al.:
Spatial chromosome folding and active transcription drive DNA fragility and formation of oncogenic MLL translocations.

The scripts were turned into command-line tools and can be executed in a Unix and MAC OS environments.
This Mini-tutorial provides a Quick-Start to use the available command-line tools.

## Installation

The first step is to have Python and R installed. Most Unix and MAC distributions already have this
programming languagess installed. If you do not have them installed please visit:

- https://www.python.org/

- https://www.r-project.org/

### Dependencies

When installing this software, all its dependencies will be automatically installed.
If you encounter any problem. Please find bellow you will find the quickest way to install all the dependencies:

```
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py
python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
pip install --user pysam
pip install --user pyBigWig
pip install deeptools
sudo apt-get install unzip
curl https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xjf samtools-1.9.tar.bz2
cd samtools-1.9
./configure --prefix=/where/to/install
make
make install
```

The program also requires a couple of R packages. After a successful installation of R please type:

```
R
install.packages(c("MASS", "OneR", "RColorBrewer", "colorspace", "ggplot2", "ggthemes", "gplots", "lattice", "plot3D", "plotly", "plotrix", "reshape", "scatterplot3d"))
```

### Gothe et al toolkit

To install this software simply type:

```
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py
cd
git clone https://github.com/eggduzao/Gothe_et_al.git
cd Gothe_et_al
pip install --user .
```

## Usage Example

In this section we will follow step by step the logic of this toolkit to
create the Figure 4C of the paper (and its corresponding numeric table; which can be easily
opened in Excel to make further plots).

This Figure shows the triple (and pairwise) correlation between double strand breaks (DSBs),
expression and the gene's distance from its closest loop anchor.

Depending on the versions of the dependencies installed,
minor numerical fluctuations are expected, since dynamic programming was used at some stages to
accelerate the execution. However, this will not modify the general trend of such a Figure.

The tool associated with this Figure is called **corr-dsb-dist-exp** . By calling from the command-line:

```
corr-dsb-dist-exp --help
```

You will see the list of arguments that this tool requires:

- max-dist: The maximum distance a gene can be from its nearest loop anchor. You should set this to an integer twice as much as you expect in the final plot.
- alias-file: The alias file is given with the package and simply contains multiple collected aliases for the same gene.
- chrom-sizes: This is a tab-separated file containing the chromosome names and their lengths.
- genes-file: A bed file containing all genes in which the analysis will be based on.
- expression-file: It can be a tab-separated file containing the gene names and their expression; or a bam file in which the expression will be calculated.
- dsb-file: A bed or a bam file containing the double strand breaks' positions.
- distance-file: At the moment, the tool only supports loops called from HiCCUPS. CTCF annotation is not required but highly recommended as in Rao et al.'s file in the example bellow.
- temp-loc: A path in which the program will store all temporary files. It can be erased after the execution; however the program itself won't erase this path as it might be useful for troubleshooting.
- output-location: The output location in which the tables and figures will be created in.

To run the tool you can simply create *input*, *output* and *temporary* folders and
copy or download the relevant data to the input folder.

We have already created some initial files for this script. The intention is to
make the execution faster. However, you are free to try with your own data!

Supposing gothe_dir is the directory of the repository

```
mkdir -p ~/input/ ~/output/ ~/temp/
cd gothe_dir
ln -s gothe_dir/data/genome/alias_hg19.txt ~/input/alias_hg19.txt
ln -s gothe_dir/data/genome/all_genes.bed ~/input/all_genes.bed
ln -s gothe_dir/data/genome/chrom.sizes.hg19.txt ~/input/chrom.sizes.hg19.txt
ln -s gothe_dir/data/expression/K562_expression.txt ~/input/K562_expression.txt
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_K562_HiCCUPS_looplist_with_motifs.txt.gz -o ~/input/K562_loops.txt.gz
```

Also, please log in the password-protected GEO repository and place the sBLISS file named *GSM3444989_K562_ETO.bed.gz* in ~/input/

Finally, you will be able to create the Figures by applying the following command:

```
corr-dsb-dist-exp --max-dist 200 --alias-file ~/input/alias_hg19.txtt --chrom-sizes ~/input/chrom.sizes.hg19.txt --genes-file ~/input/all_genes.bed --expression-file ~/input/K562_expression.txt --dsb-file ~/input/GSM3444989_K562_ETO.bed.gz --distance-file ~/input/K562_loops.txt.gz --temp-loc ~/temp/ --output-location ~/output/
```

