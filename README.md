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
```

### Gothe et al toolkit

To install this software simply type:

```
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python get-pip.py
pip install --user XXXXXXXXX
```

## Usage Example

In this section we will follow step by step the logic of this toolkit to
create the Figure 4C of the paper (and its corresponding numeric table; which can be easily
opened in Excel to make further plots). Depending on the versions of the dependencies installed,
minor numerical fluctuations are expected, since dynamic programming was used at some stages to
accelerate the execution. However, this will not modify the general trend of such a Figure.

The tool associated with this Figure is called **corr-dsb-dist-exp** . By calling from the command-line:

```
corr-dsb-dist-exp --help
```

You will see the list of arguments that this tool requires:

- max-dist: 
- alias-file: 
- chrom-sizes: 
- genes-file: 
- expression-file: 
- dsb-file: 
- distance-file: 
- temp-loc: 
- output-location: 



corr-dsb-dist-exp --max-dist 50 --alias-file /Users/egg/Projects/Test/Input/alias_hg19.txt --chrom-sizes /Users/egg/Projects/Test/Input/chrom.sizes.hg19.txt --genes-file /Users/egg/Projects/Test/Input/all_genes.bed --expression-file /Users/egg/Projects/Test/Input/K562_expression.txt --dsb-file /Users/egg/Projects/Test/Input/K562_ETO.bam --distance-file /Users/egg/Projects/Test/Input/K562_loops.txt --temp-loc /Users/egg/Projects/Test/Temp/ --output-location /Users/egg/Projects/Test/Output/

