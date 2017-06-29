These files are functions, classes and routines used to produce the 
Benford and logistic analysis in the article "Bla bla bla" Arxiv:XXX

The file Eleicoes.py contains the classes, routines, functions and some
definitions used and the file benforder.py contains the script that performs the
analysis. The file texs/begin.tex contains the "header" for producing the 
table in latex format. The output of the script are a table showing, for each
party the number of donations received by each category (CNPJ, CPF, unknown,
non-original), followed by the fitted parameters (gamma, x0=exp(csi0)) for
each party/category (for sets with more than 20 donations). The script also
produces plots (histogram, accumulated distribution and first digit distribution
compared with Benford's) for each party/category. Finally, the script produces
a pdf table from the tex file it produces.

Note that this script is written to run in a linux system.
In order to use it in any other system, adaptations may
have to be made in the script (and in Eleicoes.py). Also, the following python
modules (libraries) are used: Gnuplot (gnuplot should also be installed in the 
system); ExGUtils (https://pypi.python.org/pypi/ExGUtils) from ExGUtils, the
functions to produce random number with homogeneous distribution is used (this
could be used from another modules) and the functions used to integrate,
to produce histograms and to evaluate averages and standard deviation
(these could also be adapted from other libraries, but
would cost more work to adapt); scipy.stats functions to evaluate the chi
square distribution are used.

The other analysis performed is the logistic regression whose main functions
are in the file Logit.py and the script performing the analysis in the file
logistic.py. The functions in Logit.py use the classes matrix and the function
diag from numpy and the same chi2 functions to deal with the chi2 distribution
that are used in the Eleicoes.py file.

Note that all scripts depend on directories and file names with the data
donwloaded from Brazil's superior electoral court website. We may have
downloaded these files in different directories and used different names and,
therefore, different users may have to adapt paths and filenames in order
to use our scripts. Variables with these names and paths are properly indicated
in the scripts.

