# CoMEt #

CoMEt is a stochastic algorithm for identifying collections of mutually exclusive alterations in cohorts of sequenced tumor samples. CoMEt is written in Python, with required extensions written in C and Fortran. It was developed by the [Raphael research group](http://compbio.cs.brown.edu) in the [Department of Computer Science](http://cs.brown.edu) and [Center for Computational Molecular Biology](http://brown.edu/ccmb) at [Brown University](http://brown.edu).

CoMEt identifies a collection **M** of *t* alteration sets, each of size *k*, from a binary alteration matrix. CoMEt uses a Markov chain Monte Carlo (MCMC) algorithm to sample collections in proportion to their weight &phi;(**M**). The output of CoMEt is a list of collections, each with their sampling frequency, weight, and the weight &phi;(M) of each alteration set M &isin; **M**.

## Requirements ##

CoMEt requires the following Python modules. For each module, the latest version tested with CoMEt is given in parantheses:

1. [NetworkX](https://networkx.github.io/) (1.9.1)
2. [SciPy](http://www.scipy.org/) (0.14.1)
3. [NumPy](http://www.numpy.org/) (1.9.1).
4. [Multi-Dendrix](http://github.com/raphael-group/multi-dendrix) [optional].

CoMEt requires [Bower](http://bower.io/) to create web output.

## Setup ##

The C and Fortran extensions must be compiled before running CoMEt. To compile the extensions, run the following commands in your terminal:

    cd comet/
    python setup.py build

This will generate two compiled Python modules -- `comet/cComet.so` and `comet/permute_matrix.so` -- which can be imported directly into Python.

## Usage ##

### Input ###

The input data for CoMET consists of a:

1. **Alteration matrix**. This tab-separated file lists alterations in your dataset. Each row lists the alterations for a single sample. In each row, the first column lists a sample ID, and the remaining columns list genes that are altered in that sample. Note that the matrix is not necessarily symmetric, as different samples will have different numbers of alterations.
2. **Alteration whitelist** [optional]. If provided, the alteration matrix is restricted to only those alterations that also appear in this file.
3. **Patient whitelist** [optional]. If provided, the alteration matrix is restricted to only those samples that also appear in this file.

In all files, lines starting with `'#'` are ignored.

We provide example data in `example_datasets/`.

### Run CoMEt ###

CoMEt is run as a pipeline consisting of three steps:

1. **Run COMEt MCMC algorithm on real data**. Use the `run_comet.py` script to run the Markov chain Monte Carlo (MCMC) algorithm on the given mutation matrix. `run_comet.py` outputs a [JSON](http://json.org/) file that stores the parameters of the run, and a tab-separated file that lists the collections identified by CoMEt (sorted descending by sampling frequency).
2. **Run CoMEt MCMC algorithm on permuted data**. Use the `run_comet_permutation.py` script to find the maximum weight &phi;(**M**) across permuted datasets. This step is required for computing statistical significance and identifying the consensus modules. The output of this step is a JSON file that lists the parameters of the run, and the maximum &phi;(**M**) identified on the permuted data.
3. **Compute significance and create output web site**. Use the `compute_significance.py` script to compute the statistical significance of your results on real data. `compute_significance.py` outputs a website that can be used to visualize the results. To view the results website, download the required Javascript files (see Requirements above) and start a Python web server:

        cd OUTPUT_DIRECTORY # the output directory you provided to compute_significance.py
        bower install
        python -m SimpleHTTPServer 8000

  Then direct your browser to `http://localhost:8000`.

### Compute weights exhaustively ###

We also provide the script `run_exhaustive.py` as a simple way to compute the weight &phi;(M) for *all* gene sets *M* in a given dataset (using the same input format as above). The output of `run_exhaustive.py` is a tab-separated file that lists the weight &phi;(M) for all gene sets in the dataset (sorted ascending by &phi;(M)).

## Support ##

Please visit [our Google Group](https://groups.google.com/forum/#!forum/dendrix) to post questions and view discussions from other users, or contact us through [our research group's website](http://compbio.cs.brown.edu).

## Testing ##

To test CoMEt, run the following commands:

    cd test
    python test.py

The tests are successful if the last line of the text printed to the terminal is `"PASS"`.

## Reference ##

[Mark D.M. Leiserson](http://maxleiserson.com)\*, [Hsin-Ta Wu](http://cs.brown.edu/~bournewu/)\*, [Fabio Vandin](http://www.imada.sdu.dk/~vandinfa/), [Benjamin J. Raphael](http://compbio.cs.brown.edu). CoMEt: A Statistical Approach to Identify Combinations of Mutually Exclusive Alterations in Cancer. In *Proceedings of the 19th Annual Conference on Research in Computational Molecular Biology (RECOMB)* 2015. [Extended abstract](http://link.springer.com/chapter/10.1007%2F978-3-319-16706-0_19#page-1) and [preprint](http://arxiv.org/abs/1503.08224).

\* equal contribution
