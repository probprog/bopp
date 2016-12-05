# BOPP: Bayesian Optimization for Probabilistic Programs

- Last Bitbucket commit before this repo: [43930ea4cb08026b1d050a496262e8a6e923ce85](https://bitbucket.org/tuananhle/bopp/commits/43930ea4cb08026b1d050a496262e8a6e923ce85).
- Depedencies: [Deodorant](http://github.com/twgr/deodorant) and a temporary snapshot [Anglican](https://clojars.org/org.clojars.tuananhle/anglican) (this will be switched to generic depedency on Anglican once the corresponding pull request is processed)
- This is still a working version and is different to the code used for the paper.  We have made so people can use it during the conference, but intend to keep working on it through the conference and after.


BOPP is a package for automated marginal maximum a posteriori inference (MMAP) based around the
probabilistic programming system [Anglican](http://www.robots.ox.ac.uk/~fwood/anglican).  The
user only needs to write their model in the same manner as existing Anglican programs and by using
the `defopt` construct instead of `defquery`, select the variables to be optimized, with the
rest marginalized out.  It can also be used as a means of exploiting the target source code
to improve Bayesian optimization, delivering things such as automatic domain scaling,
unbounded optimization, and implicit constraint satisfaction including equality constraints.
The key idea is to use a series of code transformations to extract from the original program
all the things that are needed to carry out the MMAP problem, such as the target function itself
and a program for optimizing the acquisition function subject to the implicit constraints.  These
are then passed to our other package [Deodorant](http://github.com/twgr/deodorant), which uses
these to solve the problem probabilistic programs creates for BO.  The following paper should
be referred to for full algorithmic details and we ask that you cite this paper if you
use BOPP in your work.

Rainforth, T., Le, T. A., van de Meent, J.-W., Osborne, M. A., & Wood, F. (2016). Bayesian Optimization for Probabilistic Programs. In Advances in Neural Information Processing Systems.

```
@incollection{rainforth2016bayesian,
    title = {Bayesian Optimization for Probabilistic Programs},
    author = {Rainforth, Tom and Le, Tuan Anh and van de Meent, Jan-Willem and Osborne, Michael A and Wood, Frank},
    booktitle = {Advances in Neural Information Processing Systems 29},
    editor = {D. D. Lee and M. Sugiyama and U. V. Luxburg and I. Guyon and R. Garnett},
    pages = {280--288},
    year = {2016},
    publisher = {Curran Associates, Inc.},
    url = {http://papers.nips.cc/paper/6421-bayesian-optimization-for-probabilistic-programs.pdf}
}
```

## Usage ##

### Defining target programs ###

BOPP targets are specified using the macro defopt.  This is identical to defquery in [Anglican](http://www.robots.ox.ac.uk/~fwood/anglican) except that it takes as an extra input of the target variables to be optimized.  
BOPP will optimize the evidence of the program with respect to the variables, whilst marginalizing out over the others.
For example:

```lisp
(defopt simple-bimodal [y] [theta]
  (let [x (sample (normal 0 1))
        theta (sample (normal x 0.5))]
    (observe (normal (sqrt (* theta theta)) 0.5) y)))
```

specifies a model where we wish to optimize theta, marginalizing out x.  

There are a couple of small restrictions
on the programs that can be specified to ensure that they constitute valid target programs, namely:

1. Each target variable must be bound to a value directly by a sample statement with fixed measure-type distribution argument (i.e. not a weird `defdist` distribution object).

2. The program must be written such that any possible execution trace binds each optimization variable exactly once.

3. Although any target variable may be lexically multiply bound, it must have the same base measure in all possible execution traces.

Don't worry, if this doesn't make too much sense, BOPP catches violations of these automatically and gives you
and error telling you what you violated. Please see the [paper](http://papers.nips.cc/paper/6421-bayesian-optimization-for-probabilistic-programs) and [supplementary material](http://www.robots.ox.ac.uk/~twgr/assets/pdf/rainforth2016BOPP.pdf) for more information.

### Running BOPP ###

Calling BOPP is super simple, just call

`(doopt algorithm opt-query opt-query-args num-samples args)`

where
- `algorithm` = Inference algorithm used for estimating the marginal [:smc, :pcascade, :importance]
- `opt-query` = Query to be optimized, defined by defopt.
- `opt-query-args` = Fixed inputs of opt-query-args (y in earlier example)
- `num-samples` = Number of samples for constructing the estimator of the marginal.
- `args` = Optional arguments as key-value pairs.  See docstring for doopt.

This returns a lazy infinite sequence of samples.  Call (take N .) on this will convert it
to a fixed number of outputs.  Note that as BOPP is a GP based BO scheme, scaling in N is quite poor
and practically limited to around 400-500.

### Examples ###

A number of worksheets are provided to give example usage in different cases.  These can be accessed in a Gorilla REPL by running
`lein gorilla` from the base folder and opening worksheets in the `worksheets/` folder.

## Installation ##

Though BOPP currently runs of a snapshot of Anglican that means you don't need to install Anglican
explicitly, is does have the same requirements in terms of java, Leiningen etc  and so we refer the reader to http://www.robots.ox.ac.uk/~fwood/anglican/usage/index.html and recommend that users follow section 2 in the user start up guide.  We also recommend that you familiarize yourself on the syntax of Anglican through the
provided link, as this is the same syntax as Anglican, with the addition of the forms `defopt` and `doopt`.

## License ##

Copyright © Tom Rainforth, Tuan Anh Le, Jan-Willem van de Meent

Deodorant is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Deodorant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the [GNU General Public License](gpl-3.0.txt) along with Deodorant.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
