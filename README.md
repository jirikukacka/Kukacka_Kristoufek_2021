# Code repository for Kukacka and Kristoufek (2021)

---

Kukacka, J., Kristoufek, L. (2021), Does parameterization affect the complexity of agent-based models?, SSRN Working Paper, [DOI](https://dx.doi.org/10.2139/ssrn.3654362), version July 30, 2021.

This codebase is using the R language, and it is initialized as a Git.

It is authored by [Jiri Kukacka](https://ies.fsv.cuni.cz/en/staff/kukacka). Question can be directed to [jiri.kukacka@fsv.cuni.cz](jiri.kukacka@fsv.cuni.cz).

---

### Repository structure

* `simulation_gh07_cd.R`         master file #1. It provides R code for an illustrative replication of the results for the [Gaunersdorfer and Hommes (2007)](https://doi.org/10.1007/978-3-540-34625-8_9) model: subsection 5.2.2. Computational setup is trivialized for personal computers. The original setup for a cluster computation to replicate results in the paper is also included. Run this file first.

* `heatmaps_gh07_cd.R`           master file #2. It produces the presented heat map graphics: Fig. 4, panels (c) and (d) and Fig. B.11, panel (d). It loads the output of `simulation_gh07_cd.R`.

* `src`                            includes additional source scripts:

   * `gh07_cd.R`                    function generating the model output. This part of the code was downloaded in MATLAB language and translated to R from the [agentFin 0.1 documentation web page](\href{http://people.brandeis.edu/~blebaron/classes/agentfin/GaunersdorferHommes.html) [accessed 2019-10-28] by Prof. Blake LeBaron (Brandeis University, USA). Only inevitable modifications were made to link the model script to the multifractality analysis. Credit: Patrick Herb, Blake LeBaron, Axel Szmulewiez.

   * `gh07_cd_MF.R`                 function supporting evaluation of independent for-loop iterations in parallel to estimate multifractality measures.

   * `MF_DFA.R`                     set of functions needed for estimation of the multifractal spectrum, finally returning $\Delta H$ and $\Delta \alpha$ of estimated generalized Hurst exponents $H(q)$ for a given range of $q$s based on MF-DFA.

   * `heatmap.2`                    support function for `heatmaps_gh07_cd.R`, downloaded from [https://rdrr.io/cran/gplots/src/R/heatmap.2.R](rdrr.io/cran/gplots/src/R/heatmap.2.R) [accessed 2019-06-24].  

* `replication_of_results_gh07_cd` includes the files for replication of Fig. 4, panels (c) and (d) and Fig. B.11, panel (d):

   * `gh07_cd.Rda`                  data set for replication of the results in subsection 5.2.2. It must be copy-pasted to the main root to be loaded by `heatmaps_gh07_cd.R` to produce the following graphical outputs (for verification purposes):

      * `heatmap_gh07_c_mfst.pdf`          actual Fig. 4, panel (c).
   
      * `heatmap_gh07_d_conf.pdf`          actual Fig. 4, panel (d).

      * `heatmap_gh07_diff.pdf`            actual Fig. B.11, panel (d).

---

### Computational setup

* multifractal detrended fluctuation analysis follows the setting in subsection 2.2.

* general Monte Carlo setup follows subsection 3.1.

* `simulation_gh07_cd` further follows the simulation setup in subsection 3.2. and the model benchmark setup in subsection 4.2.2.

* `heatmaps_gh07_cd.R` follows the multifractality evaluation logic in subsection 3.3.

---

### Required packages

* install all by `Pkg.add(["doParallel", "foreach", "iterators", "paralell", "gplots", "gtools"])`
   * `doParallel`
   * `foreach`
   * `iterators`
   * `paralell`
   * `gplots`
   * `gtools`

---

### Variables and function arguments

* general simulation setup:
   * `runs`    number of independent Monte Carlo runs for each experiment.
   * `burn`    number of burn-in periods to discard.
   * `obs`     number of periods desired.
   * `metrics` actual number of estimates generated in `gh07_cd_MF.R`.
* parameterization:
   * `k_val`   benchmark parameterization for $\psi$, associated with the vertical axis of heatmaps.
   * `j_val`   benchmark parameterization for $\gamma$, associated with the horizontal axis of heatmaps.
   * `kk,kkk`  vertical position in the grid.
   * `jj,jjj`  horizontal position in the grid.
* multifractal spectrum estimation setup:
   * `smin`    minimum scale to consider for the scaling law.
   * `smax`    maximum scale to consider for the scaling law.
   * `sstep`   step length for scales between smin and smax.
   * `qmax`    maximum moment $q$ for estimation of the multifractal spectrum.
   * `qmin`    minimum moment $q$ for estimation of the multifractal spectrum.
   * `qstep`   step length for moments between qmin and qmax.