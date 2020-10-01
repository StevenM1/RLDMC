# RLDMC

While all code is freely available, the current code is not written with a broad user group in mind. Hence, it isn't very flexible (in terms of the experimental paradigms it supports) and not everything is well documented. We're committed to including learning/adaptation effects in Dynamic Models of Choice in the future to make use this code-base more approachable.


Model fitting can be run from run.R, which requires a modelName (all options specified in models.R) and a dataName (exp1, exp2, exp3).
- doSample() samples, making heavy use of h.RUN.dmc();
- runDiagnostics() generates a directory with files containing info on sampling diagnostics (effective sample sizes, gelman's diagnostics, plots of traces, plots of posterior predictives / "fits");
- getPriors(), as defined in utils.R, contains all priors.

All code that generates the plots of the figures can be found in the makeFigure-*.R files.

Please note that Experiment 2 in the paper (SAT experiment) is experiment 3 (exp3, data/exp3.RData) here and vice versa.


--Steven Miletic
