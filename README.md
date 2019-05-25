Reproducibility package for the article:

**A practical introduction to the bootstrap: a versatile method to make inferences by using data-driven simulations**
Rousselet G.A., Pernet C.R., Wilcox R.R.
*in preparation*

OSF repository = ...
PsyArXiv Preprint = ...

The repository contains all of the [R](https://www.r-project.org/) code  used in the article. The code is best seen by running the RMarkdown notebooks, within [RStudio](https://www.rstudio.com/). 

# Content
|folder|description|
|-----|-----|
|`code`|R `.Rmd` files to run simulations and create figures|
|`notebooks`|pdf versions of the code, with embedded figures|
|`data`|simulation results needed to run the code|
|`figures`|all the figures used in the article, in pdf format|
|`functions`|extra R functions defined in text files|
|`docs`|html versions of the notebooks|
 
 # Notebooks

The notebooks contain code to reproduce the figures and analyses presented in the article. They also contain extra resources, figures and analyses.

|Notebook|Description|Figures|
|-----|-----|-----|
|[pb](/docs/pb.md)|Description of the percentile bootstrap|Figure 1 (figure_pb.pdf)|
|[pc](/docs/pc.md)|Percent correct example|Figure 2 (figure_pc.pdf)|
|[sampdist](/docs/sampdist.md)|Illustrate bootstrap sampling distributions|Figures 3-4 (figure_lognormal_pop.pdf, figure_sampdist_lognormal.pdf)|
|[coverage](/docs/coverage.md)|Simulations of the coverage, width and power of one-sample confidence intervals|Figures 5, 7, 8 (figure_onesamp_50ci.pdf, figure_onesamp_coverage.pdf, figure_onesamp_stability.pdf)|
|[notrobust](/docs/notrobust.md)|On its own, the bootstrap does not guarantee robustness|Figure 6 (figure_notrobust.pdf)|
|[2indgps](/docs/2indgps.md)|Compare 2 independent groups|Figure 9 (figure_2indgps.pdf)|
|[compcorr](/docs/compcorr.md)|Comparison of correlation coefficients|Figure 10 (figure_compcorr.pdf)|
|[2depgps](/docs/2depgps.md)|Illustrate hierarchical bootstrap sampling|Figure 11 (figure_hpb.pdf)|
|[ptb](/docs/ptb.md)|Percentile-t bootstrap technique|Figures 12-13 (figure_ptb.pdf, figure_ptb_explained.pdf)|
   
# R packages needed
If you want to run the code in RStudio, you will need to install a few packages. 

To reproduce the figures only, you can install the required packages by typing this in the console:

`install.packages(c("ggplot2", "tibble"))` 

Or you can navigate in the GUI to Tools > Install Packages...

To install `rogme` and `facetscales`, first you will need to install `devtools`:

`install.packages("devtools")`

then:

`devtools::install_github("GRousselet/rogme")`

`devtools::install_github("zeehio/facetscales")`

To reproduce the summary figures, you also need `cowplot` to combine panels:

`install.packages("cowplot")` 

Finally, if you decide to run the simulations, you will need `beepr` to get a little auditory reward:

`install.packages("beepr")` 
  
# Additional R functions 
Here we highlight a few R functions relevant to the tutorial. Most of them are listed in the RMarkdown notebooks, with example syntax. Each notebook will install the appropriate functions for you; otherwise, in the console you can type `source(file.choose())` and select the relevant .txt file. 

## Robust estimation and hypothesis testing
To get all the statistical functions from Rand Wilcox, select the [Rallfun-v35.txt](https://github.com/GRousselet/articles/blob/master/bootstrap/functions/Rallfun-v35.txt) file. See details on this [webpage](https://dornsife.usc.edu/labs/rwilcox/software/). The full description of the functions is available in the book [Introduction to Robust Estimation and Hypothesis Testing](https://books.google.co.uk/books/about/Introduction_to_Robust_Estimation_and_Hy.html?id=8f8nBb4__EYC&printsec=frontcover&source=kp_read_button&redir_esc=y#v=onepage&q&f=false). Here are some of the functions used or mentioned in the notebooks. 

### One-sample 
|Name|Description|
|-----|-----|
|`onesampb`|one-sample percentile bootstrap for any estimator|
|`sint`|parametric inference on the median|
|`trimci`|one-sample test on trimmed means|
|`trimpb`|percentile bootstrap inferences on trimmed means|
|`trimcibt`|bootstrap-t on trimmed means|
|`hdpb`|percentile bootstrap inferences on the Harrell-Davis quantile estimator|

### Two independent groups
|Name|Description|
|-----|-----|
|`yuen`|t-test on trimmed means|
|`yuenbt`|bootstrap-t on trimmed means|
|`pb2gen`|percentile bootstrap to compare any estimators|
|`medpb2`|same as pb2gen but only to compare medians|
|`trimpb2`|same as pb2gen but only to compare trimmed means|
|`comvar2`|parametric test to compare variances|

### Two dependent groups
|Name|Description|
|-----|-----|
|`yuend`|t-test on dependent trimmed means|
|`ydbt`|bootstrap-t on trimmed means|
|`comdvar`|parametric test of variances|
|`bootdpci`|percentile bootstrap using any estimator|

### Correlations
|Name|Description|
|-----|-----|
|`pcorb`|percentile bootstrap confidence interval for Pearson's correlation|
|`corb`|percentile bootstrap confidence interval for any robust correlation|
|`wincor`|winsorised correlation|
|`pbcor`|percentage bend correlation|
|`mscor`|skipped correlations using Pearson's or Spearman's correlation|

#### Compare independent correlations
|Name|Description|
|-----|-----|
|`twopcor`|percentile bootstrap comparison of two independent Pearson's correlations|
|`twocor`|percentile bootstrap comparison of two independent robust correlations|

#### Compare dependent correlations
|Name|Description|
|-----|-----|
|`TWOpov`|compare two overlapping Pearson's correlations|
|`twoDcorR`|compare two overlapping robust correlations|
|`TWOpNOV`|compare two non-overlapping Pearson's correlations|
|`twoDNOV`|compare two non-overlapping robust correlations|

## Other custom functions
[functions.txt](https://github.com/GRousselet/articles/blob/master/bootstrap/functions/functions.txt) and [theme_gar.txt](https://github.com/GRousselet/articles/blob/master/bootstrap/functions/theme_gar.txt) contain custom code to set some ggplot2 parameters and to compute a few things. Other custom functions are defined in the notebooks.

# Extra resources

## R packages for bootstrap inferences
- [`boot`](https://www.statmethods.net/advstats/bootstrapping.html) 
- [`resample`](https://cran.r-project.org/web/packages/resample/index.html)
- [`bootstrap`](https://cran.r-project.org/web/packages/bootstrap/index.html)
- [`WRS2`](https://cran.r-project.org/web/packages/WRS2/index.html)

## Interactive demo
[Frequentist inference: confidence interval & bootstrap](https://seeing-theory.brown.edu/frequentist-inference/index.html#section2)


## Books
Suggested books on bootstrap methods, robust statistics and simulations.

[An Introduction to the Bootstrap](https://books.google.co.uk/books?id=gLlpIUxRntoC&printsec=frontcover&dq=bootstrap+efron&hl=en&sa=X&ved=0ahUKEwjiv676orHiAhVIRxUIHYkgAckQ6AEIKDAA#v=onepage&q=bootstrap%20efron&f=false)

[Robust Statistics](https://books.google.co.uk/books?id=yAmZsxWSEWgC&dq=mathematical+foundations+of+robust+statistics&hl=en&sa=X&ved=0ahUKEwj3tP6Oo7HiAhU_UhUIHftBCnQQ6AEILjAB)

[Introduction to Robust Estimation and Hypothesis Testing](https://books.google.co.uk/books/about/Introduction_to_Robust_Estimation_and_Hy.html?id=8f8nBb4__EYC&printsec=frontcover&source=kp_read_button&redir_esc=y#v=onepage&q&f=false)

[Computer Age Statistical Inference](https://books.google.co.uk/books?id=Sj1yDAAAQBAJ&printsec=frontcover&dq=computer+inference+efron&hl=en&sa=X&ved=0ahUKEwiXy7vjorHiAhU4SxUIHUm7A3kQ6AEIKDAA#v=onepage&q=computer%20inference%20efron&f=false)

[Statistics: Unlocking the Power of Data](https://books.google.co.uk/books?id=EpBEDwAAQBAJ&printsec=frontcover&dq=Statistics:+Unlocking+the+Power+of+Data&hl=en&sa=X&ved=0ahUKEwil1unSorHiAhUPTxUIHfCwBjEQ6AEIKDAA#v=onepage&q=Statistics%3A%20Unlocking%20the%20Power%20of%20Data&f=false)

[Introduction to Statistical Investigations](https://books.google.co.uk/books?id=FsvVwQEACAAJ&dq=Introduction+to+Statistical+Investigations&hl=en&sa=X&ved=0ahUKEwigt_2worHiAhUySBUIHcZqCnAQ6AEIKDAA)

[Mathematical Statistics with Resampling and R](https://books.google.co.uk/books?id=_2hvDwAAQBAJ&printsec=frontcover&dq=Mathematical+Statistics+with+Resampling+and+R&hl=en&sa=X&ved=0ahUKEwiEj_-gorHiAhXkQhUIHRauC-IQ6AEILzAB#v=onepage&q=Mathematical%20Statistics%20with%20Resampling%20and%20R&f=false)
