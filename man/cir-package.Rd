\name{cir-package}
\alias{cir-package}
\alias{cir}
\docType{package}
\title{Isotonic Regression, Centered Isotonic Regression, and Dose-Response Utilities}
\description{This package revolves around centered isotonic regression (CIR), an improvement to isotonic regression (IR). However, it also includes a flexible, useful implementation of IR, confidence-interval estimates for both CIR and IR, and additional utilities for dose-response and dose-finding data.}

\details{

Isotonic regression (IR) is a standard nonparametric estimation method for monotone data. We have developed an improvement to univariate IR, named centered isotonic regression (CIR). There are heuristic and theoretical justifications to prefer CIR over IR, but first and foremost, in most simulations it produces substantially smaller estimation error. More details appear in Oron and Flournoy (2017).

This package implements CIR, but "along the way" an enhanced interface to univariate IR is also available. IR's base-R  implementation \code{isoreg} is very limited, as its own help page admits. A few other packages provide versions of IR, but to my knowledge the \code{cir} implementation offers some unique conveniences.

In addition, Oron and Flournoy (2017) also develop theoretically-backed confidence intervals applicable to both CIR and IR. The package's convenience wrapper \code{\link{quickIsotone}} executes CIR (or IR if one chooses \code{estfun = oldPAVA}), and returns both point and interval estimates at the specified x values. 

Since our motivation for studying IR comes from dose-finding designs such as Up-and-Down, there's analogous functionality for dose-finding ("inverse") estimation of x given specified y values. In particular, \code{\link{quickInverse}} offers inverse point and interval estimates in a single call. The package now also includes an optional bias-correction shrinkage method for such designs, informed by more recent research (Flournoy and Oron, 2020).

The package's focus is dose-response data with the response assumed binary (coded as 0 or 1). Some functions might work for any input data, but others will not. In particular, the confidence intervals are only applicable to binary-response data.

The package also includes two S3 classes, \code{\link{doseResponse}} and \code{\link{DRtrace}}. The former which is more heavily used, is a data frame with elements \code{x, y, wt}, summarizing the dose-response information. The latter is a "trace" or a running description of raw dose-response data, with \code{x, y, cohort} provided at the resolution of single observations. Each class has a \code{plot} method.

If you intend to use \code{cir} mostly for analysis of an Up-and-Down experiment, note that the newer package \code{upndown} contains more convenient wrapper utilities for such usage. These wrappers use \code{cir} functions to carry out the basic estimation and visualization tasks.

Enjoy!
}
\author{
Assaf P. Oron.

Maintainer: Assaf P. Oron \code{<assaf.oron.at.gmail.com>}
}

\references{Oron, A.P. and Flournoy, N., 2017. Centered Isotonic Regression: Point and Interval Estimation for Dose-Response Studies. Statistics in Biopharmaceutical Research 9, 258-267. (author's public version available on arxiv.org).

Flournoy, N. and Oron, A.P., 2020. Bias Induced by Adaptive Dose-Finding Designs. Journal of Applied Statistics 47, 2431-2442.}


