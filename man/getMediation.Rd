% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/mediate.R
\name{addMediation}
\alias{addMediation}
\alias{getMediation}
\alias{addMediation,SummarizedExperiment-method}
\alias{getMediation,SummarizedExperiment-method}
\title{Perform mediation analysis}
\usage{
addMediation(x, ...)

getMediation(x, ...)

\S4method{addMediation}{SummarizedExperiment}(
  x,
  outcome,
  treatment,
  name = "mediation",
  mediator = NULL,
  assay.type = NULL,
  dimred = NULL,
  family = gaussian(),
  covariates = NULL,
  p.adj.method = "holm",
  add.metadata = TRUE,
  verbose = TRUE,
  ...
)

\S4method{getMediation}{SummarizedExperiment}(
  x,
  outcome,
  treatment,
  mediator = NULL,
  assay.type = NULL,
  dimred = NULL,
  family = gaussian(),
  covariates = NULL,
  p.adj.method = "holm",
  add.metadata = TRUE,
  sort = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{a
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.}

\item{...}{additional parameters that can be passed to
\code{\link[mediation:mediate]{mediate}}.}

\item{outcome}{\code{Character scalar}. Indicates the colData variable used
as outcome in the model.}

\item{treatment}{\code{Character scalar}. Indicates the colData variable
used as treatment in the model.}

\item{name}{\code{Character scalar}. A name for the column of the
\code{colData} where results will be stored. (Default: \code{"mediation"})}

\item{mediator}{\code{Character scalar}. Indicates the colData variable used
as mediator in the model. (Default: \code{NULL})}

\item{assay.type}{\code{Character scalar}. Specifies the assay used for
feature-wise mediation analysis. (Default: \code{NULL})}

\item{dimred}{\code{Character scalar}. Indicates the reduced dimension
result in \code{reducedDims(object)} for component-wise mediation analysis.
(Default: \code{NULL})}

\item{family}{\code{Character scalar}. A specification for the outcome model
link function. (Default: \code{gaussian("identity")})}

\item{covariates}{\code{Character scalar} or \code{character vector}.
Indicates the colData variables used as covariates in the model.
(Default: \code{NULL})}

\item{p.adj.method}{\code{Character scalar}. Selects adjustment method
of p-values. Passed to \code{p.adjust} function.
(Default: \code{"holm"})}

\item{add.metadata}{\code{Logical scalar}. Should the model metadata be
returned. (Default: \code{TRUE})}

\item{verbose}{\code{Logical scalar}. Should execution messages be printed.
(Default: \code{TRUE})}

\item{sort}{\code{Logical scalar}. Should the results be sorted by decreasing
significance in terms of ACME_pval. (Default: \code{FALSE})}
}
\value{
\code{getMediation} returns a \code{data.frame} of adjusted p-values and
effect sizes for the ACMEs and ADEs of every mediator given as input, whereas
\code{addMediation} returns an updated
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
instance with the same \code{data.frame} stored in the metadata as
"mediation" or as specified in the \code{name} argument. Its columns include:

\describe{
\item{mediator}{the mediator variable}
\item{acme}{the Average Causal Mediation Effect (ACME) estimate}
\item{acme_pval}{the original p-value for the ACME estimate}
\item{acme_lower}{the lower bound of the CI for the ACME estimate}
\item{acme_upper}{the upper bound of the CI for the ACME estimate}
\item{ade}{the Average Direct Effect (ADE) estimate}
\item{ade_pval}{the original p-value for the ADE estimate}
\item{ade_lower}{the lower bound of the CI for the ADE estimate}
\item{ade_upper}{the upper bound of the CI for the ADE estimate}
\item{total}{the Total Effect estimate}
\item{total_lower}{the lower bound of the CI for the Total Effect estimate}
\item{total_upper}{the upper bound of the CI for the Total Effect estimate}
\item{total_pval}{the original p-value for the Total Effect estimate}
\item{acme_padj}{the adjusted p-value for the ACME estimate}
\item{ade_padj}{the adjusted p-value for the ADE estimate}
\item{total_padj}{the adjusted p-value for the Total Effect estimate}
}

The original output of \code{\link[mediation:mediate]{mediate}} for each
analysed mediator is stored as the "model_metadata" attribute of the
resulting \code{data.frame} and can be accessed via the \code{attr} function.
}
\description{
\code{getMediation} and \code{addMediation} provide a wrapper of
\code{\link[mediation:mediate]{mediate}} for
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}.
}
\details{
This wrapper of \code{\link[mediation:mediate]{mediate}} for
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
provides a simple method to analyse the effect of a treatment variable on an
outcome variable found in \code{colData(se)} through the mediation of either
another variable in colData (argument \code{mediator}) or an assay
(argument \code{assay.type}) or a reducedDim (argument \code{dimred}).
Importantly, those three arguments are mutually exclusive.
}
\examples{
\dontrun{
# Import libraries
library(mia)
library(miaViz)
library(scater)

# Load dataset
data(hitchip1006, package = "miaTime")
tse <- hitchip1006

# Agglomerate features by family (merely to speed up execution)
tse <- agglomerateByRank(tse, rank = "Phylum")
# Convert BMI variable to numeric
tse$bmi_group <- as.numeric(tse$bmi_group)

# Analyse mediated effect of nationality on BMI via alpha diversity
# 100 permutations were done to speed up execution, but ~1000 are recommended
med_df <- getMediation(
    tse,
    outcome = "bmi_group",
    treatment = "nationality",
    mediator = "diversity",
    covariates = c("sex", "age"),
    treat.value = "Scandinavia",
    control.value = "CentralEurope",
    boot = TRUE, sims = 100)

# Visualise model statistics
plotMediation(med_df)

# Apply clr transformation to counts assay
tse <- transformAssay(tse, method = "clr", pseudocount = 1)

# Analyse mediated effect of nationality on BMI via clr-transformed features
# 100 permutations were done to speed up execution, but ~1000 are recommended
tse <- addMediation(
    tse, name = "assay_mediation",
    outcome = "bmi_group",
    treatment = "nationality",
    assay.type = "clr",
    covariates = c("sex", "age"),
    treat.value = "Scandinavia",
    control.value = "CentralEurope",
    boot = TRUE, sims = 100,
    p.adj.method = "fdr")

# Show results for first 5 mediators
head(metadata(tse)$assay_mediation, 5)

# Perform ordination
tse <- addMDS(
    tse, name = "MDS", method = "euclidean", assay.type = "clr",
    ncomponents = 3)

# Analyse mediated effect of nationality on BMI via NMDS components
# 100 permutations were done to speed up execution, but ~1000 are recommended
tse <- addMediation(
    tse, name = "reddim_mediation",
    outcome = "bmi_group",
    treatment = "nationality",
    dimred = "MDS",
    covariates = c("sex", "age"),
    treat.value = "Scandinavia",
    control.value = "CentralEurope",
    boot = TRUE, sims = 100,
    p.adj.method = "fdr")

# Show results for first 5 mediators
head(metadata(tse)$reddim_mediation, 5)

# Access model metadata
attr(metadata(tse)$reddim_mediation, "model_metadata")
}

}
