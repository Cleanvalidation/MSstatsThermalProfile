% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_cetsa_concentrations.R
\name{read_cetsa_concentrations}
\alias{read_cetsa_concentrations}
\title{Read in protein data from proteome discoverer}
\usage{
read_cetsa_concentrations(
  protein_path,
  peptide_path,
  Frac = TRUE,
  concentrations = set_concentration(16, concentrations =
    unique(CETSA_OnePot_annotation$Condition), replicates = 3),
  solvent = "DMSO",
  CARRIER = TRUE
)
}
\arguments{
\item{f.}{Path of Excel file output from Proteome Discoverer}
}
\value{
a dataframe containing extracted information
}
\description{
Read in Excel file and apply minimal pre-processing
}
