# **phenomis**: An R package for phenomics data analysis

[![Travis build status](https://travis-ci.org/SciDoPhenIA/phenomis.svg?branch=master)](https://travis-ci.org/SciDoPhenIA/phenomis)

## Description

This package provides methods to perform the statistical analysis of phenomics datasets (e.g. in proteomics and metabolomics). These methods include the reading of datasets (as 3 table *dataMatrix*, *sampleMetadata* and *variableMetadata* .tsv files) into an **ExpressionSet** object (**reading**), quality control (**inspecting**) and transformation (**transforming**) of the dataMatrix, and univariate hypothesis testing (**hypotesting**). Multivariate analysis and feature selection can be further performed with the [**ropls**](http://bioconductor.org/packages/release/bioc/html/ropls.html) and [**biosigner**](http://bioconductor.org/packages/release/bioc/html/biosigner.html) packages, respectively. Finally, features can be annotated based on their m/z (and rt) values against public or local databases (**annotating**; based on the [**biodb**](https://github.com/pkrog/biodb) package). See the [**phenomis vignette**](vignettes/phenomis.Rmd) for a detailed example of the analysis of a metabolomics dataset.

## Main contributors (SciDoPhenIA team)

Alyssa Imbert, Camille Roquencourt, Camilo Broc, Alexis Delabrière, Philippe Rinaudo, Krystyna Biletska, Natacha Lenuzza, Pierrick Roger and Etienne Thévenot.

## Maintainer

[Etienne A. Thevenot](https://etiennethevenot.pagesperso-orange.fr/)

## Methods

![](vignettes/figures/permanent/phenomis_methods.png)

## Formats

### 3 tabular file format used for import/export

Input (i.e. preprocessed) data consists of a 'samples times variables' matrix of intensities (**datMatrix** numeric matrix), in addition to sample and variable metadata (**sampleMetadata** and **variableMetadata** data frames). Theses 3 tables can be conveniently imported to/exported from R as tabular files:

![](vignettes/figures/permanent/phenomis_3table_format.png)

### **ExpressionSet** class used within the data analysis workflow

Within the R workflow, the **ExpressionSet** class perfectly handles these 3 tables (for additional information about *ExpressionSet* class, see the ['An introduction to Biobase and ExpressionSets'](https://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf) documentation from the [**Biobase**](https://doi.org/doi:10.18129/B9.bioc.Biobase) package).

## Installation

The package can be installed from github with `devtools::install_github("SciDoPhenIA/phenomis")`.

## Tutorial

See the [**phenomis vignette**](vignettes/phenomis.Rmd) for a detailed example of the analysis of a metabolomics dataset.

## Please cite

* Thévenot, E.A., Roux, A., Xu, Y., Ezan, E., and Junot, C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research* 14, 3322–3335. [doi:10.1021/acs.jproteome.5b00354](https://doi.org/10.1021/acs.jproteome.5b00354)

## Funding

* [CEA](http://www.cea.fr/english)
* [MetaboHUB](https://www.metabohub.fr/home.html) [ANR-11-INBS-0010]
* [IFB](https://www.france-bioinformatique.fr/en) [ANR-11-INBS-0013]
* [PhenoMeNal](http://phenomenal-h2020.eu/home/) [EINFRA H2020]

