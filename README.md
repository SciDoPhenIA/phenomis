# **phenomis**: An R package for phenotyping data sciences

[![Travis build status](https://travis-ci.org/SciDoPhenIA/phenomis.svg?branch=master)](https://travis-ci.org/SciDoPhenIA/phenomis)

## Description

This package provides methods to perform the statistical analysis of phenomics datasets (e.g. in proteomics and metabolomics). These methods include the reading of datasets (as 3 table *dataMatrix*, *sampleMetadata* and *variableMetadata* .tsv files) into an **ExpressionSet** object (**reading**), quality control (**inspecting**) and transformation (**transforming**) of the dataMatrix, reduction of chemical redundancy (**reducing**), and univariate hypothesis testing (**hypotesting**). Multivariate analysis and feature selection can be further performed with the [**ropls**](http://bioconductor.org/packages/release/bioc/html/ropls.html) and [**biosigner**](http://bioconductor.org/packages/release/bioc/html/biosigner.html) packages, respectively. Finally, features can be annotated based on their m/z (and rt) values against public or local databases (**annotating**; based on the [**biodb**](https://github.com/pkrog/biodb) package). See the [**phenomis vignette**](vignettes/phenomis.Rmd) for a detailed example of the analysis of a metabolomics dataset.

## Contributors

Natacha Lenuzza, Alyssa Imbert, Marie Tremblay-Franco, Pierrick Roger, Eric Venot, Sylvain Dechaumet, and Etienne A. Thevenot.

## Maintainer

Etienne A. Thévenot ([etienne.thevenot@cea.fr](mailto:etienne.thevenot@cea.fr))   
Data Sciences for Molecular Phenotyping and Precision Medicine team (SciDoPhenIA)
Medicines and Healthcare Technologies Department  
CEA, INRAE, Paris Saclay University, MetaboHUB  
91191 Gif-sur-Yvette Cedex, France  
Web: [https://scidophenia.github.io](https://scidophenia.github.io/)  

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

* Imbert, A., Rompais, M., Selloum, M., Castelli, F., Mouton-Barbosa, E., Brandolini-Bunlon, M., Chu-Van, E., Joly, C., Hirschler, A., Roger, P., Burger, T., Leblanc, S., Sorg, T., Ouzia, S., Vandenbrouck, Y., Médigue, C., Junot, C., Ferro, M., Pujos-Guillot, E., de Peredo, A. G., Fenaille, F., Carapito, C., Herault, Y., & Thévenot, E. A. (2021). ProMetIS, deep phenotyping of mouse models by combined proteomics and metabolomics analysis. *Scientific Data*, 8(1). [doi:10.1038/s41597-021-01095-3](https://doi.org/10.1038/s41597-021-01095-3)

* Thévenot, E.A., Roux, A., Xu, Y., Ezan, E., and Junot, C. (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *Journal of Proteome Research* 14, 3322–3335. [doi:10.1021/acs.jproteome.5b00354](https://doi.org/10.1021/acs.jproteome.5b00354)

## Funding

* [CEA](http://www.cea.fr/english)
* [INRAE](https://www.inrae.fr/en)
* [MetaboHUB](https://www.metabohub.fr/home.html) [ANR-11-INBS-0010]

