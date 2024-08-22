# R2easyviz is a package intended to provide an easy and quick functions for plotting single cell data in Seurat.

#### The movation for this package comes from types of plots that I create everyday while working with single cell projects. Rather than going through the work to manipulate the data of interest each time, R2easyviz provides an easy way to access data of interest, organize it in a meaningful way, and plot the results, all from the basic Seurat structure.

## Installation
#### Install Via GitHub

```r
install.packages("devtools")
devtools::install_github("ryanreis333/R2easyviz")
```

![Logo](images/HGNC_image.png)

#### The package relies on the ggplot2, dplyr, stringr, and Seurat packages.

#### If there is any issues with loading these packages alongside R2easyviz, you can download them all by running the following:
```r
install.packages(c("ggplot2", "dplyr", "stringr", "Seurat")
```

#### Alternatively, you can download them each individually with the following:
```r
install.packages("ggplot2")
install.packages("dplyr")
install.packages("stringr")
install.packages("Seurat")
```

## Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/ryanreis333/R2easyviz/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

#### [Pull Requests](https://github.com/ryanreis333/R2easyviz/pulls) are welcome for bug fixes, new features, or enhancements.
