# AcidGenomes

[![Install with Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/r-acidgenomes/README.html)

Toolkit for downloading and processing genome annotations.

## Installation

This is an [R][] package.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
install.packages(
    pkgs = "AcidGenomes",
    repos = c(
        "https://r.acidgenomics.com",
        BiocManager::repositories()
    ),
    dependencies = TRUE
)
```

### [Conda][] method

Configure [Conda][] to use the [Bioconda][] channels.

```sh
# Don't install recipe into base environment.
conda create --name='r-acidgenomes@0.4.8' 'r-acidgenomes==0.4.8'
conda activate 'r-acidgenomes@0.4.8'
R
```

[bioconda]: https://bioconda.github.io/
[conda]: https://docs.conda.io/en/latest/
[r]: https://www.r-project.org/
