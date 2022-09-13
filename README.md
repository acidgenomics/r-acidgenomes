# AcidGenomes

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
name='r-acidgenomes'
conda create --name="$name" "$name"
conda activate "$name"
R
```

### [Docker][] method

```sh
image='acidgenomics/r-packages:acidgenomes'
workdir='/mnt/work'
docker pull "$image"
docker run -it \
    --volume="${PWD}:${workdir}" \
    --workdir="$workdir" \
    "$image" \
    R
```

[bioconda]: https://bioconda.github.io/
[conda]: https://docs.conda.io/en/latest/
[docker]: https://www.docker.com/
[r]: https://www.r-project.org/
