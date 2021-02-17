# AcidGenomes

Toolkit for downloading and processing genome annotations.

## Installation

### [R][] method

```r
install.packages(
    pkgs = "AcidGenomes",
    repos = c(
        "https://r.acidgenomics.com",
        getOption("repos")
    )
)
```

### [Docker][] method

```sh
image="acidgenomics/r-acidgenomes"
workdir="/mnt/work"
docker pull "$image"
docker run -it \
    --volume="${PWD}:${workdir}" \
    --workdir="$workdir" \
    "$image" \
    R
```

[docker]: https://www.docker.com/
[r]: https://www.r-project.org/
