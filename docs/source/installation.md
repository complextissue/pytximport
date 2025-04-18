# Installation

:::{card}
:class-card: sd-bg-warning
:class-body: sd-bg-text-warning
**pytximport** only supports Python versions greater than or equal to **3.10**.
:::

:::{card} Recommendation
Installation via `Bioconda` is the recommended way to include `pytximport` in your projects. A `pip`-installabe package
is also provided.
:::

:::{card} Performance
While not required, we recommend users also install `pyarrow` for faster import of tab-separated value-based quantification files.
:::

## Installation Options

Choose an option to install this package.

::::{tab-set}

:::{tab-item} Bioconda
Install `pytximport` from `Bioconda` using `mamba` or `conda`:

```bash
mamba install -c bioconda pytximport
mamba install -c conda-forge pyarrow-core
```

:::

:::{tab-item} PyPi
Install `pytximport` package using `pip`:

```bash
python3 -m pip install pytximport pyarrow
```

:::

:::{tab-item} GitHub
Install `pytximport` from GitHub using `pip`:

```bash
python3 -m pip install git+https://github.com/complextissue/pytximport.git
```

:::

:::{tab-item} Source
This option is only recommended for potential contributors and installs additional developement dependencies.
Install `pytximport` from source:

```bash
git clone --depth 1 -b dev https://github.com/complextissue/pytximport.git
cd pytximport
pyenv local 3.12
make create-venv
source .venv/source/activate
make install-dev
```

:::

::::
