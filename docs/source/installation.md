# Installation

:::{card}
:class-card: sd-bg-warning
:class-body: sd-bg-text-warning
**pytximport** only supports Python versions greater than or equal to **3.9**.
:::

:::{card} Recommendation
Installation via `Bioconda` is the recommended way to include `pytximport` in your projects. A `pip`-installabe package
is also provided.
:::

## Installation Options

Choose an option to install this package.

::::{tab-set}

:::{tab-item} Bioconda
Install `pytximport` from `Bioconda` using `mamba` or `conda`:

```bash
mamba install -c bioconda pytximport
```

:::

:::{tab-item} PyPi
Install `pytximport` package using `pip`:

```bash
python3 -m pip install pytximport
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
git clone https://github.com/complextissue/pytximport.git
cd pytximport
pyenv local 3.9
make create-venv
source .venv/source/activate
make install-dev
```

:::

::::
