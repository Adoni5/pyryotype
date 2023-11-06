# pyryotype

[![PyPI - Version](https://img.shields.io/pypi/v/pyryotype.svg)](https://pypi.org/project/pyryotype)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pyryotype.svg)](https://pypi.org/project/pyryotype)

-----

**Table of Contents**

- [pyryotype](#pyryotype)
  - [Installation](#installation)
  - [Example usage](#example-usage)
    - [Coverage plotting](#coverage-plotting)
    - [PAF plotting](#paf-plotting)
  - [License](#license)
  - [Cytoband data](#cytoband-data)

## Installation

```console
    pip install ideogram
```

## Example usage

```python
from pyryotype import GENOME, plot_ideogram
from matplotlib import pyplot as plt
fig, axes = plt.subplots(
    ncols=1,
    nrows=22,
    figsize=(11, 11),
    facecolor="white",
)
genome = GENOME.CHM13
for ax, contig_name in zip(axes, range(1, 23)):
    chromosome = f"chr{contig_name}"
    plot_ideogram(ax, target=chromosome, genome=genome)
fig.savefig("ideogram.png", dpi=300)
```

Will output:
![Example ideogram](https://raw.githubusercontent.com/Adoni5/pyryotype/d724012befec0b56351d0db5125f8d9cf4df1816/example_outputs/ideogram.png?raw=true)

### Coverage plotting
Coverage plotting is designed to be used with the output of [Mosdepth](https://github.com/brentp/mosdepth). The following example uses the output of `mosdepth` to plot the coverage of chromosome 1. The region representing the first 100Mb of the chromosome is highlighted in black.

```python
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.ticker import EngFormatter
import pandas as pd

from pyryotype.coverage import plot_coverage


test_cov_file = Path(__file__).parent / "static" / "chr1_cov.regions.bed.gz"

fig, ax = plt.subplots(
    ncols=1,
    nrows=1,
    figsize=(11, 2),
)
df = pd.read_csv(test_cov_file, sep="\t", names=["chromosome", "start", "end", "value"])
ax = plot_coverage(ax, df, "chr1", regions=[(0, 100000000)], ylabel="Coverage", color="black")
ax.set_xlabel("Genomic Position (bp)")

# ax.yaxis.set_visible(False)
ax.set_yticks([])
ax.set_yticklabels([])
for side in ("right", "top", "left"):
    ax.spines[side].set_visible(False)
formatter = EngFormatter(unit='b', places=1)
# Set formatter for the x-axis
ax.xaxis.set_major_formatter(formatter)

fig.savefig("example_outputs/test_coverage.png", dpi=300, bbox_inches="tight")

```
Will output something like:
![Example Coverage](https://raw.githubusercontent.com/Adoni5/pyryotype/main/example_outputs/test_coverage.png)

### PAF plotting

Designed primarily to plot the alignment of assemblies to a reference genome. Must provide a valid PAF file.

```python
from matplotlib import pyplot as plt
from readpaf import parse_paf
from pathlib import Path
from pyryotype.paf_plotting import PAFProtocol, PlotMode, plot_paf_alignments
test_paf = Path(__file__).parent / "static" / "test.paf"
fig, ax = plt.subplots(
    ncols=1,
    nrows=1,
    figsize=(11, 1),
)

ax = plot_paf_alignments(
    ax,
    parse_paf(test_paf.open()),
    target="chr1",
    mapq_filter=0,
    strict=PlotMode.STRICT,
    contig_colours=PlotMode.UNIQUE_COLOURS,
)
ax.set_xlabel("")

# ax.yaxis.set_visible(False)
ax.set_yticks([])
ax.set_yticklabels([])
for side in ("right", "top", "left"):
    ax.spines[side].set_visible(False)
fig.savefig("tests/test_paf_plotting.png", dpi=300, bbox_inches="tight")
```

Will output the following image:
![Example PAF plotting](https://raw.githubusercontent.com/Adoni5/pyryotype/d724012befec0b56351d0db5125f8d9cf4df1816/example_outputs/test_paf_plotting.png?raw=true)
The colours assigned to each alignment can be changed to either be based on the Alignment Strand, or unique for each record. See the PlotMode [docstring](https://github.com/Adoni5/pyryotype/blob/0517a8805aac7b00bdddc7d612c2c80c56b6891c/src/pyryotype/paf_plotting.py#L366) for the options that can be applied.

Multiple mappings on the same chromosome can be collapsed into a single line by setting `strict=PlotMode.STRICT`.
This can be seen in the above mappings, where the large block is comprised of 3 separate alignments, from the same read. IF `strict=PlotMode.CHILL` then each alignment will be plotted separately, even if these alignments are from the same read. This looks like:
![Example PAF plotting chill](https://raw.githubusercontent.com/Adoni5/pyryotype/d724012befec0b56351d0db5125f8d9cf4df1816/example_outputs/test_paf_plotting_chill.png?raw=true)


## License

`pyryotype` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.

## Cytoband data
* HG38 - Nushell, will have to be adapted for bash `curl -L "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz" | gzip -d - | rg -Ne "^chr\\d+\t" | save cytoBand_HG38.tsv`
* CHM13 - bash yay `curl -L http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v2.0/download/chm13v2.0_cytobands_allchrs.bed.gz | gzip -d - > cytobands_chm13.bed`
