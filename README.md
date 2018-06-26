# UniDip Python Port

See reference paper:  http://www.kdd.org/kdd2016/subtopic/view/skinny-dip-clustering-in-a-sea-of-noise

UniDip is a noise robust clustering algorithm for 1 dimensional numeric data. It recursively extracts peaks of density in the data utilizing the Hartigan Dip-test of Unimodality.

## Install

coming soon...
```
pip3.6 install unidip
```

## Examples

### Basic Usage

```python
from unidip import UniDip

# create bi-modal distribution
dat = np.concatenate([np.random.randn(200)-3, np.random.randn(200)+3])

# sort data so returned indices are meaningful
dat = np.msort(dat)

# get start and stop indices of peaks 
intervals = UniDip(dat).run()
```

### Advanced Options

* **alpha**: control sensitivity as p-value. Default is 0.05. increase to isolate more peaks with less confidence. Or, decrease to isolate only peaks that are least likely to be noise.
* **mrg_dst**: Defines how close intervals must be before they are merged.
* **ntrials**: how many trials are run in Hartigan Dip Test more trials adds confidance but takes longer.

```python
intervals = UniDip(dat, alpha=0.001, ntrials=1000, mrg_dst=5).run()
```