# snATAC Utilities
A collection of functions for snATAC-seq analysis built around snapATAC2.
Includes 3 modules:

1. Preprocessing
   - Identify genomic bins overlapping with features from provided bed files
   - Count features overlapping genomic bins
   - Count & remove doublets

2. Tools
   - Downsample counts in cells
   - Specteral embedding to umap in 1 function

3. Plotting
   - Plot TSSe count plots for TSSe thresholds
   - Plot cell counts for observations
   - Plot scatter plots of observations

## Installation

Clone this repository, and follow the commands below.

```bash
git clone https://github.com/Ditchburn27/snATAC_utilities.git

cd snATAC_utilities

pip install .
```