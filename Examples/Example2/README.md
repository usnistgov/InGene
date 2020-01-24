This example takes in a family of conditional output distributions (expressions*.csv). Each expressions<n>_<k>.csv file refers to the k-th replicate of the output distribution obtained using 1/n subsample (with replacement) of the total output data. The sweep_landscapes.py computes the mutual information landscapes for each replicate of the each subsampled conditional output distributions.

The LimitingLandscape.py file in postprocessors uses this family of mutual information landscapes to obtain the correct one without the effect of finite-sampling.
