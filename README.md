Boba-seq is a barcoded approach for large-scale overexpression
screens. See

https://doi.org/10.1101/2022.10.10.511384

bobaseq.R in this code base defines the R functions to turn counts per
barcode into fitness values for each fragment and to identify
biologically consistent hits. There is also a function to plot the
per-fragment fitness values.

To use long reads to associate barcodes with insert locations, see

https://github.com/OGalOz/Boba-seq

To count barcodes from BarSeq data, see bin/MultiCodes.pl in

https://bitbucket.org/berkeleylab/feba/src
