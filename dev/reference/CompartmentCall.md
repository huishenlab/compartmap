# Compartment Call objects for analysis (experimental)

These S7 classes hold compartment call values, and their genomic region
and resolution metadata as well as methods to aid analysis and
visualization. `CompartmentCall()` is primarily a parent class from
which `CompartmapCall()` and `MultiCompartmapCall()` are derived.

## Usage

``` r
CompartmentCall(pc, res, gr, name = NULL, unitarized = FALSE)

CompartmapCall(gr, res, name = NULL, unitarized = FALSE)

MultiCompartmapCall(
  ccalls,
  name,
  unitarize = FALSE,
  filter = FALSE,
  filter_threshold = 0.02
)

scCompartmapCall(re, res, name, unitarized = FALSE, unitarize = FALSE)
```

## Arguments

- pc:

  The singular values from a compartment call

- res:

  The binning resolution used

- gr:

  The GRanges of the bins or the output of `scCompartments` or
  `getArrayCompartments` containing the 'pc' column

- name:

  An identifier for this set of compartment calls

- unitarized:

  Whether the singular values have been unitarized

- ccalls:

  A list of `CompartmapCalls` to combine

- unitarize:

  Whether to unitarize the singular values for each of the inputs calls

- filter:

  Whether to filter the singular values

- filter_threshold:

  The threshold at which to filter the singular values

- re:

  A `RaggedExperiment` of single-cell compartment calls

## Details

### Objects

`CompartmentCall` is constucted from a vector of the call values, a
[`GenomicRanges::GRanges()`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object of the bins, and the resolution of the bins. Given a `GRanges` of
bin coordinates, it can be used to store singular values/eigenvectors
from Hi-C compartment analysis. Along with `CompartmapCall` objects,
they can be combined in a `MultiCompartmapCall` to compare against each
other. All `CompartmentCall` methods work on every other *compartmap* S7
object.

`CompartmapCall()` is constructed from the `GRanges` output of
`scCompartments(group = TRUE)`. `MultiCompartmapCall()` holds multiple
`CompartmapCall()` objects at the same resolution, and is constructed
with a list of `CompartmentCall` or `CompartmapCall` objects.

For single-cell level compartment inferences use `scCompartmapCall()`,
which inherits from `MultiCompartmapCall`. This class holds multiple
single-cell level compartment calls, taking the
[`RaggedExperiment::RaggedExperiment()`](https://bioconductor.github.io/RaggedExperiment/reference/RaggedExperiment-class.html)
from
[`scCompartments()`](https://huishenlab.github.io/compartmap/dev/reference/scCompartments.md)
as its constructor. Since it inherits from `MultiCompartmapCall`, all
`MultiCompartmapCall` methods work on `scCompartmapCall` objects.

All objects may be subset to required indices with
`object[row indices, column indices]`. Their dimensions can be accessed
with [`dim()`](https://rdrr.io/r/base/dim.html),
[`nrow()`](https://rdrr.io/r/base/nrow.html) and
[`ncol()`](https://rdrr.io/r/base/nrow.html).

### Properties and accessors

These properties may be accessed with `@` like S4 slots, or using their
accessor functions. Some functions share methods with other Bioconductor
classes like `GRanges` and `SummarizedExperiment`.

- `name`,
  [`get_name()`](https://huishenlab.github.io/compartmap/dev/reference/get_name.md):
  The name of the object, used to identify it, useful when making
  `MultiCompartmapCall` objects

- `gr`, `granges()`: A `GRanges` object of the compartment bins

- `df`,
  [`DF()`](https://huishenlab.github.io/compartmap/dev/reference/DF.md):
  a `data.table` of bin indices in column `n` and compartment call
  singular values in column `pc`. For `MultiCompartmapCall` and
  `scCompartmapCall` objects, this is in a tidy format, with an
  additional `name` column.

- `res`,
  [`resolution()`](https://huishenlab.github.io/compartmap/dev/reference/resolution.md):
  The genomic resolution at which the compartments were called

- `unitarized`,
  [`is_unitarized()`](https://huishenlab.github.io/compartmap/dev/reference/is_unitarized.md):
  Whether the singular values have been unitarized

- `seqinfo`, `seqinfo()`: Seqinfo for the object's `GRanges`

- `mat`,
  [`mat()`](https://huishenlab.github.io/compartmap/dev/reference/mat.md):
  A matrix of genomic bins by singular values - a wide matrix format of
  the `df` slot. This is used to calculate correlation and agreement
  between cells in a `scCompartmapCall` and groups in a
  `MultiCompartmapCall` object.
