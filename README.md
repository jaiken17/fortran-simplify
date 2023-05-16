# fortran-simplify
Fortran module that reduces the over-sampled resolution of a polyline. This process can be useful when working with data that is noisy but usable with a significantly reduced resolution.


## How To Use
The `simplify` module is all that is necessary to use the polyline simplification algorithms. The `simplify` module does depend on a `precision` module that at least needs to contain a real kind named `dp`. This can easily be modified to accept a real kind with another name in a different module by changing the line in `simplify.f90` from
```Fortran
use precision
```
to
```Fortran
use real_kinds_module, only : dp => real_kind
```
All that is left is to add a `use` statement for the module in whichever program units you need the algorithms and to compile the module as normal.

## `polytest` Test Program



## Function Signatures And Descriptions