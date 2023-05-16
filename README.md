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

For a better explanation of each algorithm see [Polyline Simplification](https://www.codeproject.com/Articles/114797/Polyline-Simplification) which also includes graphical representations to aid in understanding.


### `nth_point(curve, n)`
Implementation of the nth point algorithm. Takes two `intent(in)` parameters:  
- `curve` - a `real(dp),dimension(:)` or `real(dp),dimension(:,:)` array containing the points of the curve to be simplified. If `real(dp), dimension(:,:)`, then `curve(i,:)` denotes each point of the curve (column major).
- `n` - a default `integer` that is interpreted as every nth element to be kept.  

This function is of the same type as `curve`.

### `radial_distance(curve, tolerance)`
Implementation of the radial distance algorithm. Takes two `intent(in)` parameters:
- `curve` - a `real(dp),dimension(:)` or `real(dp),dimension(:,:)` array containing the points of the curve to be simplified. If `real(dp), dimension(:,:)`, then `curve(i,:)` denotes each point of the curve (column major).
- `tolerance` - a `real(dp)` value interpreted as the minimum distance required between any two consecutive points in the resulting curve.

This function is of the same type as `curve`.

### `perpendicular_distance(curve, tolerance, repeat)`
Implementation of the perpendicular distance algorithm. Takes up to three `intent(in)` parameters:
- `curve` - a `real(dp),dimension(:,:)` array containing the points of the curve to be simplified. `curve(i,:)` denotes each point of the curve (column major).
- `tolerance` - a `real(dp)` value interpreted as the minimum distance a point can be from the line defined by the most recent confirmed key and the next point in the original curve.
- `repeat` - an optional, default `integer` that describes how many times to run the algorithm. This parameter is useful as the algorithm will at most remove 50% of the points from a given curve. 

This function is of the same type as `curve`.

### `reumann_witkam(curve, tolerance)`
Implementation of the Reumann-Witkam algorithm. Takes two parameters:
- `curve` - a `real(dp),dimension(:,:)` array containing the points of the curve to be simplified. `curve(i,:)` denotes each point of the curve (column major).
- `tolerance` - a `real(dp)` value interpreted as the minimum distance points must be from any line segment.

This function is of the same type as `curve`.