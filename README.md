# fortran-simplify
Fortran module that reduces the over-sampled resolution of a polyline. This process can be useful when working with data that is noisy but usable with a significantly reduced resolution.

</br>

## How To Use
The `simplify` module is all that is necessary to use the polyline simplification algorithms. The `simplify` module does depend on a `smpl_precision` module that at least needs to contain a real kind named `dp`. In this case, it uses the double precision real kind from the `stdlib_kinds` module found in https://github.com/fortran-lang/stdlib. This can easily be modified to accept a real kind with another name in a different module by changing the line in `simplify.f90` from
```Fortran
use smpl_precision
```
to
```Fortran
use real_kinds_module, only : dp => real_kind
```
All that is left is to add a `use` statement for the module in whichever program units you need the algorithms and to compile the module as normal.

To add fortran-simplify to your fpm project, simply add the following line to the `fpm.toml` file under the `[dependencies]` tag:
```toml
fortran-simplify.git = "https://github.com/jaiken17/fortran-simplify"
```


</br>

## `simplify_example` Example Program
Included is a simple test program that runs the algorithms on a couple of curves. The two dimensional curve is also written to a file `curve.data` and the output from the perpendicular distance and Reumann-Witkam algorithms are written to `perp_simple_curve.data` and `ruemann_witkam_simple_curve.data`, respectively. 

With the Fortran Package Manager ([fpm](https://github.com/fortran-lang/fpm)) installed, compile and run the example program in a shell by excuting the following command in the top-level project directory:
```bash
fpm run --example
```

The data files are then created in the example directory and can then be plotted with gnuplot by running the plot bash script (within the example directory):
```bash
./plot
```


</br>

## Function Signatures And Descriptions

For a better explanation of each algorithm see [Polyline Simplification](https://www.codeproject.com/Articles/114797/Polyline-Simplification) which also includes graphical representations to aid in understanding.

</br>

### `nth_point(curve, n)`
Implementation of the nth point algorithm. Takes two `intent(in)` parameters:  
- `curve` - a `real(dp),dimension(:)` or `real(dp),dimension(:,:)` array containing the points of the curve to be simplified. If `real(dp), dimension(:,:)`, then `curve(i,:)` denotes each point of the curve (column major).
- `n` - a default `integer` that is interpreted as every nth element to be kept.  

This function is of the same type as `curve`.

</br>

### `radial_distance(curve, tolerance)`
Implementation of the radial distance algorithm. Takes two `intent(in)` parameters:
- `curve` - a `real(dp),dimension(:)` or `real(dp),dimension(:,:)` array containing the points of the curve to be simplified. If `real(dp), dimension(:,:)`, then `curve(i,:)` denotes each point of the curve (column major).
- `tolerance` - a `real(dp)` value interpreted as the minimum distance required between any two consecutive points in the resulting curve.

This function is of the same type as `curve`.

</br>

### `perpendicular_distance(curve, tolerance, repeat)`
Implementation of the perpendicular distance algorithm. Takes up to three `intent(in)` parameters:
- `curve` - a `real(dp),dimension(:,:)` array containing the points of the curve to be simplified. `curve(i,:)` denotes each point of the curve (column major).
- `tolerance` - a `real(dp)` value interpreted as the minimum distance a point can be from the line defined by the most recent confirmed key and the next point in the original curve.
- `repeat` - an optional, default `integer` that describes how many times to run the algorithm. This parameter is useful as the algorithm will at most remove 50% of the points from a given curve. 

This function is of the same type as `curve`.

</br>

### `reumann_witkam(curve, tolerance)`
Implementation of the Reumann-Witkam algorithm. Takes two parameters:
- `curve` - a `real(dp),dimension(:,:)` array containing the points of the curve to be simplified. `curve(i,:)` denotes each point of the curve (column major).
- `tolerance` - a `real(dp)` value interpreted as the minimum distance points must be from any line segment.

This function is of the same type as `curve`.
