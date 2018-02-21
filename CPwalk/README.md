Function that allows interpolation of Range BLRMS filter outputs in 2D as a function of static misalignments of CPs in PIT and YAW.
 - Requirements are in `reqs.txt`
 - Uses a `yaml` parameter file to specify all preferences.
 - Supports plotting from data from an ASCII file, else the data can be downloaded from frames using `cdsutils` (not tested).
 - For the interpolation, `scipy.interpolate.griddata` or `scipy.interpolate.Rbf` is used.
 - The former interpolates only inside the convex hull, while the latter extrapolates for points outside the convex hull.
 - In general, I find that the approach requires considerable **INTERPRETATION** of the results.
 - Linear interpolation is required to guarantee non-negative interpolated values (which is the only physical possibility, since we are interpolating a BLRMS value).
