
coordconv
=========

coordconv is a manual port of the Geotrans 3.7 MGRS conversion code from C++ to
Go. It also makes public the UTM and UPS converters as well as the underlying
transverse mercator and polar stereographic projections.

It has been tested with 4million+ points uniformly spread across the globe and
returns identical results to the standard Geotrans library.  The unit tests
included contain an abbreviated test table at testdata/geotrans.log with 40,000
points which was generated with the included geotrans_table_cpp.  It can be
modified and used to regenerate the included geotrans.log, or a more
comprehensive one if desired.

Usage
=====

Constructors are provided so you can use custom ellipsoid parameters, but defaults are
provided for WGS84:

```go
  mgrs, _ := coordconv.DefaultMGRSConverter.ConvertFromGeodetic(s2.LatLngFromDegrees(0, 0), 5)
  fmt.Println(mgrs) // 31NAA6602100000
  geo, _ := coordconv.DefaultMGRSConverter.ConvertToGeodetic("16SGC3855124838")
  fmt.Println(geo) // [33.6366624, -84.4280571]
```

License
=======

The [geotrans code](http://earth-info.nga.mil/GandG/update/index.php?action=home#tab_wgs84-data) is public domain, so in spirit this is public domain as well.  I only made minor changes along the way to make it behave a bit more like like standard Go code (e.g. returning an error when the C++ code threw an exception)
