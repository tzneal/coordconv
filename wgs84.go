package coordconv

import "fmt"

// DefaultMGRSConverter is a WGS84 ellipsoid based MGRS converter.
var DefaultMGRSConverter *MGRS

// DefaultUTMConverter is a WGS84 ellipsoid based UTM converter.
var DefaultUTMConverter *UTM

// DefaultUPSConverter is a WGS84 ellipsoid based UPS converter.
var DefaultUPSConverter *UPS

func init() {
	const semiMajorAxis = 6378137
	const flattening = 1 / 298.257223563
	var err error
	DefaultMGRSConverter, err = NewMGRS(semiMajorAxis, flattening, "WE")
	if err != nil {
		panic(fmt.Sprintf("error constructing WGS84 MGRS converter: %s", err))
	}
	DefaultUTMConverter, err = NewUTM()
	if err != nil {
		panic(fmt.Sprintf("error constructing WGS84 UTM converter: %s", err))
	}
	DefaultUPSConverter, err = NewUPS(semiMajorAxis, flattening)
	if err != nil {
		panic(fmt.Sprintf("error constructing WGS84 UPS converter: %s", err))
	}
}
