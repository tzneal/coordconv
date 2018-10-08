package coordconv

import (
	"errors"
	"math"

	"github.com/golang/geo/s1"
	"github.com/golang/geo/s2"
)

// UTMCoord is a UTM coordinate
type UTMCoord struct {
	Zone       int
	Hemisphere Hemisphere
	Easting    float64
	Northing   float64
}

// UTM is a UTM coordinate converter
type UTM struct {
	CoordinateSystem
	ellipsCode            string
	utmOverride           int
	transverseMercatorMap [61]*TransverseMercator
}

const utmMinLat = ((-80.5 * math.Pi) / 180.0) // -80.5 degrees in radians
const utmMaxLat = ((84.5 * math.Pi) / 180.0)  //  84.5 degrees in radians
const utmMinEasting = 100000.0
const utmMaxEasting = 900000.0
const utmMinNorthing = 0.0
const utmMaxNorthing = 10000000.0

// NewUTM constructs a new UTM converter for the WGS84 ellipsoid
func NewUTM() (*UTM, error) {
	u := &UTM{
		ellipsCode: "WE",
	}

	u.semiMajorAxis = 6378137.0
	u.flattening = 1 / 298.257223563
	u.utmOverride = 0

	var centralMeridian float64
	originLatitude := 0.0
	falseEasting := 500000.0
	falseNorthing := 0.0
	scale := 0.9996

	for zone := 1; zone <= 60; zone++ {
		if zone >= 31 {
			centralMeridian = (float64(6*zone-183) * math.Pi / 180)
		} else {
			centralMeridian = (float64(6*zone+177) * math.Pi / 180)
		}

		var err error
		u.transverseMercatorMap[zone], err = NewTransverseMercator(u.semiMajorAxis, u.flattening, centralMeridian, originLatitude,
			falseEasting, falseNorthing, scale, u.ellipsCode)
		if err != nil {
			return nil, err
		}
	}
	return u, nil
}

// NewUTM2 receives the ellipsoid parameters and UTM zone override parameter as
// inputs, and sets the corresponding state variables.  override is the UTM
// override zone, 0 indicates no override.
func NewUTM2(ellipsoidSemiMajorAxis, ellipsoidFlattening float64,
	ellipsoidCode string, override int) (*UTM, error) {
	u := &UTM{
		ellipsCode: ellipsoidCode,
	}
	invF := 1 / ellipsoidFlattening

	if ellipsoidSemiMajorAxis <=
		0.0 {
		return nil, errors.New("Semi-major axis must be greater than zero")
	}
	if (invF < 250) ||
		(invF > 350) {
		return nil, errors.New("Inverse flattening must be between 250 and 350")
	}
	if (override < 0) || (override > 60) {
		return nil, errors.New("zone override out of range")
	}

	u.semiMajorAxis = ellipsoidSemiMajorAxis
	u.flattening = ellipsoidFlattening

	u.utmOverride = override

	var centralMeridian float64
	originLatitude := 0.0
	falseEasting := 500000.0
	falseNorthing := 0.0
	scale := 0.9996

	for zone := 1; zone <= 60; zone++ {
		if zone >= 31 {
			centralMeridian = (float64(6*zone-183) * math.Pi / 180)
		} else {
			centralMeridian = (float64(6*zone+177) * math.Pi / 180)
		}

		var err error
		u.transverseMercatorMap[zone], err = NewTransverseMercator(
			u.semiMajorAxis, u.flattening, centralMeridian, originLatitude,
			falseEasting, falseNorthing, scale, u.ellipsCode)
		if err != nil {
			return nil, err
		}
	}
	return u, nil
}

// ConvertFromGeodetic converts geodetic (latitude and longitude) coordinates
// to UTM projection (zone, hemisphere, easting and northing) coordinates
// according to the current ellipsoid and UTM zone override parameters.
func (u *UTM) ConvertFromGeodetic(geodeticCoordinates s2.LatLng, utmZoneOverride int) (UTMCoord, error) {
	FalseNorthing := 0.0
	longitude := geodeticCoordinates.Lng.Radians()
	latitude := geodeticCoordinates.Lat.Radians()
	if (latitude < (utmMinLat - epsilonRadians)) ||
		(latitude >= (utmMaxLat + epsilonRadians)) {
		return UTMCoord{}, errors.New("latitude out of range ")
	}
	if (longitude < (-math.Pi - epsilonRadians)) ||
		(longitude > (2*math.Pi + epsilonRadians)) {
		return UTMCoord{}, errors.New("longitude out of range")
	}

	if (latitude > -1.0e-9) && (latitude < 0) {
		latitude = 0.0
	}

	if longitude < 0 {
		longitude += (2 * math.Pi)
	}

	LatDegrees := int(latitude * 180.0 / math.Pi)
	LongDegrees := int(longitude * 180.0 / math.Pi)

	var tempZone int
	if longitude < math.Pi {
		tempZone = int(31 + (((longitude + 1.0e-10) * 180.0 / math.Pi) / 6.0))
	} else {
		tempZone = int((((longitude + 1.0e-10) * 180.0 / math.Pi) / 6.0) - 29)
	}

	if tempZone > 60 {
		tempZone = 1
	} else if tempZone < 0 {
		return UTMCoord{}, errors.New("longitude out of range")
	}

	// allow UTM zone override up to +/- one zone of the calculated zone
	if utmZoneOverride != 0 {
		if (tempZone == 1) && (utmZoneOverride == 60) {
			tempZone = utmZoneOverride
		} else if (tempZone == 60) && (utmZoneOverride == 1) {
			tempZone = utmZoneOverride
		} else if ((tempZone - 1) <= utmZoneOverride) &&
			(utmZoneOverride <= (tempZone + 1)) {
			tempZone = utmZoneOverride
		} else {
			return UTMCoord{}, errors.New("zone out of range")
		}
	} else if u.utmOverride != 0 {
		if (tempZone == 1) && (u.utmOverride == 60) {
			tempZone = u.utmOverride
		} else if (tempZone == 60) && (u.utmOverride == 1) {
			tempZone = u.utmOverride
		} else if ((tempZone - 1) <= u.utmOverride) &&
			(u.utmOverride <= (tempZone + 1)) {
			tempZone = u.utmOverride
		} else {
			return UTMCoord{}, errors.New("zone out of range")
		}
	} else { // not UTM zone override
		// check for special zone cases over southern Norway and Svalbard
		if (LatDegrees > 55) && (LatDegrees < 64) && (LongDegrees > -1) &&
			(LongDegrees < 3) {
			tempZone = 31
		}
		if (LatDegrees > 55) && (LatDegrees < 64) && (LongDegrees > 2) &&
			(LongDegrees < 12) {
			tempZone = 32
		}
		if (LatDegrees > 71) && (LongDegrees > -1) && (LongDegrees < 9) {
			tempZone = 31
		}
		if (LatDegrees > 71) && (LongDegrees > 8) && (LongDegrees < 21) {
			tempZone = 33
		}
		if (LatDegrees > 71) && (LongDegrees > 20) && (LongDegrees < 33) {
			tempZone = 35
		}
		if (LatDegrees > 71) && (LongDegrees > 32) && (LongDegrees < 42) {
			tempZone = 37
		}
	}

	transverseMercator := u.transverseMercatorMap[tempZone]
	var hemisphere Hemisphere
	if latitude < 0 {
		FalseNorthing = 10000000
		hemisphere = HemisphereSouth
	} else {
		hemisphere = HemisphereNorth
	}
	tempGeodeticCoordinates := s2.LatLng{Lng: s1.Angle(longitude), Lat: s1.Angle(latitude)}
	transverseMercatorCoordinates, err := transverseMercator.convertFromGeodetic(tempGeodeticCoordinates)
	if err != nil {
		return UTMCoord{}, err
	}
	easting := transverseMercatorCoordinates.Easting
	northing := transverseMercatorCoordinates.Northing + FalseNorthing
	if (easting < utmMinEasting) || (easting > utmMaxEasting) {
		return UTMCoord{}, errors.New("easting out of range")
	}

	if (northing < utmMinNorthing) || (northing > utmMaxNorthing) {
		return UTMCoord{}, errors.New("northing out of range")
	}

	return UTMCoord{
		Zone:       tempZone,
		Hemisphere: hemisphere,
		Easting:    easting,
		Northing:   northing,
	}, nil
}

// ConvertToGeodetic converts UTM projection (zone, hemisphere, easting and
// northing) coordinates to geodetic(latitude and  longitude) coordinates,
// according to the current ellipsoid parameters.
func (u *UTM) ConvertToGeodetic(utmCoordinates UTMCoord) (s2.LatLng, error) {
	FalseNorthing := 0.0
	zone := utmCoordinates.Zone
	hemisphere := utmCoordinates.Hemisphere
	easting := utmCoordinates.Easting
	northing := utmCoordinates.Northing

	if (zone < 1) || (zone > 60) {
		return s2.LatLng{}, errors.New("zone out of range")
	}
	if (hemisphere != HemisphereSouth) && (hemisphere != HemisphereNorth) {
		return s2.LatLng{}, errors.New("hemisphere out of range")
	}
	if (easting < utmMinEasting) || (easting > utmMaxEasting) {
		return s2.LatLng{}, errors.New("easting out of range")
	}
	if (northing < utmMinNorthing) || (northing > utmMaxNorthing) {
		return s2.LatLng{}, errors.New("northing out of range")
	}

	transverseMercator := u.transverseMercatorMap[zone]

	if hemisphere == HemisphereSouth {
		FalseNorthing = 10000000
	}

	transverseMercatorCoordinates := MapCoords{Easting: easting, Northing: northing - FalseNorthing}
	geodeticCoordinates, err := transverseMercator.convertToGeodetic(transverseMercatorCoordinates)
	if err != nil {
		return s2.LatLng{}, err
	}

	latitude := geodeticCoordinates.Lat.Radians()
	if (latitude < (utmMinLat - epsilonRadians)) ||
		(latitude >= (utmMaxLat + epsilonRadians)) {
		return s2.LatLng{}, errors.New("latitude out of range ")
	}

	return geodeticCoordinates, nil
}
