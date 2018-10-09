package coordconv

import (
	"errors"
	"math"

	"github.com/golang/geo/s2"
)

// Hemisphere represents the hemisphere, north or south
type Hemisphere byte

// Hemisphere constants
const (
	HemisphereInvalid Hemisphere = iota
	HemisphereNorth
	HemisphereSouth
)

// UPSCoord is a UPS coordinate with a specified easting/northing in meters and
// hemisphere.
type UPSCoord struct {
	Hemisphere Hemisphere
	Easting    float64
	Northing   float64
}

// UPS is a UPS coordinate converter
type UPS struct {
	semiMajorAxis          float64
	flattening             float64
	UPSOriginLatitude      float64
	polarStereographicMapN *PolarStereographic
	polarStereographicMapS *PolarStereographic
}

const epsilonRadians = 1.75e-7 // approx 1.0e-5 degrees (~1 meter) in radians

const upsFalseEasting = 2000000
const upsFalseNorthing = 2000000
const upsOriginLatitude = 0.0

const upsMaxLat = 90.0 * (math.Pi / 180.0) // 90 degrees in radians
const upsMaxOriginLat = 81.114528 * (math.Pi / 180.0)
const upsMinNorthLat = 83.5 * (math.Pi / 180.0)
const upsMaxSouthLat = -79.5 * (math.Pi / 180.0)
const upsMinEastNorth = 0.0
const upsMaxEastNorth = 4000000.0

// NewUPS construct a new UPS converter with the specified ellipsoid parameters.
func NewUPS(ellipsoidSemiMajorAxis, ellipsoidFlattening float64) (*UPS, error) {
	invF := 1 / ellipsoidFlattening
	if ellipsoidSemiMajorAxis <= 0.0 {
		return nil, errors.New("Semi-major axis must be greater than zero")
	}
	if (invF < 250) || (invF > 350) {
		return nil, errors.New("Inverse flattening must be between 250 and 350")
	}

	u := &UPS{
		UPSOriginLatitude: upsMaxOriginLat,
	}

	u.semiMajorAxis = ellipsoidSemiMajorAxis
	u.flattening = ellipsoidFlattening
	u.polarStereographicMapN, _ = NewPolarStereographicScaleFactor(u.semiMajorAxis, u.flattening, upsOriginLatitude,
		.994, HemisphereNorth, upsFalseEasting, upsFalseNorthing)

	u.polarStereographicMapS, _ = NewPolarStereographicScaleFactor(u.semiMajorAxis, u.flattening, upsOriginLatitude,
		.994, HemisphereSouth, upsFalseEasting, upsFalseNorthing)
	return u, nil
}

// ConvertFromGeodetic converts a geodetic coordinate to a UPS coordinate.
func (u *UPS) ConvertFromGeodetic(geodeticCoordinates s2.LatLng) (UPSCoord, error) {
	longitude := geodeticCoordinates.Lng.Radians()
	latitude := geodeticCoordinates.Lat.Radians()

	if (latitude < -upsMaxLat) ||
		(latitude > upsMaxLat) {
		return UPSCoord{}, errors.New("latitude out of range")
	} else if (latitude < 0) && (latitude >= (upsMaxSouthLat + epsilonRadians)) {
		return UPSCoord{}, errors.New("latitude out of range")
	} else if (latitude >= 0) && (latitude < (upsMinNorthLat - epsilonRadians)) {
		return UPSCoord{}, errors.New("latitude out of range")
	}
	if (longitude < -math.Pi) ||
		(longitude > (2 * math.Pi)) {
		return UPSCoord{}, errors.New("longitude out of range")
	}

	var polarStereographic *PolarStereographic
	var hemisphere Hemisphere
	if latitude < 0 {
		u.UPSOriginLatitude = -upsMaxOriginLat
		hemisphere = HemisphereSouth
		polarStereographic = u.polarStereographicMapS
	} else {
		u.UPSOriginLatitude = upsMaxOriginLat
		hemisphere = HemisphereNorth
		polarStereographic = u.polarStereographicMapN
	}

	polarStereographicCoordinates, _ := polarStereographic.ConvertFromGeodetic(geodeticCoordinates)
	easting := polarStereographicCoordinates.Easting
	northing := polarStereographicCoordinates.Northing

	return UPSCoord{
		Hemisphere: hemisphere,
		Easting:    easting,
		Northing:   northing,
	}, nil

}

// ConvertToGeodetic converts UPS (hemisphere, easting, and northing)
// coordinates to geodetic (latitude and longitude) coordinates according to the
// current ellipsoid parameters.
func (u *UPS) ConvertToGeodetic(upsCoordinates UPSCoord) (s2.LatLng, error) {
	hemisphere := upsCoordinates.Hemisphere
	easting := upsCoordinates.Easting
	northing := upsCoordinates.Northing

	if (hemisphere != HemisphereNorth) && (hemisphere != HemisphereSouth) {
		return s2.LatLng{}, errors.New("hemisphere invalid")
	}

	if (easting < upsMinEastNorth) || (easting > upsMaxEastNorth) {
		return s2.LatLng{}, errors.New("easting out of range")
	}
	if (northing < upsMinEastNorth) || (northing > upsMaxEastNorth) {
		return s2.LatLng{}, errors.New("northing out of range")
	}

	if hemisphere == HemisphereNorth {
		u.UPSOriginLatitude = upsMaxOriginLat
	} else {
		u.UPSOriginLatitude = -upsMaxOriginLat
	}

	polarStereographicCoordinates := MapCoords{
		Easting:  easting,
		Northing: northing,
	}
	var polarStereographic *PolarStereographic
	if hemisphere == HemisphereNorth {
		polarStereographic = u.polarStereographicMapN
	} else {
		polarStereographic = u.polarStereographicMapS
	}
	geodeticCoordinates, err := polarStereographic.ConvertToGeodetic(polarStereographicCoordinates)
	if err != nil {
		return s2.LatLng{}, err
	}

	latitude := geodeticCoordinates.Lat.Radians()

	if (latitude < 0) && (latitude >= (upsMaxSouthLat + epsilonRadians)) {
		return s2.LatLng{}, errors.New("resulting latitude out of range")
	}
	if (latitude >= 0) && (latitude < (upsMinNorthLat - epsilonRadians)) {
		return s2.LatLng{}, errors.New("resulting latitude out of range")
	}

	return geodeticCoordinates, nil
}
