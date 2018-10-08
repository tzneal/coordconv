package coordconv

import (
	"errors"
	"math"

	"github.com/golang/geo/s1"
	"github.com/golang/geo/s2"
)

type PolarStereographic struct {
	semiMajorAxis        float64
	flattening           float64
	es                   float64 // Eccentricity of ellipsoid
	esOverTwo            float64 // es / 2.0
	isSouthernHemisphere bool    // Flag variable
	polarTC              float64
	polarK90             float64
	polaraMc             float64 // Polar_a * mc
	twoPolarA            float64 // 2.0 * Polar_a

	// Polar Stereographic projection Parameters
	polarStandardParallel float64 // Latitude of origin in radians
	polarCentralMeridian  float64 // Longitude of origin in radians
	polarFalseEasting     float64 // False easting in meters
	polarFalseNorthing    float64 // False northing in meters

	// Maximum variance for easting and northing values for WGS 84.
	polarDeltaEasting  float64
	polarDeltaNorthing float64

	polarScaleFactor float64
}

// MakePolarStereographic1 receives the ellipsoid parameters and Polar
// Stereograpic (Standard Parallel) projection parameters as inputs, and sets
// the corresponding state variables.
func MakePolarStereographic1(ellipsoidSemiMajorAxis,
	ellipsoidFlattening,
	centralMeridian,
	standardParallel,
	falseEasting,
	falseNorthing float64) (*PolarStereographic, error) {

	p := &PolarStereographic{
		es:                    (0.08181919084262188000),
		esOverTwo:             (.040909595421311),
		isSouthernHemisphere:  (false),
		polarTC:               (1.0),
		polarK90:              (1.0033565552493),
		polaraMc:              (6378137.0),
		twoPolarA:             (12756274.0),
		polarCentralMeridian:  (0.0),
		polarStandardParallel: ((math.Pi * 90) / 180),
		polarFalseEasting:     (0.0),
		polarFalseNorthing:    (0.0),
		polarScaleFactor:      (1.0),
		polarDeltaEasting:     (12713601.0),
		polarDeltaNorthing:    (12713601.0),
	}

	invF := 1 / ellipsoidFlattening
	if ellipsoidSemiMajorAxis <= 0.0 {
		return nil, errors.New("Semi-major axis must be greater than zero")
	}
	if (invF < 250) ||
		(invF > 350) {
		return nil, errors.New("Inverse flattening must be between 250 and 350")
	}
	if (standardParallel < -math.Pi/2) ||
		(standardParallel > math.Pi/2) {
		return nil, errors.New("Origin Latitude out of range")
	}
	if (centralMeridian < -math.Pi) ||
		(centralMeridian > 2*math.Pi) {
		return nil, errors.New("Origin Longitude out of range")
	}

	p.semiMajorAxis = ellipsoidSemiMajorAxis
	p.flattening = ellipsoidFlattening

	p.twoPolarA = 2.0 * p.semiMajorAxis

	if centralMeridian > math.Pi {
		centralMeridian -= 2 * math.Pi
	}
	if standardParallel < 0 {
		p.isSouthernHemisphere = true
		p.polarStandardParallel = -standardParallel
		p.polarCentralMeridian = -centralMeridian
	} else {
		p.isSouthernHemisphere = false
		p.polarStandardParallel = standardParallel
		p.polarCentralMeridian = centralMeridian
	}
	p.polarFalseEasting = falseEasting
	p.polarFalseNorthing = falseNorthing

	es2 := 2*p.flattening - p.flattening*p.flattening
	p.es = math.Sqrt(es2)
	p.esOverTwo = p.es / 2.0

	if math.Abs(math.Abs(p.polarStandardParallel)-math.Pi/2) > 1.0e-10 {
		sinolat := math.Sin(p.polarStandardParallel)
		essin := p.es * sinolat
		powEs := p.polarPow(essin)
		cosolat := math.Cos(p.polarStandardParallel)
		mc := cosolat / math.Sqrt(1.0-essin*essin)
		p.polaraMc = p.semiMajorAxis * mc
		p.polarTC = math.Tan(math.Pi/4-p.polarStandardParallel/2.0) / powEs
	}

	onePlusEs := 1.0 + p.es
	oneMinusEs := 1.0 - p.es
	p.polarK90 = math.Sqrt(math.Pow(onePlusEs, onePlusEs) * math.Pow(oneMinusEs, oneMinusEs))

	slat := math.Sin(math.Abs(standardParallel))
	onePlusEsSinoLat := 1.0 + p.es*slat
	oneMinusEsSinoLat := 1.0 - p.es*slat
	p.polarScaleFactor = ((1 + slat) / 2) *
		(p.polarK90 / math.Sqrt(math.Pow(onePlusEsSinoLat, onePlusEs)*
			math.Pow(oneMinusEsSinoLat, oneMinusEs)))

	// Calculate Radius
	tempGeodeticCoordinates := s2.LatLng{Lng: s1.Angle(centralMeridian), Lat: 0}
	tempCoordinates, err := p.ConvertFromGeodetic(tempGeodeticCoordinates)
	if err != nil {
		return nil, err
	}
	p.polarDeltaNorthing = tempCoordinates.Northing

	if p.polarFalseNorthing != 0 {
		p.polarDeltaNorthing -= p.polarFalseNorthing
	}
	if p.polarDeltaNorthing < 0 {
		p.polarDeltaNorthing = -p.polarDeltaNorthing
	}
	p.polarDeltaNorthing *= 1.01

	p.polarDeltaEasting = p.polarDeltaNorthing

	return p, nil
}

// MakePolarStereographicScaleFactor ellipsoid parameters and Polar Stereograpic
// (Scale Factor) projection parameters as inputs, and sets the corresponding
// state variables.
func MakePolarStereographicScaleFactor(ellipsoidSemiMajorAxis,
	ellipsoidFlattening,
	centralMeridian,
	scaleFactor float64, hemisphere Hemisphere,
	falseEasting,
	falseNorthing float64) (*PolarStereographic, error) {
	p := &PolarStereographic{
		//	 	 coordinateType: (CoordinateType::polarStereographicScaleFactor),
		es:                    (0.08181919084262188000),
		esOverTwo:             (.040909595421311),
		isSouthernHemisphere:  false,
		polarTC:               (1.0),
		polarK90:              (1.0033565552493),
		polaraMc:              (6378137.0),
		twoPolarA:             (12756274.0),
		polarCentralMeridian:  (0.0),
		polarStandardParallel: ((math.Pi * 90) / 180),
		polarFalseEasting:     (0.0),
		polarFalseNorthing:    (0.0),
		polarScaleFactor:      (1.0),
		polarDeltaEasting:     (12713601.0),
		polarDeltaNorthing:    (12713601.0),
	}

	tolerance := 1.0e-15
	count := 30
	invF := 1 / ellipsoidFlattening

	const minScaleFactor = 0.1
	const maxScaleFactor = 3.0

	if ellipsoidSemiMajorAxis <=
		0.0 {
		return nil, errors.New("Semi-major axis must be greater than zero")
	}
	if (invF < 250) ||
		(invF > 350) {
		return nil, errors.New("Inverse flattening must be between 250 and 350")
	}
	if (scaleFactor < minScaleFactor) || (scaleFactor > maxScaleFactor) {
		return nil, errors.New("Scale factor out of range")
	}
	if (centralMeridian < -math.Pi) ||
		(centralMeridian > 2*math.Pi) {
		return nil, errors.New("Origin Longitude out of range")
	}
	if (hemisphere != HemisphereNorth) && (hemisphere != HemisphereSouth) {
		return nil, errors.New("Hemisphere out of range")
	}

	p.semiMajorAxis = ellipsoidSemiMajorAxis
	p.flattening = ellipsoidFlattening
	p.polarScaleFactor = scaleFactor
	p.polarFalseEasting = falseEasting
	p.polarFalseNorthing = falseNorthing

	p.twoPolarA = 2.0 * p.semiMajorAxis
	es2 := 2*p.flattening - p.flattening*p.flattening
	p.es = math.Sqrt(es2)
	p.esOverTwo = p.es / 2.0

	onePlusEs := 1.0 + p.es
	oneMinusEs := 1.0 - p.es
	p.polarK90 =
		math.Sqrt(math.Pow(onePlusEs, onePlusEs) * math.Pow(oneMinusEs, oneMinusEs))

	sk := 0.0
	skPlus1 := -1 + 2*p.polarScaleFactor
	for math.Abs(skPlus1-sk) > tolerance && count != 0 {
		sk = skPlus1
		onePlusEsSk := 1.0 + p.es*sk
		oneMinusEsSk := 1.0 - p.es*sk
		skPlus1 = ((2 * p.polarScaleFactor *
			math.Sqrt(math.Pow(onePlusEsSk, onePlusEs)*
				math.Pow(oneMinusEsSk, oneMinusEs))) /
			p.polarK90) - 1
		count--
	}

	if count == 0 {
		return nil, errors.New("origin latitude error")
	}

	standardParallel := 0.0
	if skPlus1 >= -1.0 && skPlus1 <= 1.0 {
		standardParallel = math.Asin(skPlus1)
	} else {
		return nil, errors.New("origin latitude error")
	}

	if hemisphere == HemisphereSouth {
		standardParallel *= -1.0
	}

	if centralMeridian > math.Pi {
		centralMeridian -= 2 * math.Pi
	}
	if standardParallel < 0 {
		p.isSouthernHemisphere = true
		p.polarStandardParallel = -standardParallel
		p.polarCentralMeridian = -centralMeridian
	} else {
		p.isSouthernHemisphere = false
		p.polarStandardParallel = standardParallel
		p.polarCentralMeridian = centralMeridian
	}

	sinolat := math.Sin(p.polarStandardParallel)

	if math.Abs(math.Abs(p.polarStandardParallel)-math.Pi/2) > 1.0e-10 {
		essin := p.es * sinolat
		powEs := p.polarPow(essin)
		cosolat := math.Cos(p.polarStandardParallel)
		mc := cosolat / math.Sqrt(1.0-essin*essin)
		p.polaraMc = p.semiMajorAxis * mc
		p.polarTC = math.Tan(math.Pi/4-p.polarStandardParallel/2.0) / powEs
	}

	// Calculate Radius
	tempGeodeticCoordinates := s2.LatLng{Lng: s1.Angle(centralMeridian), Lat: 0}
	tempCoordinates, err := p.ConvertFromGeodetic(tempGeodeticCoordinates)
	if err != nil {
		return nil, err
	}
	p.polarDeltaNorthing = tempCoordinates.Northing

	if p.polarFalseNorthing != 0 {
		p.polarDeltaNorthing -= p.polarFalseNorthing
	}
	if p.polarDeltaNorthing < 0 {
		p.polarDeltaNorthing = -p.polarDeltaNorthing
	}
	p.polarDeltaNorthing *= 1.01

	p.polarDeltaEasting = p.polarDeltaNorthing
	return p, nil
}

// ConvertFromGeodetic converts geodetic coordinates (latitude and longitude) to
// Polar Stereographic coordinates (easting and northing), according to the
// current ellipsoid and Polar Stereographic projection parameters.
func (p *PolarStereographic) ConvertFromGeodetic(geodeticCoordinates s2.LatLng) (MapCoords, error) {
	longitude := geodeticCoordinates.Lng.Radians()
	latitude := geodeticCoordinates.Lat.Radians()

	if (latitude < -math.Pi/2) || (latitude > math.Pi/2) {
		return MapCoords{}, errors.New("latitide out of range")
	} else if (latitude < 0) && (!p.isSouthernHemisphere) {
		return MapCoords{}, errors.New("latitude and Origin Latitude in different hemispheres")
	} else if (latitude > 0) && (p.isSouthernHemisphere) {
		return MapCoords{}, errors.New("latitude and Origin Latitude in different hemispheres")
	}
	if (longitude < -math.Pi) || (longitude > 2*math.Pi) {
		return MapCoords{}, errors.New("longitude out of range")
	}

	var easting, northing float64
	if math.Abs(math.Abs(latitude)-math.Pi/2) < 1.0e-10 {
		easting = p.polarFalseEasting
		northing = p.polarFalseNorthing
	} else {
		if p.isSouthernHemisphere {
			longitude *= -1.0
			latitude *= -1.0
		}
		dlam := longitude - p.polarCentralMeridian
		if dlam > math.Pi {
			dlam -= 2 * math.Pi
		}
		if dlam < -math.Pi {
			dlam += 2 * math.Pi
		}
		slat := math.Sin(latitude)
		essin := p.es * slat
		powEs := p.polarPow(essin)
		t := math.Tan(math.Pi/4-latitude/2.0) / powEs

		var rho float64
		if math.Abs(math.Abs(p.polarStandardParallel)-math.Pi/2) > 1.0e-10 {
			rho = p.polaraMc * t / p.polarTC
		} else {
			rho = p.twoPolarA * t / p.polarK90
		}

		if p.isSouthernHemisphere {
			easting = -(rho*math.Sin(dlam) - p.polarFalseEasting)
			northing = rho*math.Cos(dlam) + p.polarFalseNorthing
		} else {
			easting = rho*math.Sin(dlam) + p.polarFalseEasting
			northing = -rho*math.Cos(dlam) + p.polarFalseNorthing
		}
	}
	return MapCoords{Easting: easting, Northing: northing}, nil
}

// ConvertToGeodetic converts Polar Stereographic coordinates (easting and
// northing) to geodetic coordinates (latitude and longitude) according to the
// current ellipsoid and Polar Stereographic projection Parameters.
func (p *PolarStereographic) ConvertToGeodetic(mapProjectionCoordinates MapCoords) (s2.LatLng, error) {
	easting := mapProjectionCoordinates.Easting
	northing := mapProjectionCoordinates.Northing

	minEasting := p.polarFalseEasting - p.polarDeltaEasting
	maxEasting := p.polarFalseEasting + p.polarDeltaEasting
	minNorthing := p.polarFalseNorthing - p.polarDeltaNorthing
	maxNorthing := p.polarFalseNorthing + p.polarDeltaNorthing

	if easting > maxEasting ||
		easting < minEasting {
		return s2.LatLng{}, errors.New("easting out of range")
	}
	if northing > maxNorthing ||
		northing < minNorthing {
		return s2.LatLng{}, errors.New("northing out of range")
	}

	dy := northing - p.polarFalseNorthing
	dx := easting - p.polarFalseEasting

	// Radius of point with origin of false easting, false northing
	rho := math.Sqrt(dx*dx + dy*dy)

	deltaRadius := math.Sqrt(p.polarDeltaEasting*p.polarDeltaEasting +
		p.polarDeltaNorthing*p.polarDeltaNorthing)

	if rho > deltaRadius {
		return s2.LatLng{}, errors.New("Point is outside of projection area")
	}

	var latitude, longitude float64
	if (dy == 0.0) && (dx == 0.0) {
		latitude = math.Pi / 2
		longitude = p.polarCentralMeridian
	} else {
		if p.isSouthernHemisphere {
			dy *= -1.0
			dx *= -1.0
		}

		var t float64
		if math.Abs(math.Abs(p.polarStandardParallel)-math.Pi/2) > 1.0e-10 {
			t = rho * p.polarTC / (p.polaraMc)
		} else {
			t = rho * p.polarK90 / (p.twoPolarA)
		}
		PHI := math.Pi/2 - 2.0*math.Atan(t)
		tempPHI := 0.0
		for math.Abs(PHI-tempPHI) > 1.0e-10 {
			tempPHI = PHI
			sinPhi := math.Sin(PHI)
			essin := p.es * sinPhi
			powEs := p.polarPow(essin)
			PHI = math.Pi/2 - 2.0*math.Atan(t*powEs)
		}
		latitude = PHI
		longitude = p.polarCentralMeridian + math.Atan2(dx, -dy)

		if longitude > math.Pi {
			longitude -= 2 * math.Pi
		} else if longitude < -math.Pi {
			longitude += 2 * math.Pi
		}

		if latitude > math.Pi/2 { // force distorted values to 90, -90 degrees
			latitude = math.Pi / 2
		} else if latitude < -math.Pi/2 {
			latitude = -math.Pi / 2
		}

		if longitude > math.Pi { // force distorted values to 180, -180 degrees
			longitude = math.Pi
		} else if longitude < -math.Pi {
			longitude = -math.Pi
		}
	}
	if p.isSouthernHemisphere {
		latitude *= -1.0
		longitude *= -1.0
	}

	return s2.LatLng{Lat: s1.Angle(latitude), Lng: s1.Angle(longitude)}, nil
}

func (p *PolarStereographic) polarPow(esSin float64) float64 {
	return math.Pow((1.0-esSin)/(1.0+esSin), p.esOverTwo)
}
