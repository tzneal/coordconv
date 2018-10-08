package coordconv

import (
	"bytes"
	"errors"
	"fmt"
	"math"
	"unicode"

	"github.com/golang/geo/s1"
	"github.com/golang/geo/s2"
)

// MGRS is an coordinate converter to and from MGRS coordinates.
type MGRS struct {
	semiMajorAxis     float64
	flattening        float64
	ups               *UPS
	utm               *UTM
	MGRSEllipsoidCode string
}

const espilon2 = 4.99e-4
const mgrsMaxPrecision = 5                             // Maximum precision of easting & northing
const minMGRSNonPolarLat = (-80.0 * (math.Pi / 180.0)) // -80 deg in rad
const maxMGRSNonPolarLat = (84.0 * (math.Pi / 180.0))  //  84 deg in rad
const mgrsMinEasting = 100000.0
const mgrsMaxEasting = 900000.0
const mgrsMinNorthing = 0.0
const mgrsMaxNorthing = 10000000.0

const letterA = 0
const letterB = 1
const letterC = 2
const letterD = 3
const letterE = 4
const letterF = 5
const letterG = 6
const letterH = 7
const letterI = 8
const letterJ = 9
const letterK = 10
const letterL = 11
const letterM = 12
const letterN = 13
const letterO = 14
const letterP = 15
const letterQ = 16
const letterR = 17
const letterS = 18
const letterT = 19
const letterU = 20
const letterV = 21
const letterW = 22
const letterX = 23
const letterY = 24
const letterZ = 25

const clarke1886 = "CC"
const clarke1880 = "CD"
const bessel1841 = "BR"
const bessel1841Namibia = "BN"

type latitudeBand struct {
	letter         int     // letter representing latitude band
	minNorthing    float64 // minimum northing for latitude band
	north          float64 // upper latitude for latitude band
	south          float64 // lower latitude for latitude band
	northingOffset float64 // latitude band northing offset
}

var latitudeBands = [20]latitudeBand{
	{letterC, 1100000.0, -72.0, -80.5, 0.0},
	{letterD, 2000000.0, -64.0, -72.0, 2000000.0},
	{letterE, 2800000.0, -56.0, -64.0, 2000000.0},
	{letterF, 3700000.0, -48.0, -56.0, 2000000.0},
	{letterG, 4600000.0, -40.0, -48.0, 4000000.0},
	{letterH, 5500000.0, -32.0, -40.0, 4000000.0},
	{letterJ, 6400000.0, -24.0, -32.0, 6000000.0},
	{letterK, 7300000.0, -16.0, -24.0, 6000000.0},
	{letterL, 8200000.0, -8.0, -16.0, 8000000.0},
	{letterM, 9100000.0, 0.0, -8.0, 8000000.0},
	{letterN, 0.0, 8.0, 0.0, 0.0},
	{letterP, 800000.0, 16.0, 8.0, 0.0},
	{letterQ, 1700000.0, 24.0, 16.0, 0.0},
	{letterR, 2600000.0, 32.0, 24.0, 2000000.0},
	{letterS, 3500000.0, 40.0, 32.0, 2000000.0},
	{letterT, 4400000.0, 48.0, 40.0, 4000000.0},
	{letterU, 5300000.0, 56.0, 48.0, 4000000.0},
	{letterV, 6200000.0, 64.0, 56.0, 6000000.0},
	{letterW, 7000000.0, 72.0, 64.0, 6000000.0},
	{letterX, 7900000.0, 84.5, 72.0, 6000000.0}}

type upsConstant struct {
	letter        int     // letter representing latitude band
	ltr2LowValue  int     // 2nd letter range - low number
	ltr2HighValue int     // 2nd letter range - high number
	ltr3HighValue int     // 3rd letter range - high number (UPS)
	falseEasting  float64 // False easting based on 2nd letter
	falseNorthing float64 // False northing based on 3rd letter
}

var upsConstants = [4]upsConstant{
	{letterA, letterJ, letterZ, letterZ, 800000.0, 800000.0},
	{letterB, letterA, letterR, letterZ, 2000000.0, 800000.0},
	{letterY, letterJ, letterZ, letterP, 800000.0, 1300000.0},
	{letterZ, letterA, letterJ, letterP, 2000000.0, 1300000.0}}

// NewMGRS constructs an MGRS converter with ellipsoid parameters.
func NewMGRS(ellipsoidSemiMajorAxis, ellipsoidFlattening float64,
	ellipsoidCode string) (*MGRS, error) {
	m := &MGRS{}

	invF := 1 / ellipsoidFlattening
	if ellipsoidSemiMajorAxis <=
		0.0 {
		return nil, errors.New("Semi-major axis must be greater than zero")
	}
	if (invF < 250) || (invF > 350) {
		return nil, errors.New("Inverse flattening must be between 250 and 350")
	}

	m.semiMajorAxis = ellipsoidSemiMajorAxis
	m.flattening = ellipsoidFlattening
	m.MGRSEllipsoidCode = ellipsoidCode

	var err error
	m.ups, err = NewUPS(m.semiMajorAxis, m.flattening)
	if err != nil {
		return nil, err
	}

	m.utm, err = NewUTM2(m.semiMajorAxis, m.flattening, m.MGRSEllipsoidCode, 0)
	if err != nil {
		return nil, err
	}
	return m, nil
}

// ConvertFromGeodetic converts Geodetic (latitude and longitude) coordinates to
// an MGRS coordinate string, according to the current ellipsoid parameters.
func (m *MGRS) ConvertFromGeodetic(geodeticCoordinates s2.LatLng, precision int) (string, error) {
	latitude := geodeticCoordinates.Lat.Radians()
	longitude := geodeticCoordinates.Lng.Radians()

	if (latitude < -math.Pi/2) ||
		(latitude > math.Pi/2) {
		return "", errors.New("latitude out of range")
	}
	if (longitude < (-math.Pi - epsilonRadians)) ||
		(longitude > (2*math.Pi + epsilonRadians)) {
		return "", errors.New("longitude out of range")
	}
	if (precision < 0) || (precision > mgrsMaxPrecision) {
		return "", errors.New("precision out of range")
	}

	var mgrsCoords string
	// If the latitude is within the valid mgrs non polar range [-80, 84),
	// convert to mgrs using the utm path,
	// otherwise convert to mgrs using the ups path
	if (latitude >= minMGRSNonPolarLat-epsilonRadians) &&
		(latitude < maxMGRSNonPolarLat+epsilonRadians) {
		utmCoordinates, err := m.utm.ConvertFromGeodetic(geodeticCoordinates, 0)
		if err != nil {
			return "", err
		}
		mgrsCoords, err = m.fromUTM(utmCoordinates, longitude, latitude, precision)
		if err != nil {
			return "", err
		}
	} else {

		upsCoordinates, err := m.ups.ConvertFromGeodetic(geodeticCoordinates)
		if err != nil {
			return "", err
		}
		mgrsCoords, err = m.fromUPS(upsCoordinates, precision)
		if err != nil {
			return "", err
		}
	}
	return mgrsCoords, nil
}

// ConvertFromUTM converts UTM (zone, easting, and northing) coordinates to an
// MGRS coordinate string, according to the current ellipsoid parameters.  If
// any errors occur, an exception is thrown with a description of the error.
func (m *MGRS) ConvertFromUTM(utmCoordinates UTMCoord, precision int) (string, error) {
	zone := utmCoordinates.Zone
	hemisphere := utmCoordinates.Hemisphere
	easting := utmCoordinates.Easting
	northing := utmCoordinates.Northing

	if (zone < 1) || (zone > 60) {
		return "", errors.New("zone out of range")
	}
	if (hemisphere != HemisphereSouth) && (hemisphere != HemisphereNorth) {
		return "", errors.New("hemisphere out of range")
	}
	if (easting < mgrsMinEasting) || (easting > mgrsMaxEasting) {
		return "", errors.New("esating out of range")
	}
	if (northing < mgrsMinNorthing) || (northing > mgrsMaxNorthing) {
		return "", errors.New("northing out of range")
	}
	if (precision < 0) || (precision > mgrsMaxPrecision) {
		return "", errors.New("precision out of range")
	}

	geodeticCoordinates, err := m.utm.ConvertToGeodetic(utmCoordinates)
	if err != nil {
		return "", err
	}

	// If the latitude is within the valid mgrs non polar range [-80, 84),
	// convert to mgrs using the utm path,
	// otherwise convert to mgrs using the ups path
	latitude := geodeticCoordinates.Lat.Radians()

	var mgrsCoord string
	if (latitude >= (minMGRSNonPolarLat - epsilonRadians)) &&
		(latitude < (maxMGRSNonPolarLat + epsilonRadians)) {
		mgrsCoord, err = m.fromUTM(utmCoordinates, geodeticCoordinates.Lng.Radians(),
			latitude, precision)
		if err != nil {
			return "", err
		}
	} else {
		upsCoordinates, err := m.ups.ConvertFromGeodetic(geodeticCoordinates)
		if err != nil {
			return "", err
		}
		mgrsCoord, err = m.fromUPS(upsCoordinates, precision)
		if err != nil {
			return "", err
		}
	}

	return mgrsCoord, nil
}

// fromUPS converts UPS (hemisphere, easting, and northing) coordinates to an
// MGRS coordinate string according to the current ellipsoid parameters.
func (m *MGRS) fromUPS(upsCoordinates UPSCoord, precision int) (string, error) {
	hemisphere := upsCoordinates.Hemisphere
	easting := upsCoordinates.Easting
	northing := upsCoordinates.Northing

	divisor := computeScale(precision)

	easting = float64(int((easting+espilon2)/divisor)) * divisor
	northing = float64(int((northing+espilon2)/divisor)) * divisor

	var letters [3]byte
	var falseEasting float64  // False easting for 2nd letter
	var falseNorthing float64 // False northing for 3rd letter
	var ltr2LowValue int      // 2nd letter range - low number
	if hemisphere == HemisphereNorth {
		if easting >= 2000000 {
			letters[0] = letterZ
		} else {
			letters[0] = letterY
		}

		index := letters[0] - 22
		ltr2LowValue = upsConstants[index].ltr2LowValue
		falseEasting = upsConstants[index].falseEasting
		falseNorthing = upsConstants[index].falseNorthing
	} else {
		if easting >= 2000000 {
			letters[0] = letterB
		} else {
			letters[0] = letterA
		}

		ltr2LowValue = upsConstants[letters[0]].ltr2LowValue
		falseEasting = upsConstants[letters[0]].falseEasting
		falseNorthing = upsConstants[letters[0]].falseNorthing
	}

	gridNorthing := northing
	gridNorthing = gridNorthing - falseNorthing
	letters[2] = byte(gridNorthing / 100000)

	if letters[2] > letterH {
		letters[2] = letters[2] + 1
	}

	if letters[2] > letterN {
		letters[2] = letters[2] + 1
	}

	gridEasting := easting
	gridEasting = gridEasting - falseEasting
	letters[1] = byte(ltr2LowValue + (int(gridEasting / 100000)))

	if easting < 2000000 {
		if letters[1] > letterL {
			letters[1] = letters[1] + 3
		}
		if letters[1] > letterU {
			letters[1] = letters[1] + 2
		}
	} else {
		if letters[1] > letterC {
			letters[1] = letters[1] + 2
		}
		if letters[1] > letterH {
			letters[1] = letters[1] + 1
		}
		if letters[1] > letterL {
			letters[1] = letters[1] + 3
		}
	}

	MGRSString, err := makeMGRSString(0, letters[:], easting, northing, precision)
	if err != nil {
		return "", err
	}

	return MGRSString, nil
}

// fromUTM calculates an MGRS coordinate string based on the zone, latitude,
// easting and northing.
func (m *MGRS) fromUTM(utmCoordinates UTMCoord, longitude, latitude float64, precision int) (string, error) {
	var letters [3]byte
	zone := utmCoordinates.Zone
	easting := utmCoordinates.Easting
	northing := utmCoordinates.Northing

	var err error
	letters[0], err = getLatitudeLetter(latitude)
	if err != nil {
		return "", nil
	}

	const Lat6 = (6.0 * (math.Pi / 180.0))
	// Check if the point is within it's natural zone
	// If it is not, put it there
	pad := espilon2 / 6378137.0
	var naturalZone int
	if longitude < math.Pi {
		naturalZone = int(31 + ((longitude + pad) / Lat6))
	} else {
		naturalZone = int(((longitude + pad) / Lat6) - 29)
	}

	if naturalZone > 60 {
		naturalZone = 1
	}
	if zone != naturalZone { // reconvert to override zone
		utmOverride, err := NewUTM2(m.semiMajorAxis, m.flattening, m.MGRSEllipsoidCode, naturalZone)
		if err != nil {
			return "", err
		}
		geodeticCoordinates := s2.LatLng{Lng: s1.Angle(longitude), Lat: s1.Angle(latitude)}
		utmCoordinatesOverride, err := utmOverride.ConvertFromGeodetic(geodeticCoordinates, 0)
		if err != nil {
			return "", err
		}
		zone = utmCoordinatesOverride.Zone
		easting = utmCoordinatesOverride.Easting
		northing = utmCoordinatesOverride.Northing
	}

	// UTM special cases
	var override int
	if letters[0] == letterV {
		// V latitude band
		if (zone == 31) && (easting >= 500000.0) {
			override = 32 // extension of zone 32V
		}
	} else if letters[0] == letterX {
		if (zone == 32) && (easting < 500000.0) { // extension of zone 31X
			override = 31
		} else if ((zone == 32) &&
			(easting >= 500000.0)) || // western extension of zone 33X
			((zone == 34) &&
				(easting < 500000.0)) { // eastern extension of zone 33X
			override = 33
		} else if ((zone == 34) &&
			(easting >= 500000.0)) || // western extension of zone 35X
			((zone == 36) &&
				(easting < 500000.0)) { // eastern extension of zone 35X
			override = 35
		} else if (zone == 36) &&
			(easting >= 500000.0) { // western extension of zone 37X
			override = 37
		}
	}

	if override != 0 { // reconvert to override zone
		utmOverride, err := NewUTM2(m.semiMajorAxis, m.flattening, m.MGRSEllipsoidCode, override)
		if err != nil {
			return "", err
		}
		geodeticCoordinates := s2.LatLng{Lng: s1.Angle(longitude), Lat: s1.Angle(latitude)}
		utmCoordinatesOverride, err := utmOverride.ConvertFromGeodetic(geodeticCoordinates, 0)
		if err != nil {
			return "", err
		}

		zone = utmCoordinatesOverride.Zone
		easting = utmCoordinatesOverride.Easting
		northing = utmCoordinatesOverride.Northing
	}

	divisor := computeScale(precision)

	easting = float64(int((easting+espilon2)/divisor)) * divisor
	northing = float64(int((northing+espilon2)/divisor)) * divisor

	if latitude <= 0.0 && northing == 1.0e7 {
		latitude = 0.0
		northing = 0.0
	}

	ltr2LowValue, _, patternOffset := m.getGridValues(zone)

	// Northing used to derive 3rd letter of MGRS
	gridNorthing := northing

	for gridNorthing >= 2000000 {
		gridNorthing = gridNorthing - 2000000
	}
	gridNorthing = gridNorthing + patternOffset
	if gridNorthing >= 2000000 {
		gridNorthing = gridNorthing - 2000000
	}

	letters[2] = byte(gridNorthing / 100000)
	if letters[2] > letterH {
		letters[2] = letters[2] + 1
	}

	if letters[2] > letterN {
		letters[2] = letters[2] + 1
	}

	letters[1] = byte(ltr2LowValue + (int(easting/100000) - 1))
	if (ltr2LowValue == letterJ) && (letters[1] > letterN) {
		letters[1] = letters[1] + 1
	}

	MGRSString, err := makeMGRSString(zone, letters[:], easting, northing, precision)
	if err != nil {
		return "", err
	}
	return MGRSString, nil
}

func computeScale(prec int) float64 {
	scale := 1.0e5
	switch prec {
	case 0:
		scale = 1.0e5
	case 1:
		scale = 1.0e4
	case 2:
		scale = 1.0e3
	case 3:
		scale = 1.0e2
	case 4:
		scale = 1.0e1
	case 5:
		scale = 1.0e0
	}
	return scale
}

// makeMGRSString constructs an MGRS string from its component parts
func makeMGRSString(zone int, letters []byte,
	easting, northing float64, precision int) (string, error) {

	const alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

	buf := bytes.Buffer{}
	if zone != 0 {
		fmt.Fprintf(&buf, "%2.2d", zone)
	}

	for j := 0; j < 3; j++ {
		if letters[j] < 0 || letters[j] >= byte(len(alphabet)) {
			return "", errors.New("invalid letters")
		}
		buf.WriteByte(alphabet[letters[j]])
	}

	divisor := computeScale(precision)

	easting = math.Mod(easting, 100000.0)
	if easting >= 99999.5 {
		easting = 99999.0
	}
	const espilon2 = 4.99e-1
	east := int((easting + espilon2) / divisor)
	fmt.Fprintf(&buf, "%*.*d", precision, precision, east)
	northing = math.Mod(northing, 100000.0)
	if northing >= 99999.5 {
		northing = 99999.0
	}
	north := int((northing + espilon2) / divisor)
	fmt.Fprintf(&buf, "%*.*d", precision, precision, north)
	return buf.String(), nil
}

// breakMGRSString breaks down an MGRS coordinate string into its component
// parts.
func breakMGRSString(MGRSString string) (zone int, letters []byte,
	easting, northing float64, precision int, err error) {

	// remove any spaces from MGRS string
	buf := bytes.Buffer{}
	for _, b := range MGRSString {
		// check for invalid character
		if !isdigit(byte(b)) && !isalpha(byte(b)) {
			err = errors.New("invalid character")
			return
		}
		buf.WriteRune(b)
	}
	tempMGRSString := buf.String()
	i := 0
	for isdigit(tempMGRSString[i]) {
		i++
	}
	numDigits := i
	if numDigits <= 2 {
		if numDigits > 0 {
			// get zone
			zoneString := tempMGRSString[:2]
			fmt.Sscanf(zoneString, "%d", &zone)
			if (zone < 1) || (zone > 60) {
				err = errors.New("zone out of range")
				return
			}
		} else {
			zone = 0
		}
	} else {
		err = errors.New("too few digits")
		return
	}
	j := i
	for i < len(tempMGRSString) && isalpha(tempMGRSString[i]) {
		i++
	}
	numLetters := i - j
	if numLetters == 3 {
		letters = make([]byte, 3)
		// get letters
		letters[0] = (toupper(tempMGRSString[j]) - 'A')
		if (letters[0] == letterI) || (letters[0] == letterO) {
			err = errors.New("invalid letter 0")
			return
		}
		letters[1] = (toupper(tempMGRSString[j+1]) - 'A')
		if (letters[1] == letterI) || (letters[1] == letterO) {
			err = errors.New("invalid letter 1")
			return
		}
		letters[2] = (toupper(tempMGRSString[j+2]) - 'A')
		if (letters[2] == letterI) || (letters[2] == letterO) {
			err = errors.New("invalid letter 2")
			return
		}
	} else {
		err = errors.New("wrong number of letters")
		return
	}
	j = i
	for i < len(tempMGRSString) && isdigit(tempMGRSString[i]) {
		i++
	}
	numDigits = i - j
	if (numDigits <= 10) && (numDigits%2 == 0) {
		// get easting & northing
		var east, north int
		n := numDigits / 2
		precision = n
		if n > 0 {
			eastString := tempMGRSString[j : j+n]
			fmt.Sscanf(eastString, "%d", &east)
			northString := tempMGRSString[j+n : j+2*n]
			fmt.Sscanf(northString, "%d", &north)
			multiplier := computeScale(n)

			easting = float64(east) * multiplier   // + (multiplier/2.0); // move to
			northing = float64(north) * multiplier // + (multiplier/2.0); // grid center
		} else {
			easting = 0.0
			northing = 0.0
		}
	} else {
		err = errors.New("wrong number of digits")
		return
	}
	return
}

// getGridValues sets the letter range used for the 2nd letter in the MGRS
// coordinate string, based on the set number of the utm zone. It also sets the
// pattern offset using a value of A for the second letter of the grid square,
// based on the grid pattern and set number of the utm zone.
func (m *MGRS) getGridValues(zone int) (ltr2LowValue, ltr2HighValue int, patternOffset float64) {
	// Set number (1-6) based on UTM zone number
	setNumber := zone % 6
	if setNumber == 0 {
		setNumber = 6
	}

	// Pattern based on ellipsoid code
	var aaPattern bool
	if m.MGRSEllipsoidCode == clarke1886 ||
		m.MGRSEllipsoidCode == clarke1880 ||
		m.MGRSEllipsoidCode == bessel1841 ||
		m.MGRSEllipsoidCode == bessel1841Namibia {
		aaPattern = false
	} else {
		aaPattern = true
	}

	if (setNumber == 1) || (setNumber == 4) {
		ltr2LowValue = letterA
		ltr2HighValue = letterH
	} else if (setNumber == 2) || (setNumber == 5) {
		ltr2LowValue = letterJ
		ltr2HighValue = letterR
	} else if (setNumber == 3) || (setNumber == 6) {
		ltr2LowValue = letterS
		ltr2HighValue = letterZ
	}

	// False northing at A for second letter of grid square
	if aaPattern {
		if (setNumber % 2) == 0 {
			patternOffset = 500000.0
		} else {
			patternOffset = 0.0
		}
	} else {
		if (setNumber % 2) == 0 {
			patternOffset = 1500000.0
		} else {
			patternOffset = 1000000.00
		}
	}
	return
}

// ConvertToGeodetic converts an MGRS coordinate string to Geodetic (latitude
// and longitude) coordinates according to the current ellipsoid parameters.
func (m *MGRS) ConvertToGeodetic(mgrsorUSNGCoordinates string) (s2.LatLng, error) {
	zone, letters, mgrsEasting, mgrsNorthing, precision, err := breakMGRSString(mgrsorUSNGCoordinates)
	if err != nil {
		return s2.LatLng{}, err
	}
	var geodeticCoordinates s2.LatLng
	if zone != 0 {
		utmCoordinates, err := m.toUTM(zone, letters, mgrsEasting, mgrsNorthing, precision)
		if err != nil {
			return s2.LatLng{}, err
		}

		geodeticCoordinates, err = m.utm.ConvertToGeodetic(utmCoordinates)
		if err != nil {
			return s2.LatLng{}, err
		}
	} else {
		upsCoordinates, err := m.toUPS(letters, mgrsEasting, mgrsNorthing)
		if err != nil {
			return s2.LatLng{}, err
		}
		geodeticCoordinates, err = m.ups.ConvertToGeodetic(upsCoordinates)
		if err != nil {
			return s2.LatLng{}, err
		}
	}
	return geodeticCoordinates, err
}

// toUTM converts an MGRS coordinate string to UTM projection (zone, hemisphere,
// easting and northing) coordinates according to the current ellipsoid
// parameters. The string return value is a possible warning about conditions
// that occurred during conversion.
func (m *MGRS) toUTM(zone int, letters []byte, easting, northing float64, precision int) (UTMCoord, error) {

	var hemisphere Hemisphere
	if (letters[0] == letterX) && ((zone == 32) || (zone == 34) || (zone == 36)) {
		return UTMCoord{}, errors.New("invalid letters")
	} else if (letters[0] == letterV) && (zone == 31) && (letters[1] > letterD) {
		return UTMCoord{}, errors.New("invalid letters")
	}

	if letters[0] < letterN {
		hemisphere = HemisphereSouth
	} else {
		hemisphere = HemisphereNorth
	}

	ltr2LowValue, ltr2HighValue, patternOffset := m.getGridValues(zone)

	// Check that the second letter of the MGRS string is within the range of
	// valid second letter values. Also check that the third letter is valid
	if (letters[1] < byte(ltr2LowValue)) ||
		(letters[1] > byte(ltr2HighValue)) ||
		(letters[2] > letterV) {
		return UTMCoord{}, errors.New("invalid letters")
	}

	gridEasting := float64((letters[1])-byte(ltr2LowValue)+1) * 100000
	if (ltr2LowValue == letterJ) && (letters[1] > letterO) {
		gridEasting = gridEasting - 100000
	}

	rowLetterNorthing := (float64)(letters[2]) * 100000
	if letters[2] > letterO {
		rowLetterNorthing = rowLetterNorthing - 100000
	}

	if letters[2] > letterI {
		rowLetterNorthing = rowLetterNorthing - 100000
	}

	if rowLetterNorthing >= 2000000 {
		rowLetterNorthing = rowLetterNorthing - 2000000
	}

	minNorthing, northingOffset, err := m.getLatitudeBandMinNorthing(letters[0])
	if err != nil {
		return UTMCoord{}, err
	}

	gridNorthing := rowLetterNorthing - patternOffset
	if gridNorthing < 0 {
		gridNorthing += 2000000
	}

	gridNorthing += northingOffset

	if gridNorthing < minNorthing {
		gridNorthing += 2000000
	}

	easting = gridEasting + easting
	northing = gridNorthing + northing

	utmCoordinates := UTMCoord{
		Zone:       zone,
		Hemisphere: hemisphere,
		Easting:    easting,
		Northing:   northing,
	}

	// check that point is within Zone Letter bounds
	geodeticCoordinates, err := m.utm.ConvertToGeodetic(utmCoordinates)
	if err != nil {
		return UTMCoord{}, err
	}
	latitude := geodeticCoordinates.Lat.Radians()

	divisor := 100000 / computeScale(precision)

	inRange, err := m.inLatitudeRange(letters[0], latitude, math.Pi/180/divisor)
	if err != nil {
		return UTMCoord{}, err
	}

	if !inRange {
		// check adjacent bands
		prevBand := letters[0] - 1
		nextBand := letters[0] + 1

		if letters[0] == letterC { // if last band, do not go off list
			prevBand = letters[0]
		}

		if letters[0] == letterX {
			nextBand = letters[0]
		}

		if prevBand == letterI || prevBand == letterO {
			prevBand--
		}

		if nextBand == letterI || nextBand == letterO {
			nextBand++
		}

		prevInRange, err := m.inLatitudeRange(prevBand, latitude, math.Pi/180/divisor)
		if err != nil {
			return UTMCoord{}, err
		}
		nextInRange, err := m.inLatitudeRange(nextBand, latitude, math.Pi/180/divisor)
		if err != nil {
			return UTMCoord{}, err
		}
		if prevInRange && nextInRange {
			// TODO: do we care?
			// warning - "Latitude band boundary cuts across 100km square"
		} else {
			return UTMCoord{}, errors.New("MGRS invalid")
		}
	}
	return utmCoordinates, nil
}

// getLatitudeBandMinNorthing receives a latitude band letter and uses the
// latitudeBands to determine the minimum northing and northing offset for
// that latitude band letter.
func (m *MGRS) getLatitudeBandMinNorthing(letter byte) (minNorthing, northingOffset float64, err error) {
	if (letter >= letterC) && (letter <= letterH) {
		minNorthing = latitudeBands[letter-2].minNorthing
		northingOffset = latitudeBands[letter-2].northingOffset
	} else if (letter >= letterJ) && (letter <= letterN) {
		minNorthing = latitudeBands[letter-3].minNorthing
		northingOffset = latitudeBands[letter-3].northingOffset
	} else if (letter >= letterP) && (letter <= letterX) {
		minNorthing = latitudeBands[letter-4].minNorthing
		northingOffset = latitudeBands[letter-4].northingOffset
	} else {
		err = errors.New("invalid MGRS")
	}
	return
}

// toUPS converts an MGRS coordinate string to UPS (hemisphere, easting, and northing)
// coordinates, according to the current ellipsoid parameters.
func (m *MGRS) toUPS(
	letters []byte,
	easting,
	northing float64) (UPSCoord, error) {
	var ltr2HighValue int     // 2nd letter range - high number
	var ltr3HighValue int     // 3rd letter range - high number (UPS)
	var ltr2LowValue int      // 2nd letter range - low number
	var falseEasting float64  // False easting for 2nd letter
	var falseNorthing float64 // False northing for 3rd letter
	var gridEasting float64   // easting for 100,000 meter grid square
	var gridNorthing float64  // northing for 100,000 meter grid square
	var hemisphere Hemisphere
	var index int

	if (letters[0] == letterY) || (letters[0] == letterZ) {
		hemisphere = HemisphereNorth

		index = int(letters[0]) - 22
		ltr2LowValue = upsConstants[index].ltr2LowValue
		ltr2HighValue = upsConstants[index].ltr2HighValue
		ltr3HighValue = upsConstants[index].ltr3HighValue
		falseEasting = upsConstants[index].falseEasting
		falseNorthing = upsConstants[index].falseNorthing
	} else if (letters[0] == letterA) || (letters[0] == letterB) {
		hemisphere = HemisphereSouth

		ltr2LowValue = upsConstants[letters[0]].ltr2LowValue
		ltr2HighValue = upsConstants[letters[0]].ltr2HighValue
		ltr3HighValue = upsConstants[letters[0]].ltr3HighValue
		falseEasting = upsConstants[letters[0]].falseEasting
		falseNorthing = upsConstants[letters[0]].falseNorthing
	} else {
		return UPSCoord{}, errors.New("Invalid MGRS string")
	}

	// Check that the second letter of the MGRS string is within the range of
	// valid second letter values Also check that the third letter is valid
	if (int(letters[1]) < ltr2LowValue) || (int(letters[1]) > ltr2HighValue) ||
		((letters[1] == letterD) || (letters[1] == letterE) ||
			(letters[1] == letterM) || (letters[1] == letterN) ||
			(letters[1] == letterV) || (letters[1] == letterW)) ||
		(int(letters[2]) > ltr3HighValue) {
		return UPSCoord{}, errors.New("Invalid MGRS string")
	}

	gridNorthing = float64(letters[2])*100000 + falseNorthing
	if letters[2] > letterI {
		gridNorthing = gridNorthing - 100000
	}

	if letters[2] > letterO {
		gridNorthing = gridNorthing - 100000
	}

	gridEasting = float64(int(letters[1])-ltr2LowValue)*100000 + falseEasting
	if ltr2LowValue != letterA {
		if letters[1] > letterL {
			gridEasting = gridEasting - 300000.0
		}

		if letters[1] > letterU {
			gridEasting = gridEasting - 200000.0
		}
	} else {
		if letters[1] > letterC {
			gridEasting = gridEasting - 200000.0
		}

		if letters[1] > letterI {
			gridEasting = gridEasting - 100000
		}

		if letters[1] > letterL {
			gridEasting = gridEasting - 300000.0
		}
	}

	easting = gridEasting + easting
	northing = gridNorthing + northing
	return UPSCoord{
		Hemisphere: hemisphere,
		Easting:    easting,
		Northing:   northing,
	}, nil
}

// inLatitudeRange receives a latitude band letter and uses the
// latitudeBands to determine the latitude band boundaries for that
// latitude band letter.
func (m *MGRS) inLatitudeRange(letter byte, latitude, border float64) (bool, error) {
	var north, south float64
	if (letter >= letterC) && (letter <= letterH) {
		north = latitudeBands[letter-2].north * math.Pi / 180
		south = latitudeBands[letter-2].south * math.Pi / 180
	} else if (letter >= letterJ) && (letter <= letterN) {
		north = latitudeBands[letter-3].north * math.Pi / 180
		south = latitudeBands[letter-3].south * math.Pi / 180
	} else if (letter >= letterP) && (letter <= letterX) {
		north = latitudeBands[letter-4].north * math.Pi / 180
		south = latitudeBands[letter-4].south * math.Pi / 180
	} else {
		return false, errors.New("invalid MGRS")
	}

	if ((south - border) <= latitude) && (latitude <= (north + border)) {
		return true, nil
	}
	return false, nil
}
func isdigit(r byte) bool {
	return r >= '0' && r <= '9'
}

func isalpha(r byte) bool {
	return r >= 'a' && r <= 'z' ||
		r >= 'A' && r <= 'Z'
}
func toupper(b byte) byte {
	return byte(unicode.ToUpper(rune(b)))
}

// getLatitudeLetter receives a latitude value and uses the latitudeBands
// to determine the latitude band letter for that latitude.
func getLatitudeLetter(latitude float64) (byte, error) {
	const Lat72 = (72.0 * (math.Pi / 180.0))
	const Lat845 = (84.5 * (math.Pi / 180.0))
	const Lat80 = (80.0 * (math.Pi / 180.0))
	const Lat805 = (80.5 * (math.Pi / 180.0))
	const Lat8 = (8.0 * (math.Pi / 180.0))

	if latitude >= Lat72 && latitude < Lat845 {
		return letterX, nil
	} else if latitude > -Lat805 && latitude < Lat72 {
		band := int(((latitude + Lat80) / Lat8) + 1.0e-12)
		if band < 0 {
			band = 0
		}
		return byte(latitudeBands[band].letter), nil
	}
	return 0, errors.New("latitude out of range")
}
