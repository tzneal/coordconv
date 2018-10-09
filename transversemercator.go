package coordconv

import (
	"errors"
	"math"

	"github.com/golang/geo/s1"

	"github.com/golang/geo/s2"
)

const nTerms = 6

// TransverseMercator provides conversions between Geodetic coordinates
// (latitude and longitude) and Transverse Mercator projection coordinates
// (easting and northing).
type TransverseMercator struct {
	// Ellipsoid Parameters
	semiMajorAxis float64
	flattening    float64
	ellipsCode    string // 2 Letter ellipsoid code

	tranMercEps float64 // Eccentricity

	tranMercK0R4    float64 // SCALE_FACTOR*R4
	tranMercK0R4inv float64 // 1/(SCALE_FACTOR*R4)

	tranMercACoeff [8]float64
	tranMercBCoeff [8]float64

	// Transverse_Mercator projection Parameters
	tranMercOriginLat     float64 // Latitude of origin in radians
	tranMercOriginLong    float64 // Longitude of origin in radians
	TranMercFalseNorthing float64 // False northing in meters
	tranMercFalseEasting  float64 // False easting in meters
	tranMercScaleFactor   float64 // Scale factor

	// Maximum variance for easting and northing values
	tranMercDeltaEasting  float64
	tranMercDeltaNorthing float64
}

// NewTransverseMercator constructs a new TransverseMercator converter.
func NewTransverseMercator(ellipsoidSemiMajorAxis, ellipsoidFlattening, centralMeridian,
	latitudeOfTrueScale, falseEasting, falseNorthing, scaleFactor float64,
	ellipsoidCode string) (*TransverseMercator, error) {
	invFlattening := 1.0 / ellipsoidFlattening

	t := &TransverseMercator{
		tranMercOriginLong:    centralMeridian,
		tranMercOriginLat:     latitudeOfTrueScale,
		tranMercFalseEasting:  falseEasting,
		TranMercFalseNorthing: falseNorthing,
		tranMercScaleFactor:   scaleFactor,
		tranMercDeltaEasting:  20000000.0,
		tranMercDeltaNorthing: 10000000.0,
	}
	t.semiMajorAxis = ellipsoidSemiMajorAxis
	t.flattening = ellipsoidFlattening

	if ellipsoidCode == "" {
		return nil, errors.New("missing ellipsoid code")
	}
	if ellipsoidSemiMajorAxis <=
		0.0 {
		return nil, errors.New("Semi-major axis must be greater than zero")
	}
	if invFlattening < 150 {
		return nil, errors.New("inverse ellipsoid flattening out of range")
	}
	if (latitudeOfTrueScale < -math.Pi/2) ||
		(latitudeOfTrueScale > math.Pi/2) {
		return nil, errors.New("latitudeOfTrueScale out of range")
	}
	if (centralMeridian < -math.Pi) ||
		(centralMeridian > (2 * math.Pi)) {
		return nil, errors.New("centralMeridian out of range")
	}

	const minScaleFactor = 0.1
	const maxScaleFactor = 10.0
	if (scaleFactor < minScaleFactor) || (scaleFactor > maxScaleFactor) {
		return nil, errors.New("scale factor out of range")
	}

	if t.tranMercOriginLong > math.Pi {
		t.tranMercOriginLong -= (2 * math.Pi)
	}

	// Eccentricity
	t.tranMercEps = math.Sqrt(2*t.flattening - t.flattening*t.flattening)

	var n1, R4oa float64
	// added ellipsoid code as part of DR30125
	t.generateCoefficients(invFlattening, &n1, t.tranMercACoeff[:], t.tranMercBCoeff[:],
		&R4oa, t.ellipsCode)

	t.tranMercK0R4 = R4oa * t.tranMercScaleFactor * ellipsoidSemiMajorAxis
	t.tranMercK0R4inv = 1.0 / t.tranMercK0R4
	return t, nil
}

func (t *TransverseMercator) generateCoefficients(invfla float64, n1 *float64,
	aCoeff []float64,
	bCoeff []float64,
	R4oa *float64,
	ellipsoidCode string) {
	/*  Generate Coefficients for Transverse Mercator algorithms
	===----- ===---------
	Algorithm developed by: C. Rollins   April 18, 2006

	INPUT
	-----
	invfla    Inverse flattening (reciprocal flattening)

	OUTPUT
	------
	n1        Helmert's "n"
	aCoeff    Coefficients for omega as a trig series in chi
	bBoeff    Coefficients for chi as a trig series in omega
	R4oa      Ratio "R4 over a", i.e. R4/a

	EXPLANATIONS
	------------
	omega is rectifying latitude
	chi is conformal latitude
	psi is geocentric latitude
	phi is geodetic latitude, commonly, "the latitude"
	R4 is the meridional isoperimetric radius
	"a" is the semi-major axis of the ellipsoid
	"b" is the semi-minor axis of the ellipsoid
	Helmert's n = (a - b)/(a + b)

	This calculation depends only on the shape of the ellipsoid and are
	independent of the ellipsoid size.

	The array Acoeff(8) stores eight coefficients corresponding
	to k = 2, 4, 6, 8, 10, 12, 14, 16 in the notation "a sub k".
	Likewise Bcoeff(8) etc.
	*/

	*n1 = 1.0 / (2*invfla - 1.0)

	n2 := *n1 * *n1
	n3 := n2 * *n1
	n4 := n3 * *n1
	n5 := n4 * *n1
	n6 := n5 * *n1
	n7 := n6 * *n1
	n8 := n7 * *n1
	n9 := n8 * *n1
	n10 := n9 * *n1

	// checks ellipsoid code and assigns values for corresponding coefficients.
	// Uses default computation if ellipsoid code isn't found. This will be
	// for user defined ellipsoids.

	switch ellipsoidCode {
	case "AA", "AM":
		aCoeff[0] = 8.3474517669594013740e-04
		aCoeff[1] = 7.554352936725572895e-07
		aCoeff[2] = 1.18487391005135489e-09
		aCoeff[3] = 2.3946872955703565e-12
		aCoeff[4] = 5.610633978440270e-15
		aCoeff[5] = 1.44858956458553e-17

		bCoeff[0] = -8.3474551646761162264e-04
		bCoeff[1] = -5.863630361809676570e-08
		bCoeff[2] = -1.65562038746920803e-10
		bCoeff[3] = -2.1340335537652749e-13
		bCoeff[4] = -3.720760760132477e-16
		bCoeff[5] = -7.08304328877781e-19
	case "EA", "EB", "EC", "ED", "EE":
		aCoeff[0] = 8.3064943111192510534e-04
		aCoeff[1] = 7.480375027595025021e-07
		aCoeff[2] = 1.16750772278215999e-09
		aCoeff[3] = 2.3479972304395461e-12
		aCoeff[4] = 5.474212231879573e-15
		aCoeff[5] = 1.40642257446745e-17

		bCoeff[0] = -8.3064976590443772201e-04
		bCoeff[1] = -5.805953517555717859e-08
		bCoeff[2] = -1.63133251663416522e-10
		bCoeff[3] = -2.0923797199593389e-13
		bCoeff[4] = -3.630200927775259e-16
		bCoeff[5] = -6.87666654919219e-19

	case "BN", "BR":
		aCoeff[0] = 8.3522527226849818552e-04
		aCoeff[1] = 7.563048340614894422e-07
		aCoeff[2] = 1.18692075307408346e-09
		aCoeff[3] = 2.4002054791393298e-12
		aCoeff[4] = 5.626801597980756e-15
		aCoeff[5] = 1.45360057224474e-17

		bCoeff[0] = -8.3522561262703079182e-04
		bCoeff[1] = -5.870409978661008580e-08
		bCoeff[2] = -1.65848307463131468e-10
		bCoeff[3] = -2.1389565927064571e-13
		bCoeff[4] = -3.731493368666479e-16
		bCoeff[5] = -7.10756898071999e-19
	case "KA", "HE", "FA":
		aCoeff[0] = 8.3761175713442343106e-04
		aCoeff[1] = 7.606346200814720197e-07
		aCoeff[2] = 1.19713032035541037e-09
		aCoeff[3] = 2.4277772986483520e-12
		aCoeff[4] = 5.707722772225013e-15
		aCoeff[5] = 1.47872454335773e-17

		bCoeff[0] = -8.3761210042019176501e-04
		bCoeff[1] = -5.904169154078546237e-08
		bCoeff[2] = -1.67276212891429215e-10
		bCoeff[3] = -2.1635549847939549e-13
		bCoeff[4] = -3.785212121016612e-16
		bCoeff[5] = -7.23053625983667e-19
	case "WD":
		aCoeff[0] = 8.3772481044362217923e-04
		aCoeff[1] = 7.608400388863560936e-07
		aCoeff[2] = 1.19761541904924067e-09
		aCoeff[3] = 2.4290893081322466e-12
		aCoeff[4] = 5.711579173743133e-15
		aCoeff[5] = 1.47992364667635e-17

		bCoeff[0] = -8.3772515386847544554e-04
		bCoeff[1] = -5.905770828762463028e-08
		bCoeff[2] = -1.67344058948464124e-10
		bCoeff[3] = -2.1647255130188214e-13
		bCoeff[4] = -3.787772179729998e-16
		bCoeff[5] = -7.23640523525528e-19
	case "WE":
		aCoeff[0] = 8.3773182062446983032e-04
		aCoeff[1] = 7.608527773572489156e-07
		aCoeff[2] = 1.19764550324249210e-09
		aCoeff[3] = 2.4291706803973131e-12
		aCoeff[4] = 5.711818369154105e-15
		aCoeff[5] = 1.47999802705262e-17

		bCoeff[0] = -8.3773216405794867707e-04
		bCoeff[1] = -5.905870152220365181e-08
		bCoeff[2] = -1.67348266534382493e-10
		bCoeff[3] = -2.1647981104903862e-13
		bCoeff[4] = -3.787930968839601e-16
		bCoeff[5] = -7.23676928796690e-19
	case "RF":
		aCoeff[0] = 8.3773182472855134012e-04
		aCoeff[1] = 7.608527848149655006e-07
		aCoeff[2] = 1.19764552085530681e-09
		aCoeff[3] = 2.4291707280369697e-12
		aCoeff[4] = 5.711818509192422e-15
		aCoeff[5] = 1.47999807059922e-17

		bCoeff[0] = -8.3773216816203523672e-04
		bCoeff[1] = -5.905870210369121594e-08
		bCoeff[2] = -1.67348268997717031e-10
		bCoeff[3] = -2.1647981529928124e-13
		bCoeff[4] = -3.787931061803592e-16
		bCoeff[5] = -7.23676950110361e-19
	case "SA", "AN":
		aCoeff[0] = 8.3775209887947194075e-04
		aCoeff[1] = 7.608896263599627157e-07
		aCoeff[2] = 1.19773253021831769e-09
		aCoeff[3] = 2.4294060763606098e-12
		aCoeff[4] = 5.712510331613028e-15
		aCoeff[5] = 1.48021320370432e-17

		bCoeff[0] = -8.3775244233790270051e-04
		bCoeff[1] = -5.906157468586898015e-08
		bCoeff[2] = -1.67360438158764851e-10
		bCoeff[3] = -2.1650081225048788e-13
		bCoeff[4] = -3.788390325953455e-16
		bCoeff[5] = -7.23782246429908e-19
	case "ID":
		aCoeff[0] = 8.3776052087969078729e-04
		aCoeff[1] = 7.609049308144604484e-07
		aCoeff[2] = 1.19776867565343872e-09
		aCoeff[3] = 2.4295038464530901e-12
		aCoeff[4] = 5.712797738386076e-15
		aCoeff[5] = 1.48030257891140e-17

		bCoeff[0] = -8.3776086434848497443e-04
		bCoeff[1] = -5.906276799395007586e-08
		bCoeff[2] = -1.67365493472742884e-10
		bCoeff[3] = -2.1650953495573773e-13
		bCoeff[4] = -3.788581120060625e-16
		bCoeff[5] = -7.23825990889693e-19
	case "IN", "HO":
		aCoeff[0] = 8.4127599100356448089e-04
		aCoeff[1] = 7.673066923431950296e-07
		aCoeff[2] = 1.21291995794281190e-09
		aCoeff[3] = 2.4705731165688123e-12
		aCoeff[4] = 5.833780550286833e-15
		aCoeff[5] = 1.51800420867708e-17

		bCoeff[0] = -8.4127633881644851945e-04
		bCoeff[1] = -5.956193574768780571e-08
		bCoeff[2] = -1.69484573979154433e-10
		bCoeff[3] = -2.2017363465021880e-13
		bCoeff[4] = -3.868896221495780e-16
		bCoeff[5] = -7.42279219864412e-19
	case "WO":
		aCoeff[0] = 8.4411652150600103279e-04
		aCoeff[1] = 7.724989750172583427e-07
		aCoeff[2] = 1.22525529789972041e-09
		aCoeff[3] = 2.5041361775549209e-12
		aCoeff[4] = 5.933026083631383e-15
		aCoeff[5] = 1.54904908794521e-17

		bCoeff[0] = -8.4411687285559594196e-04
		bCoeff[1] = -5.996681687064322548e-08
		bCoeff[2] = -1.71209836918814857e-10
		bCoeff[3] = -2.2316811233502163e-13
		bCoeff[4] = -3.934782433323038e-16
		bCoeff[5] = -7.57474665717687e-19

	case "CC":
		aCoeff[0] = 8.4703742793654652315e-04
		aCoeff[1] = 7.778564517658115212e-07
		aCoeff[2] = 1.23802665917879731e-09
		aCoeff[3] = 2.5390045684252928e-12
		aCoeff[4] = 6.036484469753319e-15
		aCoeff[5] = 1.58152259295850e-17

		bCoeff[0] = -8.4703778294785813001e-04
		bCoeff[1] = -6.038459874600183555e-08
		bCoeff[2] = -1.72996106059227725e-10
		bCoeff[3] = -2.2627911073545072e-13
		bCoeff[4] = -4.003466873888566e-16
		bCoeff[5] = -7.73369749524777e-19
	case "CG":
		aCoeff[0] = 8.5140099460764136776e-04
		aCoeff[1] = 7.858945456038187774e-07
		aCoeff[2] = 1.25727085106103462e-09
		aCoeff[3] = 2.5917718627340128e-12
		aCoeff[4] = 6.193726879043722e-15
		aCoeff[5] = 1.63109098395549e-17

		bCoeff[0] = -8.5140135513650084564e-04
		bCoeff[1] = -6.101145475063033499e-08
		bCoeff[2] = -1.75687742410879760e-10
		bCoeff[3] = -2.3098718484594067e-13
		bCoeff[4] = -4.107860472919190e-16
		bCoeff[5] = -7.97633133452512e-19

	case "CD":
		aCoeff[0] = 8.5140395445291970541e-04
		aCoeff[1] = 7.859000119464140978e-07
		aCoeff[2] = 1.25728397182445579e-09
		aCoeff[3] = 2.5918079321459932e-12
		aCoeff[4] = 6.193834639108787e-15
		aCoeff[5] = 1.63112504092335e-17

		bCoeff[0] = -8.5140431498554106268e-04
		bCoeff[1] = -6.101188106187092184e-08
		bCoeff[2] = -1.75689577596504470e-10
		bCoeff[3] = -2.3099040312610703e-13
		bCoeff[4] = -4.107932016207395e-16
		bCoeff[5] = -7.97649804397335e-19

	default:
		// computation below is for user defined ellipsoids
		// Computation of coefficient a2
		coeff := 0.0
		coeff += (-18975107.0) * n8 / 50803200.0
		coeff += (72161.0) * n7 / 387072.0
		coeff += (7891.0) * n6 / 37800.0
		coeff += (-127.0) * n5 / 288.0
		coeff += (41.0) * n4 / 180.0
		coeff += (5.0) * n3 / 16.0
		coeff += (-2.0) * n2 / 3.0
		coeff += (1.0) * *n1 / 2.0

		aCoeff[0] = coeff

		//   Computation of coefficient a4
		coeff = 0.0
		coeff += (148003883.0) * n8 / 174182400.0
		coeff += (13769.0) * n7 / 28800.0
		coeff += (-1983433.0) * n6 / 1935360.0
		coeff += (281.0) * n5 / 630.0
		coeff += (557.0) * n4 / 1440.0
		coeff += (-3.0) * n3 / 5.0
		coeff += (13.0) * n2 / 48.0

		aCoeff[1] = coeff

		//   Computation of coefficient a6
		coeff = 0.0
		coeff += (79682431.0) * n8 / 79833600.0
		coeff += (-67102379.0) * n7 / 29030400.0
		coeff += (167603.0) * n6 / 181440.0
		coeff += (15061.0) * n5 / 26880.0
		coeff += (-103.0) * n4 / 140.0
		coeff += (61.0) * n3 / 240.0

		aCoeff[2] = coeff

		//   Computation of coefficient a8
		coeff = 0.0
		coeff += (-40176129013.0) * n8 / 7664025600.0
		coeff += (97445.0) * n7 / 49896.0
		coeff += (6601661.0) * n6 / 7257600.0
		coeff += (-179.0) * n5 / 168.0
		coeff += (49561.0) * n4 / 161280.0

		aCoeff[3] = coeff

		//   Computation of coefficient a10
		coeff = 0.0
		coeff += (2605413599.0) * n8 / 622702080.0
		coeff += (14644087.0) * n7 / 9123840.0
		coeff += (-3418889.0) * n6 / 1995840.0
		coeff += (34729.0) * n5 / 80640.0

		aCoeff[4] = coeff

		//   Computation of coefficient a12
		coeff = 0.0
		coeff += (175214326799.0) * n8 / 58118860800.0
		coeff += (-30705481.0) * n7 / 10378368.0
		coeff += (212378941.0) * n6 / 319334400.0

		aCoeff[5] = coeff

		//   Computation of coefficient a14
		coeff = 0.0
		coeff += (-16759934899.0) * n8 / 3113510400.0
		coeff += (1522256789.0) * n7 / 1383782400.0

		aCoeff[6] = coeff

		//   Computation of coefficient a16
		coeff = 0.0
		coeff += (1424729850961.0) * n8 / 743921418240.0

		aCoeff[7] = coeff

		//   Computation of coefficient b2
		coeff = 0.0
		coeff += (-7944359.0) * n8 / 67737600.0
		coeff += (5406467.0) * n7 / 38707200.0
		coeff += (-96199.0) * n6 / 604800.0
		coeff += (81.0) * n5 / 512.0
		coeff += (1.0) * n4 / 360.0
		coeff += (-37.0) * n3 / 96.0
		coeff += (2.0) * n2 / 3.0
		coeff += (-1.0) * *n1 / 2.0

		bCoeff[0] = coeff

		//   Computation of coefficient b4
		coeff = 0.0
		coeff += (-24749483.0) * n8 / 348364800.0
		coeff += (-51841.0) * n7 / 1209600.0
		coeff += (1118711.0) * n6 / 3870720.0
		coeff += (-46.0) * n5 / 105.0
		coeff += (437.0) * n4 / 1440.0
		coeff += (-1.0) * n3 / 15.0
		coeff += (-1.0) * n2 / 48.0

		bCoeff[1] = coeff

		//   Computation of coefficient b6
		coeff = 0.0
		coeff += (6457463.0) * n8 / 17740800.0
		coeff += (-9261899.0) * n7 / 58060800.0
		coeff += (-5569.0) * n6 / 90720.0
		coeff += (209.0) * n5 / 4480.0
		coeff += (37.0) * n4 / 840.0
		coeff += (-17.0) * n3 / 480.0

		bCoeff[2] = coeff

		//   Computation of coefficient b8
		coeff = 0.0
		coeff += (-324154477.0) * n8 / 7664025600.0
		coeff += (-466511.0) * n7 / 2494800.0
		coeff += (830251.0) * n6 / 7257600.0
		coeff += (11.0) * n5 / 504.0
		coeff += (-4397.0) * n4 / 161280.0

		bCoeff[3] = coeff

		//   Computation of coefficient b10
		coeff = 0.0
		coeff += (-22894433.0) * n8 / 124540416.0
		coeff += (8005831.0) * n7 / 63866880.0
		coeff += (108847.0) * n6 / 3991680.0
		coeff += (-4583.0) * n5 / 161280.0

		bCoeff[4] = coeff

		//   Computation of coefficient b12
		coeff = 0.0
		coeff += (2204645983.0) * n8 / 12915302400.0
		coeff += (16363163.0) * n7 / 518918400.0
		coeff += (-20648693.0) * n6 / 638668800.0

		bCoeff[5] = coeff

		//   Computation of coefficient b14
		coeff = 0.0
		coeff += (497323811.0) * n8 / 12454041600.0
		coeff += (-219941297.0) * n7 / 5535129600.0

		bCoeff[6] = coeff

		//   Computation of coefficient b16
		coeff = 0.0
		coeff += (-191773887257.0) * n8 / 3719607091200.0

		bCoeff[7] = coeff
	}

	coeff := 0.0
	coeff += 49 * n10 / 65536.0
	coeff += 25 * n8 / 16384.0
	coeff += n6 / 256.0
	coeff += n4 / 64.0
	coeff += n2 / 4
	coeff++
	*R4oa = coeff / (1 + *n1)
}

func (t *TransverseMercator) checkLatLon(latitude, deltaLon float64) error {
	// test is based on distance from central meridian = deltaLon
	if deltaLon > math.Pi {
		deltaLon -= (2 * math.Pi)
	}
	if deltaLon < -math.Pi {
		deltaLon += (2 * math.Pi)
	}

	testAngle := math.Abs(deltaLon)

	delta := math.Abs(deltaLon - math.Pi)
	if delta < testAngle {
		testAngle = delta
	}

	delta = math.Abs(deltaLon + math.Pi)
	if delta < testAngle {
		testAngle = delta
	}

	// Away from the equator, is also valid
	delta = math.Pi/2 - latitude
	if delta < testAngle {
		testAngle = delta
	}

	delta = math.Pi/2 + latitude
	if delta < testAngle {
		testAngle = delta
	}
	const maxDeltaLong = ((math.Pi * 70) / 180.0)
	if testAngle > maxDeltaLong {
		return errors.New("longitude out of range")
	}
	return nil
}

func (t *TransverseMercator) latLonToNorthingEasting(latitude, longitude float64, northing, easting *float64) error {
	//  Convert longitude (Greenwhich) to longitude from the central meridian
	//  (-Pi, Pi] equivalent needed for checkLatLon.
	//  Compute its cosine and sine.
	lambda := longitude - t.tranMercOriginLong
	if lambda > math.Pi {
		lambda -= (2 * math.Pi)
	}
	if lambda < -math.Pi {
		lambda += (2 * math.Pi)
	}
	if err := t.checkLatLon(latitude, lambda); err != nil {
		return err
	}

	cosLam := math.Cos(lambda)
	sinLam := math.Sin(lambda)
	cosPhi := math.Cos(latitude)
	sinPhi := math.Sin(latitude)

	var c2ku, s2ku [8]float64
	var c2kv, s2kv [8]float64

	//  Ellipsoid to sphere
	//  --------- -- ------

	//  Convert geodetic latitude, Phi, to conformal latitude, Chi
	//  Only the cosine and sine of Chi are actually needed.
	P := math.Exp(t.tranMercEps * aTanH(t.tranMercEps*sinPhi))
	part1 := (1 + sinPhi) / P
	part2 := (1 - sinPhi) * P
	denom := part1 + part2
	cosChi := 2 * cosPhi / denom
	sinChi := (part1 - part2) / denom

	//  Sphere to first plane
	//  ------ -- ----- -----

	// Apply spherical theory of transverse Mercator to get (u,v) coord.s
	U := aTanH(cosChi * sinLam)
	V := math.Atan2(sinChi, cosChi*cosLam)

	// Use trig identities to compute cosh(2kU), sinh(2kU), cos(2kV), sin(2kV)
	computeHyperbolicSeries(2.0*U, c2ku[:], s2ku[:])
	computeTrigSeries(2.0*V, c2kv[:], s2kv[:])

	//  First plane to second plane
	//  Accumulate terms for X and Y
	xStar := 0.0
	yStar := 0.0

	for k := nTerms - 1; k >= 0; k-- {
		xStar += t.tranMercACoeff[k] * s2ku[k] * c2kv[k]
		yStar += t.tranMercACoeff[k] * c2ku[k] * s2kv[k]
	}

	xStar += U
	yStar += V

	// Apply isoperimetric radius, scale adjustment, and offsets
	*easting = (t.tranMercK0R4 * xStar)
	*northing = (t.tranMercK0R4 * yStar)
	return nil
}

func (t *TransverseMercator) convertFromGeodetic(geodeticCoordinates s2.LatLng) (MapCoords, error) {
	longitude := geodeticCoordinates.Lng.Radians()
	latitude := geodeticCoordinates.Lat.Radians()

	if longitude > math.Pi {
		longitude -= (2 * math.Pi)
	}
	if longitude < -math.Pi {
		longitude += (2 * math.Pi)
	}

	//  Convert longitude (Greenwhich) to longitude from the central meridian
	//  (-Pi, Pi] equivalent needed for checkLatLon.
	//  Compute its cosine and sine.
	lambda := longitude - t.tranMercOriginLong
	if lambda > math.Pi {
		lambda -= (2 * math.Pi)
	}
	if lambda < -math.Pi {
		lambda += (2 * math.Pi)
	}
	if err := t.checkLatLon(latitude, lambda); err != nil {
		return MapCoords{}, err
	}

	var easting, northing float64
	if err := t.latLonToNorthingEasting(latitude, longitude, &northing, &easting); err != nil {
		return MapCoords{}, err
	}

	// The origin may move form (0,0) and this is represented by
	// a change in the false Northing/Easting values.
	var falseEasting, falseNorthing float64
	if err := t.latLonToNorthingEasting(t.tranMercOriginLat, t.tranMercOriginLong,
		&falseNorthing, &falseEasting); err != nil {
		return MapCoords{}, err
	}

	easting += t.tranMercFalseEasting - falseEasting
	northing += t.TranMercFalseNorthing - falseNorthing

	invFlattening := 1.0 / t.flattening
	if invFlattening < 290.0 || invFlattening > 301.0 {
		//warning =  "Eccentricity is outside range that algorithm accuracy has been tested."
		// TODO: do we care about this?
	}

	return MapCoords{
		Easting:  easting,
		Northing: northing,
	}, nil
}

func (t *TransverseMercator) convertToGeodetic(mapProjectionCoordinates MapCoords) (s2.LatLng, error) {
	easting := mapProjectionCoordinates.Easting
	northing := mapProjectionCoordinates.Northing

	if (easting < (t.tranMercFalseEasting - t.tranMercDeltaEasting)) ||
		(easting > (t.tranMercFalseEasting +
			t.tranMercDeltaEasting)) {
		return s2.LatLng{}, errors.New("easting out of range")
	}

	if (northing < (t.TranMercFalseNorthing - t.tranMercDeltaNorthing)) ||
		(northing > (t.TranMercFalseNorthing +
			t.tranMercDeltaNorthing)) {
		return s2.LatLng{}, errors.New("northing out of range")
	}

	var longitude, latitude float64
	// The origin may move form (0,0) and this is represented by
	// a change in the false Northing/Easting values.
	var falseEasting, falseNorthing float64
	if err := t.latLonToNorthingEasting(t.tranMercOriginLat, t.tranMercOriginLong,
		&falseNorthing, &falseEasting); err != nil {
		return s2.LatLng{}, err
	}

	easting -= (t.tranMercFalseEasting - falseEasting)
	northing -= (t.TranMercFalseNorthing - falseNorthing)

	t.northingEastingToLatLon(northing, easting, &latitude, &longitude)

	if longitude > math.Pi {
		longitude = longitude - (2 * math.Pi)
	}
	if longitude <= -math.Pi {
		longitude = longitude + (2 * math.Pi)
	}

	if math.Abs(latitude) > (90.0 * math.Pi / 180.0) {
		return s2.LatLng{}, errors.New("northing out of range")
	}
	if (longitude) > (math.Pi) {
		longitude -= (2 * math.Pi)
		if math.Abs(longitude) > math.Pi {
			return s2.LatLng{}, errors.New("easting out of range")
		}
	} else if (longitude) < (-math.Pi) {
		longitude += (2 * math.Pi)
		if math.Abs(longitude) > math.Pi {
			return s2.LatLng{}, errors.New("easting out of range")
		}
	}

	invFlattening := 1.0 / t.flattening
	if invFlattening < 290.0 || invFlattening > 301.0 {
		// warning = "Eccentricity is outside range that algorithm accuracy has been tested."
		// TODO: do we care about the warning?
	}
	return s2.LatLng{Lat: s1.Angle(latitude), Lng: s1.Angle(longitude)}, nil
}

func (t *TransverseMercator) northingEastingToLatLon(northing,
	easting float64,
	latitude, longitude *float64) {
	var c2kx, s2kx, c2ky, s2ky [8]float64

	//  Undo offsets, scale change, and factor R4
	//  ---- -------  ----- ------  --- ------ --
	xStar := t.tranMercK0R4inv * (easting)
	yStar := t.tranMercK0R4inv * (northing)

	// Use trig identities to compute cosh(2kU), sinh(2kU), cos(2kV), sin(2kV)
	computeHyperbolicSeries(2.0*xStar, c2kx[:], s2kx[:])
	computeTrigSeries(2.0*yStar, c2ky[:], s2ky[:])

	//  Second plane (x*, y*) to first plane (u, v)
	//  ------ ----- -------- -- ----- ----- ------
	U := 0.0
	V := 0.0

	for k := nTerms - 1; k >= 0; k-- {
		U += t.tranMercBCoeff[k] * s2kx[k] * c2ky[k]
		V += t.tranMercBCoeff[k] * c2kx[k] * s2ky[k]
	}

	U += xStar
	V += yStar

	//  First plane to sphere
	//  ----- ----- -- ------
	coshU := math.Cosh(U)
	sinhU := math.Sinh(U)
	cosV := math.Cos(V)
	sinV := math.Sin(V)

	var lambda float64
	//   Longitude from central meridian
	if (math.Abs(cosV) < 10E-12) && (math.Abs(coshU) < 10E-12) {
		lambda = 0
	} else {
		lambda = math.Atan2(sinhU, cosV)
	}

	//   Conformal latitude
	sinChi := sinV / coshU
	*latitude = geodeticLat(sinChi, t.tranMercEps)

	// Longitude from Greenwich
	// --------  ---- ---------
	*longitude = t.tranMercOriginLong + lambda
}

func geodeticLat(sinChi, e float64) float64 {
	sOld := 1.0e99
	s := sinChi
	onePlusSinChi := 1.0 + sinChi
	oneMinusSinChi := 1.0 - sinChi

	for n := 0; n < 30; n++ {
		p := math.Exp(e * aTanH(e*s))
		pSq := p * p
		s = (onePlusSinChi*pSq - oneMinusSinChi) /
			(onePlusSinChi*pSq + oneMinusSinChi)

		if math.Abs(s-sOld) < 1.0e-12 {
			break
		}
		sOld = s
	}
	return math.Asin(s)
}

func computeHyperbolicSeries(twoX float64, c2kx, s2kx []float64) {
	// Use trig identities to compute
	// c2kx[k] = cosh(2kX), s2kx[k] = sinh(2kX)   for k = 0 .. 8
	c2kx[0] = math.Cosh(twoX)
	s2kx[0] = math.Sinh(twoX)
	c2kx[1] = 2.0*c2kx[0]*c2kx[0] - 1.0
	s2kx[1] = 2.0 * c2kx[0] * s2kx[0]
	c2kx[2] = c2kx[0]*c2kx[1] + s2kx[0]*s2kx[1]
	s2kx[2] = c2kx[1]*s2kx[0] + c2kx[0]*s2kx[1]
	c2kx[3] = 2.0*c2kx[1]*c2kx[1] - 1.0
	s2kx[3] = 2.0 * c2kx[1] * s2kx[1]
	c2kx[4] = c2kx[0]*c2kx[3] + s2kx[0]*s2kx[3]
	s2kx[4] = c2kx[3]*s2kx[0] + c2kx[0]*s2kx[3]
	c2kx[5] = 2.0*c2kx[2]*c2kx[2] - 1.0
	s2kx[5] = 2.0 * c2kx[2] * s2kx[2]
	c2kx[6] = c2kx[0]*c2kx[5] + s2kx[0]*s2kx[5]
	s2kx[6] = c2kx[5]*s2kx[0] + c2kx[0]*s2kx[5]
	c2kx[7] = 2.0*c2kx[3]*c2kx[3] - 1.0
	s2kx[7] = 2.0 * c2kx[3] * s2kx[3]
}

func computeTrigSeries(twoY float64, c2ky,
	s2ky []float64) {
	// Use trig identities to compute
	// c2ky[k] = cos(2kY), s2ky[k] = sin(2kY)   for k = 0 .. 8
	c2ky[0] = math.Cos(twoY)
	s2ky[0] = math.Sin(twoY)
	c2ky[1] = 2.0*c2ky[0]*c2ky[0] - 1.0
	s2ky[1] = 2.0 * c2ky[0] * s2ky[0]
	c2ky[2] = c2ky[1]*c2ky[0] - s2ky[1]*s2ky[0]
	s2ky[2] = c2ky[1]*s2ky[0] + c2ky[0]*s2ky[1]
	c2ky[3] = 2.0*c2ky[1]*c2ky[1] - 1.0
	s2ky[3] = 2.0 * c2ky[1] * s2ky[1]
	c2ky[4] = c2ky[3]*c2ky[0] - s2ky[3]*s2ky[0]
	s2ky[4] = c2ky[3]*s2ky[0] + c2ky[0]*s2ky[3]
	c2ky[5] = 2.0*c2ky[2]*c2ky[2] - 1.0
	s2ky[5] = 2.0 * c2ky[2] * s2ky[2]
	c2ky[6] = c2ky[5]*c2ky[0] - s2ky[5]*s2ky[0]
	s2ky[6] = c2ky[5]*s2ky[0] + c2ky[0]*s2ky[5]
	c2ky[7] = 2.0*c2ky[3]*c2ky[3] - 1.0
	s2ky[7] = 2.0 * c2ky[3] * s2ky[3]
}
func aTanH(x float64) float64 {
	return (0.5 * math.Log((1+x)/(1-x)))
}
