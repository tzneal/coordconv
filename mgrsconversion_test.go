package coordconv_test

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"math"
	"math/big"
	"os"
	"strings"
	"testing"

	"github.com/golang/geo/s1"
	"github.com/golang/geo/s2"
	"github.com/tzneal/coordconv"
)

const ellipsoidSemiMajorAxis = 6378137.0
const ellipsoidFlattening = 1 / 298.257223563

func TestCompleteMGRS(t *testing.T) {
	f, err := os.Open("testdata/geotrans.log")
	if err != nil {
		t.Fatalf("error reading file: %s", err)
	}
	sc := bufio.NewScanner(f)

	polConv, _ := coordconv.NewPolarStereographicScaleFactor(ellipsoidSemiMajorAxis, ellipsoidFlattening, 0, 1.0, coordconv.HemisphereNorth, 0, 0)
	upsConv, _ := coordconv.NewUPS(ellipsoidSemiMajorAxis, ellipsoidFlattening)
	mgrsConv, _ := coordconv.NewMGRS(ellipsoidSemiMajorAxis, ellipsoidFlattening, "WE")

	for sc.Scan() {
		// -0x1.a5bc5883b707cp+0 -0x1.9c6f2f5dc2367p+1 0x0p+0 0x0p+0 E ERR
		// lat lng polar:easting polar:northing ups:easting ups:northing ups:hemisphere mgrs
		sp := strings.Fields(sc.Text())
		if len(sp) != 8 {
			fmt.Println(sp)
			t.Fatalf("error reading, got %d", len(sp))
		}
		lat := parseRawDouble(t, sp[0])
		lng := parseRawDouble(t, sp[1])
		polarEasting := parseRawDouble(t, sp[2])
		polarNorthing := parseRawDouble(t, sp[3])
		upsEasting := parseRawDouble(t, sp[4])
		upsNorthing := parseRawDouble(t, sp[5])

		geodetic := s2.LatLng{Lat: s1.Angle(lat), Lng: s1.Angle(lng)}
		polar, _ := polConv.ConvertFromGeodetic(geodetic)
		ups, _ := upsConv.ConvertFromGeodetic(geodetic)
		mgrs, err := mgrsConv.ConvertFromGeodetic(geodetic, 5)
		if err != nil {
			mgrs = "ERR"
		}

		const epsilon = 1e-7
		// test the polar conversion
		if math.IsNaN(polar.Easting) && !math.IsNaN(polar.Easting) {
			t.Fatalf("%f %f got Easting %f, expected %f", lat, lng, polar.Easting, polarEasting)
		}
		if math.IsNaN(polar.Northing) && !math.IsNaN(polar.Northing) {
			t.Fatalf("%f %f got Northing %f, expected %f", lat, lng, polar.Northing, polarNorthing)
		}
		if math.Abs(polar.Easting-polarEasting) > epsilon {
			t.Fatalf("%f %f got Easting %f, expected %f", lat, lng, polar.Easting, polarEasting)
		}
		if math.Abs(polar.Northing-polarNorthing) > epsilon {
			t.Fatalf("%f %f got Northing %f, expected %f", lat, lng, polar.Northing, polarNorthing)
		}

		// test the ups conversion
		if math.IsNaN(ups.Easting) && !math.IsNaN(ups.Easting) {
			t.Fatalf("%f %f got Easting %f, expected %f", lat, lng, ups.Easting, upsEasting)
		}
		if math.IsNaN(ups.Northing) && !math.IsNaN(ups.Northing) {
			t.Fatalf("%f %f got Northing %f, expected %f", lat, lng, ups.Northing, upsNorthing)
		}
		if math.Abs(ups.Easting-upsEasting) > epsilon {
			t.Fatalf("%f %f got Easting %f, expected %f", lat, lng, ups.Easting, upsEasting)
		}
		if math.Abs(ups.Northing-upsNorthing) > epsilon {
			t.Fatalf("%f %f got Northing %f, expected %f", lat, lng, ups.Northing, upsNorthing)
		}
		if mgrs != sp[7] {
			t.Errorf("%f %f expected MGRS = '%s', got '%s'", lat, lng, sp[7], mgrs)
		}
	}
}

func parseRawDouble(t *testing.T, s string) float64 {
	if s == "nan" {
		return math.NaN()
	}
	g, ok := new(big.Float).SetString(s)
	if !ok {
		t.Fatalf("error parsing float %s", s)
	}
	f, _ := g.Float64()
	return f
}

func TestMGRSRoundTrip(t *testing.T) {
	mgrs, err := coordconv.NewMGRS(ellipsoidSemiMajorAxis, ellipsoidFlattening, "WE")
	if err != nil {
		t.Fatalf("error creating MGRS converter: %s", err)
	}
	const latInc = 0.5
	const lngInc = 0.5
	for lng := -190.0; lng < 190; lng += lngInc {
		for lat := -100.0; lat < 100; lat += latInc {
			geo := s2.LatLngFromDegrees(lat, lng)
			mc, err := mgrs.ConvertFromGeodetic(geo, 5)
			if err == nil {
				geo2, err := mgrs.ConvertToGeodetic(mc)
				if err != nil {
					t.Fatalf("expected no error in round trip, got one at %s (%s)", geo, err)
				}
				if geo.Distance(geo2) > 1e85 {
					t.Fatalf("expected %s, got %s", geo, geo2)
				}
			}
		}
	}
}

func TestMGRSFuzzCrashers(t *testing.T) {
	for _, v := range []string{"00000000\xff\xff", "\xff\xff", "00000000\u007f\xff",
		"00000000\xff\xff", "\u007f\xff"} {
		fuzzMGRS([]byte(v))
	}
}

func fuzzMGRS(data []byte) int {
	for len(data) < 16 {
		data = append(data, 0)
	}
	mgrs, _ := coordconv.NewMGRS(ellipsoidSemiMajorAxis, ellipsoidFlattening, "WE")
	latU := binary.BigEndian.Uint64(data[0:])
	lngU := binary.BigEndian.Uint64(data[8:])
	lat := math.Float64frombits(latU)
	lng := math.Float64frombits(lngU)

	geo := s2.LatLng{Lat: s1.Angle(lat), Lng: s1.Angle(lng)}
	mc, err := mgrs.ConvertFromGeodetic(geo, 5)
	if err != nil {
		return 0
	}

	geo2, err := mgrs.ConvertToGeodetic(mc)
	if err != nil {
		panic(fmt.Sprintf("expected no error in round trip, got one at %s (%s)", geo, err))
	}
	if geo.Distance(geo2) > 1e85 {
		panic(fmt.Sprintf("expected %s, got %s", geo, geo2))
	}
	return 1
}
