package coordconv_test

import (
	"testing"

	"github.com/golang/geo/s2"
	"github.com/tzneal/coordconv"
)

func TestUPSRoundTrip(t *testing.T) {
	ups, err := coordconv.NewUPS(ellipsoidSemiMajorAxis, ellipsoidFlattening)
	if err != nil {
		t.Fatalf("error creating UPS converter: %s", err)
	}
	const latInc = 0.5
	const lngInc = 0.5
	for lng := -190.0; lng < 190; lng += lngInc {
		for lat := -100.0; lat < 100; lat += latInc {
			geo := s2.LatLngFromDegrees(lat, lng)
			uc, err := ups.ConvertFromGeodetic(geo)
			if err == nil {
				geo2, err := ups.ConvertToGeodetic(uc)
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
