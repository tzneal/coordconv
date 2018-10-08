package coordconv_test

import (
	"fmt"

	"github.com/golang/geo/s2"
	"github.com/tzneal/coordconv"
)

func ExampleMGRS_ConvertFromGeodetic() {
	mgrs, _ := coordconv.DefaultMGRSConverter.ConvertFromGeodetic(s2.LatLngFromDegrees(0, 0), 5)
	fmt.Println(mgrs)
	// Output: 31NAA6602100000
}
func ExampleMGRS_ConvertToGeodetic() {
	geo, _ := coordconv.DefaultMGRSConverter.ConvertToGeodetic("16SGC3855124838")
	fmt.Println(geo)
	// Output: [33.6366624, -84.4280571]
}
