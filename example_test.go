package coordconv_test

import (
	"fmt"

	"github.com/golang/geo/s2"
	"github.com/tzneal/coordconv"
)

func ExampleToMGRS() {
	mgrs, _ := coordconv.DefaultMGRSConverter.ConvertFromGeodetic(s2.LatLngFromDegrees(0, 0), 5)
	fmt.Println(mgrs)
}
func ExampleFromMGRS() {
	geo, _ := coordconv.DefaultMGRSConverter.ConvertToGeodetic("16SGC3855124838")
	fmt.Println(geo)
}
