package glatl

import "math"

// Vol returns volume of b
func Vol(b Lattice) int64 {
	gsoB, _ := GSO(b)
	var p float64 = 1

	for i := 0; i < b.NumRows; i++ {
		p *= math.Sqrt(gsoB.At[i])
	}

	return int64(math.Round(p))
}
