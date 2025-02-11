package glatl

import (
	"math"

	"github.com/satoshin-des/glal/mat"
)

// PartSize partialy size reduce b
func PartSize(b Lattice, mu mat.Matrix, i int, j int) {
	if mu.At[i][j] > 0.5 || mu.At[i][j] < -0.5 {
		q := int64(math.Round(mu.At[i][j]))
		for k := 0; k < b.NumCols; k++ {
			b.Basis[i][k] -= q * b.Basis[j][k]
		}
		for k := 0; k < j+1; k++ {
			mu.At[i][k] -= float64(q) * mu.At[j][k]
		}
	}
}
