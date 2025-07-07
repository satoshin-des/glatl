package glatl

import (
	"github.com/satoshin-des/glal"
	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// SolveSVP returns shortest vector on lattice spaned by b
func SolveSVP(b Lattice) vec.Vector {
	gsoB, mu := GSO(b)
	basisMat := mat.ZeroMat(b.NumRows, b.NumCols)

	for i := 0; i < b.NumRows; i++ {
		for j := 0; j < b.NumCols; j++ {
			basisMat.At[i][j] = float64(b.Basis[i][j])
		}
	}

	return glal.Mul(ENUM(mu, gsoB, gsoB.At[0]), basisMat)
}
