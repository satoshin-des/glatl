package glatl

import (
	"math"

	"github.com/satoshin-des/glal"
	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// BabaiNP returns a solution of approx-CVP on the lattice b target to w using Babai's rounding
//
// panic if length of w is not equal to length of b[i]
func BabaiRound(b Lattice, w vec.Vector) vec.Vector {
	basisMat := mat.ZeroMat(b.NumRows, b.NumCols)
	for i := 0; i < b.NumRows; i++ {
		for j := 0; j < b.NumCols; j++ {
			basisMat.At[i][j] = float64(b.Basis[i][j])
		}
	}

	coeff := glal.Mul(w, mat.PseudoInv(basisMat))
	v := vec.ZeroVec(w.Length)

	for j := 0; j < b.NumCols; j++ {
		v.At[j] = 0
		for i := 0; i < b.NumCols; i++ {
			v.At[j] += math.Round(coeff.At[i]) * float64(b.Basis[i][j])
		}
	}

	return v
}
