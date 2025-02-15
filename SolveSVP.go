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

	enumVec := vec.ZeroVec(b.NumRows)
	oldEnumVec := vec.ZeroVec(b.NumRows)
	rho := vec.ZeroVec(b.NumRows + 1)
	oldRho := vec.ZeroVec(b.NumRows + 1)

	for R := gsoB.At[0]; ; R *= 0.9801 {
		oldRho.At[0] = rho.At[0]
		for i := 0; i < b.NumRows; i++ {
			oldRho.At[i+1] = rho.At[i+1]
			oldEnumVec.At[i] = enumVec.At[i]
		}

		enumVec = ENUM(mu, gsoB, R, rho)

		if vec.IsZero(enumVec) {
			return glal.Mul(oldEnumVec, basisMat)
		}
	}
}
