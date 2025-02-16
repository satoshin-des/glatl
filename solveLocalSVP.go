package glatl

import (
	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// SolveSVP returns shortest vector on lattice spaned by b
func solveLocalSVP(mu mat.Matrix, gsoB vec.Vector) vec.Vector {
	enumVec := vec.ZeroVec(mu.NumRows)
	oldEnumVec := vec.ZeroVec(mu.NumRows)

	for R := gsoB.At[0]; ; R *= 0.9801 {
		for i := 0; i < mu.NumRows; i++ {
			oldEnumVec.At[i] = enumVec.At[i]
		}

		enumVec = ENUM(mu, gsoB, R)

		if vec.IsZero(enumVec) {
			return oldEnumVec
		}
	}
}
