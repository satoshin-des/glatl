package glatl

import (
	"math"

	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// BabaiNP returns a solution of approx-CVP on the lattice b target to w using Babai's nearest plane algorithm
//
// panic if length of w is not equal to length of b[i]
//
// L. Babai. On Lovasz' lattice reduction and the nearest lattice point problem. (1986)
func BabaiNP(b Lattice, w vec.Vector) vec.Vector {
	if b.NumCols != w.Length {
		panic("size of target vector and b are not equal")
	}

	gsoMat := mat.ZeroMat(b.NumRows, b.NumCols)
	mu := mat.Identity(b.NumRows)
	v := vec.ZeroVec(w.Length)

	var c float64
	var dotTemp, normTemp float64

	for i := 0; i < b.NumRows; i++ {
		mu.At[i][i] = 1

		for j := 0; j < b.NumCols; j++ {
			gsoMat.At[i][j] = float64(b.Basis[i][j])
		}

		for j := 0; j < i; j++ {
			dotTemp = 0
			normTemp = 0
			for k := 0; k < b.NumCols; k++ {
				dotTemp += float64(b.Basis[i][k]) * gsoMat.At[j][k]
				normTemp += gsoMat.At[j][k] * gsoMat.At[j][k]
			}

			mu.At[i][j] = dotTemp / normTemp

			for k := 0; k < b.NumCols; k++ {
				gsoMat.At[i][k] -= mu.At[i][j] * gsoMat.At[j][k]
			}
		}
	}

	for i := 0; i < v.Length; i++ {
		v.At[i] = w.At[i]
	}

	for i := b.NumRows - 1; i >= 0; i-- {
		dotTemp = 0
		for j := 0; j < b.NumCols; j++ {
			dotTemp += float64(v.At[j]) * gsoMat.At[i][j]
		}

		normTemp = 0
		for j := 0; j < b.NumCols; j++ {
			normTemp += gsoMat.At[i][j] * gsoMat.At[i][j]
		}

		c = math.Round(dotTemp / normTemp)
		for j := 0; j < b.NumCols; j++ {
			v.At[j] -= c * float64(b.Basis[i][j])
		}
	}

	return vec.Sub(w, v)
}
