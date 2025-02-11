package glatl

import (
	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// GSO returns squared norm of GSO-matrix and GSO-coefficient matrix
func GSO(b Lattice) (vec.Vector, mat.Matrix) {
	gsoB := vec.ZeroVec(b.NumRows)
	gsoMat := mat.ZeroMat(b.NumRows, b.NumCols)
	mu := mat.Identity(b.NumRows)

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

		normTemp = 0
		for j := 0; j < b.NumCols; j++ {
			normTemp += gsoMat.At[i][j] * gsoMat.At[i][j]
		}
		gsoB.At[i] = normTemp
	}

	return gsoB, mu
}
