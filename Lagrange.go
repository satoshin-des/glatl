package glatl

import "math"

// Lagrange computes Lagrange reduced basis
//
// panic if dimension of b is not 2
func Lagrange(b Lattice) {
	if b.NumRows != 2 {
		panic("2 dimensional lattice only can be Lagrange-reduced")
	}

	var normTemp1, normTemp2, v int64 = 0, 0, 0

	for j := 0; j < b.NumCols; j++ {
		normTemp1 += b.Basis[0][j] * b.Basis[0][j]
		normTemp2 += b.Basis[1][j] * b.Basis[1][j]
	}

	if normTemp1 > normTemp2 {
		for j := 0; j < b.NumCols; j++ {
			v = b.Basis[0][j]
			b.Basis[0][j] = b.Basis[1][j]
			b.Basis[1][j] = v
		}
	}

	for {
		normTemp1 = 0
		normTemp2 = 0
		for j := 0; j < b.NumCols; j++ {
			normTemp1 += b.Basis[0][j] * b.Basis[1][j]
			normTemp2 += b.Basis[0][j] * b.Basis[0][j]
		}

		for j := 0; j < b.NumCols; j++ {
			v = b.Basis[1][j] - int64(math.Round(float64(normTemp1)/float64(normTemp2)))*b.Basis[0][j]
			b.Basis[1][j] = b.Basis[0][j]
			b.Basis[0][j] = v
		}

		normTemp1 = 0
		normTemp2 = 0
		for j := 0; j < b.NumCols; j++ {
			normTemp1 += b.Basis[0][j] * b.Basis[0][j]
			normTemp2 += b.Basis[1][j] * b.Basis[1][j]
		}

		if normTemp1 >= normTemp2 {
			for j := 0; j < b.NumCols; j++ {
				v = b.Basis[0][j]
				b.Basis[0][j] = b.Basis[1][j]
				b.Basis[1][j] = v
			}
			break
		}
	}
}
