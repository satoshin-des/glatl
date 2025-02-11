package glatl

import (
	"math/rand"
)

func RandBasis(n int) Lattice {
	var mat Lattice

	mat.NumRows = n
	mat.NumCols = n
	mat.Basis = make([][]int64, n)
	for i := 0; i < mat.NumRows; i++ {
		mat.Basis[i] = make([]int64, n)
		for j := 0; j < mat.NumCols; j++ {
			mat.Basis[i][j] = 0
		}
		mat.Basis[i][i] = 1
		mat.Basis[i][0] = rand.Int63n(99999) + 10000
	}

	return mat
}
