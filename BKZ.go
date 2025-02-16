package glatl

import (
	"github.com/satoshin-des/glal"
	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// BKZ computes "BKZ-reduced" basis without MLLL to abord faster
//
// // panic if delta < 1/4 or delta > 1 or if beta <= 2 or beta > dim(b)
//
// C. P. Schnorr, M. Euchner. Lattice basis reduction: Improved practical algorithms and solving subset sum problem.(1994)
func BKZ(b Lattice, beta int, delta float64) {
	if delta < 0.25 || delta > 1 {
		panic("reduction parameter must be in [1/4, 1]")
	}

	if beta <= 1 || beta > b.NumRows {
		panic("block size is invalid. block size must be in [2, dim(b)].")
	}

	var d int
	var basisMat mat.Matrix
	var shortestVec vec.Vector
	var v vec.Vector
	var vTemp vec.Vector
	var mu_ mat.Matrix
	var gsoB_ vec.Vector

	gsoB, mu := GSO(b)

	var k int = 0
	var l int
	var bkzTours int = 0

	for z := 0; z < b.NumRows-1; {
		if maxLoop > 0 && bkzTours >= maxLoop {
			break
		}

		if k == b.NumRows-1 {
			k = 0
			z++
			bkzTours++
		}
		k++

		if k+beta-1 < b.NumRows {
			l = k + beta - 1
		} else {
			l = b.NumRows
		}

		d = l - k + 1

		gsoB_ = vec.ZeroVec(d)
		mu_ = mat.Identity(d)
		for i := 0; i < d; i++ {
			gsoB_.At[i] = gsoB.At[i+k-1]
			for j := 0; j < d; j++ {
				mu_.At[i][j] = mu.At[i+k-1][j+k-1]
			}
		}

		v = solveLocalSVP(mu_, gsoB_)
		vTemp = vec.ZeroVec(v.Length)

		for i := 0; i < v.Length; i++ {
			vTemp.At[i] = v.At[i]
		}
		vTemp.At[0] -= 1

		if !vec.IsZero(vTemp) {
			z = 0

			basisMat = mat.ZeroMat(d, b.NumCols)
			for i := 0; i < d; i++ {
				for j := 0; j < b.NumCols; j++ {
					basisMat.At[i][j] = float64(b.Basis[k-1+i][j])
				}
			}
			shortestVec = glal.Mul(v, basisMat)

			for i := d - 1; i >= 0; i-- {
				if v.At[i] == 1 || v.At[i] == -1 {
					for j := 0; j < b.NumCols; j++ {
						b.Basis[i+k-1][j] = int64(shortestVec.At[j])
					}
					break
				}
			}
		} else {
			z++
		}

		LLL(b, delta)
		gsoB, mu = GSO(b)
	}
}
