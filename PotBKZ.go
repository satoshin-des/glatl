package glatl

import (
	"math"

	"github.com/satoshin-des/glal"
	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// PotBKZ computes reduced basis
//
// 佐藤 新, 安田 雅哉. 自己双対型PotBKZ基底簡約の提案とBKZとの比較. (2025)
func PotBKZ(b Lattice, beta int, delta float64) {
	var j int = 0
	var d int
	var k int
	var v vec.Vector
	var w vec.Vector
	var gsoB_ vec.Vector
	var logB_ vec.Vector
	var mu_ mat.Matrix
	var basisMat mat.Matrix

	PotLLL(b, delta)
	gsoB, mu := GSO(b)

	for z := 0; z < b.NumRows-1; {
		if j == b.NumRows-2 {
			j = 0
		}
		j++

		if j+beta-1 < b.NumRows-1 {
			k = j + beta - 1
		} else {
			k = b.NumRows - 1
		}

		d = k - j + 1

		gsoB_ = vec.ZeroVec(d)
		logB_ = vec.ZeroVec(d)
		mu_ = mat.Identity(d)
		for i := 0; i < d; i++ {
			gsoB_.At[i] = gsoB.At[i+j-1]
			logB_.At[i] = math.Log(gsoB_.At[i])
			for l := 0; l < d; l++ {
				mu_.At[i][l] = mu.At[i+j-1][l+j-1]
			}
		}

		v = PotENUM(mu_, gsoB_, logB_)

		if !vec.IsZero(v) {
			z = 0

			basisMat = mat.ZeroMat(d, b.NumCols)
			for i := 0; i < d; i++ {
				for l := 0; l < b.NumCols; l++ {
					basisMat.At[i][l] = float64(b.Basis[j-1+i][l])
				}
			}
			w = glal.Mul(v, basisMat)

			for i := d - 1; i >= 0; i-- {
				if v.At[i] == 1 || v.At[i] == -1 {
					for l := 0; l < b.NumCols; l++ {
						b.Basis[i+j-1][l] = int64(w.At[l])
					}
					DeepIns(b, j, i+j-1)
					break
				}
			}

			PotLLL(b, delta)
			gsoB, mu = GSO(b)
		} else {
			z++
		}
	}
}
