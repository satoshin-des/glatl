package glatl

import (
	"github.com/satoshin-des/glal/vec"
)

// BKZ computes "BKZ-reduced" basis without MLLL to abord faster
//
// panic if delta < 1/4 or delta > 1 or if beta <= 2 or beta > dim(b)
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
	var v vec.Vector = vec.ZeroVec(b.NumRows)
	var w vec.Vector

	LLL(b, 0.99)

	gsoB, mu := GSO(b)

	var k int = 0
	var l int
	var h int
	var bkzTours int = 0

	for z := 0; z < b.NumRows-1; {
		if maxLoop > 0 && bkzTours >= maxLoop {
			break
		}

		if k == b.NumRows-1 {
			k = 0
			bkzTours++
		}
		k++

		if k+beta-1 < b.NumRows {
			l = k + beta - 1
		} else {
			l = b.NumRows
		}

		if l+1 < b.NumRows {
			h = l + 1
		} else {
			h = b.NumRows
		}

		d = l - k + 1

		w = LocalENUM(mu, gsoB, 0.99*gsoB.At[k-1], k-1, l)

		if !vec.IsZero(w) {
			z = 0

			for i := 0; i < b.NumCols; i++ {
				v.At[i] = 0
				for j := 0; j < d; j++ {
					v.At[i] += w.At[j] * float64(b.Basis[j+k-1][i])
				}
			}

			for i := d - 1; i >= 0; i-- {
				if w.At[i] == 1 || w.At[i] == -1 {
					for j := 0; j < b.NumCols; j++ {
						b.Basis[i+k-1][j] = int64(v.At[j])
					}
					DeepIns(b, k-1, i+k-1)
					break
				}
			}

			partLLL(b, delta, k-1, h)
		} else {
			z++

			partLLL(b, delta, h-2, h)
		}

		gsoB, mu = GSO(b)
	}
}
