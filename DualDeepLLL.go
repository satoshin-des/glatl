package glatl

import (
	"math"

	"github.com/satoshin-des/glal/mat"
)

// dual version of DeepLLL
func DualDeepLLL(b Lattice, delta float64) {
	if delta < 0.25 || delta > 1 {
		panic("reduction parameter must be in [1/4, 1]")
	}

	var d float64
	var q float64
	var l int

	gsoB, mu := GSO(b)
	dualMu := mat.Identity(b.NumRows)

	for k := b.NumRows - 2; k >= 0; {
		dualMu.At[k][k] = 1.0

		for j := k + 1; j < b.NumRows; j++ {
			dualMu.At[k][j] = 0
			for i := k; i < j; i++ {
				dualMu.At[k][j] -= mu.At[j][i] * dualMu.At[k][i]
			}

			if dualMu.At[k][j] > 0.5 || dualMu.At[k][j] < -0.5 {
				q = math.Round(dualMu.At[k][j])
				for i := 0; i < b.NumCols; i++ {
					b.Basis[j][i] += int64(q) * b.Basis[k][i]
				}
				for i := j; i < b.NumRows; i++ {
					dualMu.At[k][i] -= q * dualMu.At[j][i]
				}
				for i := 0; i <= k; i++ {
					mu.At[j][i] += q * mu.At[k][i]
				}
			}
		}

		d = 0
		l = b.NumRows - 1
		for j := k; j < b.NumRows; j++ {
			d += dualMu.At[k][j] * dualMu.At[k][j] / gsoB.At[j]
		}
		for l > k {
			if gsoB.At[l]*d < delta {
				DualDeepIns(b, k, l)

				if l < b.NumRows-2 {
					k = l + 1
				} else {
					k = b.NumRows - 1
				}

				gsoB, mu = GSO(b)
			} else {
				d -= dualMu.At[k][l] * dualMu.At[k][l] / gsoB.At[l]
				l--
			}
		}
		k--
	}
}
