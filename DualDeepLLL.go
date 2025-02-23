package glatl

import (
	"math"

	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// dual version of DeepLLL
//
// panic if delta < 1/4 or delta > 1
//
// M. Yasuda, J. Yamaguchi, M. Ooka, S. Nakamura. Development of a dual version of DeepBKZ and its application to solving the LWE challenge.(2018)
func DualDeepLLL(b Lattice, delta float64) {
	if delta < 0.25 || delta > 1 {
		panic("reduction parameter must be in [1/4, 1]")
	}

	var d float64
	var q float64
	var l int
	var rho float64
	dualRho := vec.ZeroVec(b.NumRows)

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
				rho = 1.0 / gsoB.At[k]

				vec.SetZero(dualRho)
				dualRho.At[k] = rho
				for h := k + 1; h < b.NumRows; h++ {
					rho += dualMu.At[k][h] * dualMu.At[k][h] / gsoB.At[h]
					dualRho.At[h] = rho
				}

				DualDeepIns(b, k, l)
				updateDualDeepInsGSO(k, l, gsoB, mu, dualMu, dualRho)

				if l < b.NumRows-2 {
					k = l + 1
				} else {
					k = b.NumRows - 1
				}
			} else {
				d -= dualMu.At[k][l] * dualMu.At[k][l] / gsoB.At[l]
				l--
			}
		}
		k--
	}
}
