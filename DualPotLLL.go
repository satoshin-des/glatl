package glatl

import (
	"math"

	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// dual version of PotLLL
//
// panic if delta < 1/4 or delta > 1
//
// 佐藤 新, 安田 雅哉. 自己双対型PotBKZ基底簡約の提案とBKZとの比較. (2025)
func DualPotLLL(b Lattice, delta float64) {
	if delta < 0.25 || delta > 1 {
		panic("reduction parameter must be in [1/4, 1]")
	}

	var p float64
	var pMin float64
	var q float64
	var s float64
	var d float64
	var l int

	dualD := vec.ZeroVec(b.NumRows)
	dualMu := mat.Identity(b.NumRows)

	LLL(b, delta)
	gsoB, mu := GSO(b)

	for k := b.NumRows - 1; k >= 0; {
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

		p = 1.0
		pMin = 1.0
		l = b.NumRows - 1
		for j := k + 1; j < b.NumRows; j++ {
			s = 0.0
			for i := k; i <= j; i++ {
				s += dualMu.At[k][i] * dualMu.At[k][i] / gsoB.At[i]
			}
			p *= gsoB.At[j]
			p *= s

			if p < pMin {
				l = j
				pMin = p
			}
		}

		if delta > pMin {
			d = 1.0 / gsoB.At[k]
			vec.SetZero(dualD)
			dualD.At[k] = d
			for h := k + 1; h < b.NumRows; h++ {
				d += dualMu.At[k][h] * dualMu.At[k][h] / gsoB.At[h]
				dualD.At[h] = d
			}

			DualDeepIns(b, k, l)
			updateDualDeepInsGSO(k, l, gsoB, mu, dualMu, dualD)

			k = l
		} else {
			k--
		}
	}
}
