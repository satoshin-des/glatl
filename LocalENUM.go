package glatl

import (
	"math"

	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// ENUM returns a lattice vector whose norm is lesser than or equal to R on lattice spaned by b
//
// N. Gama, P. Q. Nguyen, O. Regev. Lattice enumeration using extreme pruning.(2010)
func LocalENUM(mu mat.Matrix, gsoB vec.Vector, R float64, start int, end int) vec.Vector {
	var n int = end - start
	var hasSolution bool = false
	var temp float64
	var lastNonzero int = 0

	r := vec.ZeroVec(n + 1)
	w := vec.ZeroVec(n)
	v := vec.ZeroVec(n)
	coeffVector := vec.ZeroVec(n)
	c := vec.ZeroVec(n)
	rho := vec.ZeroVec(n + 1)
	sigma := mat.ZeroMat(n+1, n)

	v.At[0] = 1
	for i := 0; i < n; i++ {
		r.At[i] = float64(i)
	}

	for k := 0; ; {
		//println("k=", k)
		temp = v.At[k] - c.At[k]
		temp *= temp
		rho.At[k] = rho.At[k+1] + temp*gsoB.At[k+start]
		if rho.At[k] <= R {
			if k == 0 {
				hasSolution = true
				for i := 0; i < coeffVector.Length; i++ {
					coeffVector.At[i] = v.At[i]
				}
				R = math.Min(0.99*rho.At[0], R)
			} else {
				k--
				if r.At[k] <= r.At[k+1] {
					r.At[k] = r.At[k+1]
				}
				for i := int(r.At[k]); i > k; i-- {
					sigma.At[i][k] = sigma.At[i+1][k] + mu.At[i+start][k+start]*v.At[i]
				}
				c.At[k] = -sigma.At[k+1][k]
				v.At[k] = math.Round(c.At[k])
				w.At[k] = 1
			}
		} else {
			k++
			if k == n {
				if !hasSolution {
					vec.SetZero(coeffVector)
				}
				return coeffVector
			} else {
				r.At[k] = float64(k)
				if k >= lastNonzero {
					lastNonzero = k
					v.At[k]++
				} else {
					if v.At[k] > c.At[k] {
						v.At[k] -= w.At[k]
					} else {
						v.At[k] += w.At[k]
					}
					w.At[k]++
				}
			}
		}
	}
}
