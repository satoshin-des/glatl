package glatl

import (
	"math"

	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

// ENUM returns a lattice vector whose norm is lesser than or equal to R on lattice spaned by b
func ENUM(mu mat.Matrix, gsoB vec.Vector, R float64, rho vec.Vector) vec.Vector {
	n := mu.NumRows
	var temp float64
	var lastNonzero int = 0

	r := vec.ZeroVec(n + 1)
	w := vec.ZeroVec(n)
	v := vec.ZeroVec(n)
	c := vec.ZeroVec(n)
	sigma := mat.ZeroMat(n+1, n)
	vec.SetZero(rho)

	v.At[0] = 1
	for i := 0; i < n; i++ {
		r.At[i] = float64(i)
	}

	for k := 0; ; {
		temp = float64(v.At[k]) - c.At[k]
		temp *= temp
		rho.At[k] = rho.At[k+1] + temp*gsoB.At[k]
		if rho.At[k] <= R {
			if k == 0 {
				return v
			} else {
				k--
				if r.At[k] <= r.At[k+1] {
					r.At[k] = r.At[k+1]
				}
				for i := int(r.At[k]); i > k; i-- {
					sigma.At[i][k] = sigma.At[i+1][k] + mu.At[i][k]*v.At[i]
				}
				c.At[k] = -sigma.At[k+1][k]
				v.At[k] = math.Round(c.At[k])
				w.At[k] = 1
			}
		} else {
			k++
			if k == n {
				vec.SetZero(v)
				return v
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
