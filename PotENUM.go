package glatl

import (
	"math"

	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

func PotENUM(mu mat.Matrix, gsoB vec.Vector, logB vec.Vector) vec.Vector {
	var lastNonzero int = 0
	var temp float64
	var P float64 = 0
	var R float64 = logB.At[0]
	n := mu.NumRows
	r := vec.ZeroVec(n + 1)
	w := vec.ZeroVec(n)
	v := vec.ZeroVec(n)
	c := vec.ZeroVec(n)
	d := vec.ZeroVec(n + 1)
	sigma := mat.ZeroMat(n+1, n)

	v.At[0] = 1

	for i := 0; i <= n; i++ {
		r.At[i] = float64(i)
	}

	for k := 0; ; {
		temp = float64(v.At[k]) - c.At[k]
		temp *= temp
		d.At[k] = d.At[k+1] + temp*gsoB.At[k]

		if float64(k+1)*math.Log(d.At[k])+P < float64(k+1)*math.Log(0.99)+R {
			if k == 0 {
				return v
			} else {
				P += math.Log(d.At[k])
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
				r.At[k-1] = float64(k)
				if k >= lastNonzero {
					lastNonzero = k
					v.At[k]++
					P = 0
					R = 0
					for i := 0; i <= lastNonzero; i++ {
						R += logB.At[i]
					}
				} else {
					if v.At[k] > c.At[k] {
						v.At[k] -= w.At[k]
					} else {
						v.At[k] += w.At[k]
					}
					w.At[k]++
					P -= math.Log(d.At[k])
				}
			}
		}
	}
}
