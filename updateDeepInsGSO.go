package glatl

import (
	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

func updateDeepInsGSO(i int, k int, gsoB vec.Vector, mu mat.Matrix, b Lattice) {
	var T float64
	var eps float64
	var P vec.Vector = vec.ZeroVec(b.NumRows)
	var D vec.Vector = vec.ZeroVec(b.NumRows)
	var S vec.Vector = vec.ZeroVec(b.NumRows)

	vec.SetZero(P)
	vec.SetZero(D)
	vec.SetZero(S)

	P.At[k] = gsoB.At[k]
	D.At[k] = gsoB.At[k]
	for j := k - 1; j >= i; j-- {
		P.At[j] = mu.At[k][j] * gsoB.At[j]
		D.At[j] = D.At[j+1] + mu.At[k][j]*P.At[j]
	}

	for j := k; j > i; j-- {
		T = mu.At[k][j-1] / D.At[j]
		for l := b.NumRows - 1; l > k; l-- {
			S.At[l] += mu.At[l][j] * P.At[j]
			mu.At[l][j] = mu.At[l][j-1] - T*S.At[l]
		}
		for l := k; l > j; l-- {
			S.At[l] += mu.At[l-1][j] * P.At[j]
			mu.At[l][j] = mu.At[l-1][j-1] - T*S.At[l]
		}
	}

	T = 1.0 / D.At[i]

	for l := b.NumRows - 1; l > k; l-- {
		mu.At[l][i] = T * (S.At[l] + mu.At[l][i]*P.At[i])
	}
	for l := k; l >= i+2; l-- {
		mu.At[l][i] = T * (S.At[l] + mu.At[l-1][i]*P.At[i])
	}

	mu.At[i+1][i] = T * P.At[i]
	for j := 0; j < i; j++ {
		eps = mu.At[k][j]
		for l := k; l > i; l-- {
			mu.At[l][j] = mu.At[l-1][j]
		}
		mu.At[i][j] = eps
	}

	for j := k; j > i; j-- {
		gsoB.At[j] = D.At[j] * gsoB.At[j-1] / D.At[j-1]
	}
	gsoB.At[i] = D.At[i]
}
