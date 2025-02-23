package glatl

import (
	"github.com/satoshin-des/glal/mat"
	"github.com/satoshin-des/glal/vec"
)

func updateDualDeepInsGSO(k int, l int, gsoB vec.Vector, mu mat.Matrix, dualMu mat.Matrix, dualD vec.Vector) {
	n := mu.NumRows
	var sum float64
	xi := mat.ZeroMat(mu.NumRows, mu.NumCols)
	for i := 0; i < mu.NumRows; i++ {
		for j := 0; j < mu.NumCols; j++ {
			xi.At[i][j] = mu.At[i][j]
		}
	}

	for i := l + 1; i < n; i++ {
		sum = 0
		for h := k; h <= l; h++ {
			sum += dualMu.At[k][h] * mu.At[i][h]
		}
		xi.At[i][l] = sum
	}

	for j := k; j < l; j++ {
		for i := j + 1; i < l; i++ {
			sum = 0
			for h := k; h <= j; h++ {
				sum += dualMu.At[k][h] * mu.At[i+1][h]
			}

			xi.At[i][j] = mu.At[i+1][j+1]*dualD.At[j]/dualD.At[j+1] - dualMu.At[k][j+1]/(dualD.At[j+1]*gsoB.At[j+1])*sum
		}
		xi.At[l][j] = -dualMu.At[k][j+1] / (dualD.At[j+1] * gsoB.At[j+1])
		for i := l + 1; i < n; i++ {
			sum = 0
			for h := k; h <= j; h++ {
				sum += dualMu.At[k][h] * mu.At[i][h]
			}

			xi.At[i][j] = mu.At[i][j+1]*dualD.At[j]/dualD.At[j+1] - dualMu.At[k][j+1]/(dualD.At[j+1]*gsoB.At[j+1])*sum
		}
	}

	for j := 0; j < k; j++ {
		for i := k; i < l; i++ {
			xi.At[i][j] = mu.At[i+1][j]
		}
		xi.At[l][j] = mu.At[k][j]
	}

	for i := 0; i < mu.NumRows; i++ {
		for j := 0; j < mu.NumCols; j++ {
			mu.At[i][j] = xi.At[i][j]
		}
	}

	for j := k; j < l; j++ {
		gsoB.At[j] = dualD.At[j+1] * gsoB.At[j+1] / dualD.At[j]
	}
	gsoB.At[l] = 1.0 / dualD.At[l]
}
