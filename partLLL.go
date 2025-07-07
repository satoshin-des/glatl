package glatl

// partLLL computes partially LLL-reduced basis
func partLLL(b Lattice, delta float64, start int, n int) {
	var tempInt int64
	var tempFloat, nu float64

	gsoB, mu := GSO(b)

	for k := start; k < n; {
		for j := k - 1; j > -1; j-- {
			PartSize(b, mu, k, j)
		}

		if k > 0 && gsoB.At[k] < (delta-mu.At[k][k-1]*mu.At[k][k-1])*gsoB.At[k-1] {
			for j := 0; j < b.NumCols; j++ {
				tempInt = b.Basis[k][j]
				b.Basis[k][j] = b.Basis[k-1][j]
				b.Basis[k-1][j] = tempInt
			}

			nu = mu.At[k][k-1]
			tempFloat = gsoB.At[k] + nu*nu*gsoB.At[k-1]
			mu.At[k][k-1] = nu * gsoB.At[k-1] / tempFloat
			gsoB.At[k] *= gsoB.At[k-1] / tempFloat
			gsoB.At[k-1] = tempFloat

			for j := 0; j < k-1; j++ {
				tempFloat = mu.At[k][j]
				mu.At[k][j] = mu.At[k-1][j]
				mu.At[k-1][j] = tempFloat
			}

			for i := k + 1; i < n; i++ {
				tempFloat = mu.At[i][k]
				mu.At[i][k] = mu.At[i][k-1] - nu*tempFloat
				mu.At[i][k-1] = tempFloat + mu.At[k][k-1]*mu.At[i][k]
			}

			k--
		} else {
			k++
		}
	}
}
