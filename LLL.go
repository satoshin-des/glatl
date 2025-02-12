package glatl

// LLL computes delta-LLL-reduced basis
//
// panic if delta < 1/4 or delta > 1
//
// A. K. Lenstra, H. W. Lenstra, L. Lovasz. Factoring polynomials with rational coefficients. 1982
func LLL(b Lattice, delta float64) {
	if delta < 0.25 || delta > 1 {
		panic("reduction parameter must be in [1/4, 1]")
	}

	var tempInt int64
	var tempFloat, nu float64

	gsoB, mu := GSO(b)

	for k := 1; k < b.NumRows; {
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

			for i := k + 1; i < b.NumRows; i++ {
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
