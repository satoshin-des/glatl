package glatl

// partDeepLLL computes partially DeepLLL-reduced basis
func partDeepLLL(b Lattice, delta float64, start int, n int) {
	if delta < 0.25 || delta > 1 {
		panic("reduction parameter must be in [1/4, 1]")
	}

	var temp float64
	gsoB, mu := GSO(b)

	for k := start; k < n; {
		for j := k - 1; j > -1; j-- {
			PartSize(b, mu, k, j)
		}

		temp = 0
		for j := 0; j < b.NumCols; j++ {
			temp += float64(b.Basis[k][j] * b.Basis[k][j])
		}

		for i := 0; i < k; {
			if temp >= delta*gsoB.At[i] {
				temp -= mu.At[k][i] * mu.At[k][i] * gsoB.At[i]
				i++
			} else {
				DeepIns(b, i, k)
				updateDeepInsGSO(i, k, gsoB, mu, b)

				if i < 1 {
					k = 0
				} else {
					k = i - 1
				}
			}
		}
		k++
	}
}
