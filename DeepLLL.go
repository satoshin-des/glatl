package glatl

func DeepLLL(b Lattice, delta float64) {
	var temp float64
	var t int64
	gsoB, mu := GSO(b)

	for k := 1; k < b.NumRows; {
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
				for l := 0; l < b.NumCols; l++ {
					t = b.Basis[k][l]
					for j := k; j > i; j-- {
						b.Basis[j][l] = b.Basis[j-1][l]
					}
					b.Basis[i][l] = t
				}

				gsoB, mu = GSO(b)

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
