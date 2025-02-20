package glatl

// PotLLL computes delta-PotLLL reduced basis
//
// panic if delta < 1/4 or delta > 1
//
// F. Fontein, M. Schneider, U. Wagner. PotLLL: A polynomial time version of LLL with deep insertions.(2014)
func PotLLL(b Lattice, delta float64) {
	if delta < 0.25 || delta > 1 {
		panic("reduction parameter must be in [1/4, 1]")
	}

	LLL(b, delta)

	gsoB, mu := GSO(b)

	var p, pMin, sum float64
	var k int

	for l := 0; l < b.NumRows; {
		for j := l - 1; j > -1; j-- {
			PartSize(b, mu, l, j)
		}

		p = 1
		pMin = 1
		k = 0

		for j := l - 1; j > -1; j-- {
			sum = 0
			for i := j; i < l; i++ {
				sum += mu.At[l][i] * mu.At[l][i] * gsoB.At[i]
			}
			p *= (gsoB.At[l] + sum) / gsoB.At[j]

			if p < pMin {
				k = j
				pMin = p
			}
		}

		if delta > pMin {
			DeepIns(b, k, l)

			updateDeepInsGSO(k, l, gsoB, mu, b)

			l = k
		} else {
			l++
		}
	}
}
