package glatl

// DeepIns computes deep-insertion of b with i, k
func DeepIns(b Lattice, i int, k int) {
	var t int64

	for l := 0; l < b.NumCols; l++ {
		t = b.Basis[k][l]
		for j := k; j > i; j-- {
			b.Basis[j][l] = b.Basis[j-1][l]
		}
		b.Basis[i][l] = t
	}
}
