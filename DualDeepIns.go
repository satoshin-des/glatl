package glatl

func DualDeepIns(b Lattice, k int, l int) {
	var t int64

	for j := 0; j < b.NumCols; j++ {
		t = b.Basis[k][j]
		for i := k; i < l; i++ {
			b.Basis[i][j] = b.Basis[i+1][j]
		}
		b.Basis[l][j] = t
	}
}
