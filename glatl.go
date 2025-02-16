package glatl

// Lattice is lattice basis matrix type
type Lattice struct {
	NumRows int
	NumCols int
	Basis   [][]int64
}

var maxLoop int = 10
