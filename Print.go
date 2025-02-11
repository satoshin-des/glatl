package glatl

import "fmt"

// Print prints lattice basis matrix
func Print(b Lattice) {
	for i := 0; i < b.NumRows; i++ {
		for j := 0; j < b.NumCols; j++ {
			fmt.Printf("%d ", b.Basis[i][j])
		}
		fmt.Println()
	}
}
