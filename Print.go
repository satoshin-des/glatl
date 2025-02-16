package glatl

import "fmt"

// Print prints lattice basis matrix
func Print(b Lattice) {
	fmt.Printf("[\n")
	for i := 0; i < b.NumRows; i++ {
		fmt.Printf("[%d", b.Basis[i][0])
		for j := 1; j < b.NumCols; j++ {
			fmt.Printf(", %d", b.Basis[i][j])
		}
		fmt.Printf("],\n")
	}
	fmt.Printf("]\n")
}
