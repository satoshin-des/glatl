package glatl

// HKZ computes "HKZ-reduced" basis without MLLL to abord faster
//
// panic if delta < 1/4 or delta > 1
//
// C. P. Schnorr, M. Euchner. Lattice basis reduction: Improved practical algorithms and solving subset sum problem.(1994)
// A. Korkine, G. Zolotareff. Sur les formes quadratiques positives ternaires. (1872)
// A. Korkine, G. Zolotareff. Sur les formes quadratiques. (1873)
func HKZ(b Lattice, delta float64) {
	BKZ(b, b.NumRows, delta)
}
