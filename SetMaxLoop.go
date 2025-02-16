package glatl

// SetMaxLoop sets limit of loops in BKZ
//
// unset limit of loops if n is not positive
func SetMaxLoop(n int) {
	maxLoop = n
}
