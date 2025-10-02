package main //Kozlovsky_13gr_5v

import (
	"fmt"
	"math"
)

func binarySearch(arr []int, n int) int {

	l := 0
	r := len(arr)
	
	for l != r {
		k := (l + r) / 2

		if arr[k] == n {
			return k
		}

		if arr[k] > n {
			r = k
		} else {
			l = k + 1
		}

	}
	return l
}

type Matrix interface{
	Add(other Matrix) Matrix
	Substract(other Matrix) Matrix
	Mult(other Matrix) Matrix
	Transpose() Matrix
	String() string
	Copy() Matrix
}

type SquareMatrix struct{
	size int
	data [][]float64
}

// func (m SquareMatrix) Add(other SquareMatrix) SquareMatrix{

// }



func main() {

	matrix := make([][]float64, 15)
	for i := range(matrix) {
		matrix[i] = make([]float64, 15)
		for j := range(matrix[i]){
			if i == j{
				matrix[i][j] = 5 * math.Sqrt(float64(i+1))
				continue
			}
			if i == j + 1{
				matrix[i][j] = -0.01 * (math.Sqrt(float64(i+1)) + float64(j+1))
			}
		}
	}

	for _, row := range(matrix){
		for _, el := range(row){
			fmt.Printf("%10.5f",el)
		}
		fmt.Println()
	}

	// метод отражений


}