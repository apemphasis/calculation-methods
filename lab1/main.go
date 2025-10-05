package main //Kozlovsky_13gr_5v

import (
	"fmt"
	"math"
)

type Matrix interface{
	Add(other Matrix) Matrix
	Substract(other Matrix) Matrix
	Mult(other Matrix) Matrix
	Transpose() Matrix
	String() string
	Copy() Matrix
}

type SquareMatrix struct{
	data [][]float64
}

func (m SquareMatrix) Add(other SquareMatrix) SquareMatrix{
	res := m
	for i := range(res.data){
		for j:= range(res.data[i]){
			res.data[i][j] += other.data[i][j]
		}
	}
	return res
}

func (m SquareMatrix) String() string{
	res := ""
	for i := range(m.data){
		for j := range(m.data[i]){
			res += fmt.Sprintf("%10.5f", m.data[i][j])
		}
		res += "\n"
	}
	return res
}

func main() {

	matrix := &SquareMatrix{data: make([][]float64, 15)} 
	for i := range(matrix.data) {
		matrix.data[i] = make([]float64, 15)
		for j := range(matrix.data[i]){
			if i == j{
				matrix.data[i][j] = 5 * math.Sqrt(float64(i+1))
				continue
			}
			if i == j + 1{
				matrix.data[i][j] = -0.01 * (math.Sqrt(float64(i+1)) + float64(j+1))
			}
		}
	}

	fmt.Println(matrix.String())

	// метод отражений


}