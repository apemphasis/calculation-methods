package main //Kozlovsky_13gr_5v

import (
	"fmt"
	"math"
)

type SquareMatrix struct{
	data [][]float64
}

func NewSquareMatrix(data [][]float64) SquareMatrix{
	n := len(data)

	new_data := make([][]float64, n)

	for i := range(new_data){
		new_data[i] = make([]float64, n)
		copy(new_data[i], data[i])
	}

	return SquareMatrix{data: new_data}
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

	matrix := SquareMatrix{data: make([][]float64, 15)} 
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

	m1 := NewSquareMatrix([][]float64 {{1.0, 4, 2},{13, 5, 3}, {0.2, 1, 5.9}})
	m2 := NewSquareMatrix([][]float64 {{1.0, 4, 2},{13, 5, 3}, {0.2, 1, 5.9}})


	fmt.Println(m1.Add(m2).String())

	// метод отражений

}