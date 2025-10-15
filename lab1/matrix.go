package main

import (
	"fmt"
	"math"
	//"math"
)

type SquareMatrix struct {
	data [][]float64
}

func NewSquareMatrixFromData(data [][]float64) *SquareMatrix {
	n := len(data)

	new_data := make([][]float64, n)

	for i := range new_data {
		new_data[i] = make([]float64, n)
		copy(new_data[i], data[i])
	}

	return &SquareMatrix{data: new_data}
}

func NewSquareMatrix(deg int) *SquareMatrix {
	new_data := make([][]float64, deg)
	for i := range(new_data){
		new_data[i] = make([]float64, deg)
	}
	return &SquareMatrix{data: new_data}
}

func (m *SquareMatrix) GetCol(idx int) *Vector{
	res := NewVector(len(m.data))

	for i := range(m.data){
		res.data[i] = m.data[i][idx]
	}

	return res
}

func (m *SquareMatrix) Add(other *SquareMatrix) *SquareMatrix {
	res := NewSquareMatrix(len(m.data))
	for i := range res.data {
		for j := range res.data[i] {
			res.data[i][j] = m.data[i][j] +  other.data[i][j]
		}
	}
	return res
}

func (m *SquareMatrix) MultMatrix(other *SquareMatrix) *SquareMatrix {
	res := NewSquareMatrix(len(m.data))

	for i := range res.data {
		for j := range res.data[i] {
			var sum float64
			for k := range m.data {
				sum += m.data[i][k] * other.data[k][j]
			}
			res.data[i][j] = sum
		}
	}

	return res
}

func (m *SquareMatrix) MultAlpha(a float64) *SquareMatrix {
	res := NewSquareMatrixFromData(m.data)
	for i := range(res.data) {
		for j := range(res.data[i]){
			res.data[i][j] *= a 
		}
	}
	return res
} 

func (m *SquareMatrix) Transpose() *SquareMatrix {
	res := NewSquareMatrix(len(m.data))

	for i := range(res.data){
		for j := range(res.data[i]){
			res.data[i][j] = m.data[j][i]
		}
	}
	return res
}

func (m *SquareMatrix) OperateVector(v *Vector) *Vector {
	res := NewVector(len(v.data))

	for i := range(m.data) {
		sum := 0.0
		for j := range(m.data[i]){
			sum += m.data[i][j] * v.data[j]
		}
		res.data[i] = sum
	}
	//res.ClearZeros()
	return res
}

func (m *SquareMatrix) ExtractMatrix(xStart, xEnd, yStart, yEnd int) *SquareMatrix{
	res := NewSquareMatrix(xEnd + 1 - xStart)
	for i:= range(res.data){
		for j := range(res.data[i]){
			res.data[i][j] = m.data[yStart + i][xStart + j]
		}
	}
	return res
}

func (m *SquareMatrix) CompleteMatrix(deg int) *SquareMatrix {
	res := NewSquareMatrix(deg)
	start := deg - len(m.data)
	for i := range(res.data){
		for j := range(res.data[i]){
			if i >= start && j >= start {
				res.data[i][j] = m.data[i - start][j - start]
			}else if i == j{
				res.data[i][j] = 1
			}
		}
	}
	return res
}

func (m *SquareMatrix) MaxRowNorm() float64 {
	var max float64 = 0.0
	for i:= range(m.data) {
		sum := 0.0
		for j := range(m.data[i]){
			sum += math.Abs(m.data[i][j])
		}
		if sum > max {
			max = sum
		}
	}
	return max
}

// func (m *SquareMatrix) ClearZeros() {
// 	for i := range(m.data){
// 		for j := range(m.data[i]){
// 			if math.Abs(m.data[i][j]) < 1e-10{
// 				m.data[i][j] = 0
// 			}
// 		}
// 	}
// } //! CLEAR ZEROS

func (m *SquareMatrix) ToString() string {
	
	res := ""
	for i := range (m.data) {
		for j := range( m.data[i]) {
			if m.data[i][j] == 0 {
				res += fmt.Sprintf("%12d", 0)
			}else{
				res += fmt.Sprintf("%12.5f", m.data[i][j])
			}
		}
		res += "\n"
	}
	return res
}