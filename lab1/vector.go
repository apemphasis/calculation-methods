package main

import (
	"fmt"
	"math"
)

type Vector struct {
	data []float64
}

func NewVector(deg int) *Vector{
	data := make([]float64, deg)
	return &Vector{data: data}
}

func NewVectorFromData(data []float64) *Vector{
	new_data := make([]float64, len(data))
	copy(new_data, data)
	return &Vector{data: new_data}
}

func (v *Vector) Abs() float64 {
	var res float64

	for _, el := range v.data {
		res += el * el
	}

	return math.Sqrt(res)
}

func (v *Vector) Add(other * Vector) *Vector {
	res := NewVector(len(v.data))
	for i := range(res.data){
		res.data[i] = v.data[i] + other.data[i]
	}
	return res
}

func (v *Vector) MultAlpha(a float64) *Vector {
	res := NewVectorFromData(v.data)
	for i := range(res.data) {
		res.data[i] *= a 
	}
	return res
} 

// func (v *Vector) ClearZeros() {
// 	for i := range(v.data){
// 		if math.Abs(v.data[i]) < 1e-10{
// 			v.data[i] = 0
// 		}
// 	}
// }

func (v *Vector) ToString() string {
	res := ""
	for i := range(v.data){
		if v.data[i] == 0 {
				res += fmt.Sprintf("%12d", 0)
		}else{
			res += fmt.Sprintf("%12.5f", v.data[i])
		}
	}

	return res
}