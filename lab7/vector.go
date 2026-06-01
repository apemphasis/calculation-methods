package main

import (
	"fmt"
	"math"
)

// Vector структура для представления вектора
type Vector struct {
	data []float64
}

// NewVector создает новый вектор заданной размерности, инициализированный нулями
func NewVector(deg int) *Vector {
	data := make([]float64, deg)
	return &Vector{data: data}
}

// NewVectorFromData создает новый вектор из массива данных (с глубоким копированием)
func NewVectorFromData(data []float64) *Vector {
	new_data := make([]float64, len(data))
	copy(new_data, data)
	return &Vector{data: new_data}
}

// Abs вычисляет евклидову норму (длину) вектора
func (v *Vector) Abs() float64 {
	var res float64

	// Суммируем квадраты всех элементов вектора
	for _, el := range v.data {
		res += el * el
	}

	// Возвращаем квадратный корень из суммы квадратов
	return math.Sqrt(res)
}

// Add выполняет поэлементное сложение двух векторов
func (v *Vector) Add(other *Vector) *Vector {
	res := NewVector(len(v.data))
	for i := range res.data {
		res.data[i] = v.data[i] + other.data[i]
	}
	return res
}

// MultAlpha умножает все элементы вектора на скаляр
func (v *Vector) MultAlpha(a float64) *Vector {
	res := NewVectorFromData(v.data)
	for i := range res.data {
		res.data[i] *= a
	}
	return res
}

// ToString возвращает строковое представление вектора для вывода
func (v *Vector) String() string {
	res := ""
	for i := range v.data {
		if v.data[i] == 0 {
			// Форматируем нули без десятичных знаков
			res += fmt.Sprintf("%12d", 0)
		} else {
			// Форматируем ненулевые элементы с 5 знаками после запятой
			res += fmt.Sprintf("%12.5f", v.data[i])
		}
	}

	return res
}

func (v * Vector) CubicNorm() float64 {
	res := 0.0
	for _, el := range v.data { 
		if res < math.Abs(el) {
			res = math.Abs(el)
		}
		
	}
	return res
}

func (v *Vector) Normalize() *Vector {
	res := make([]float64, len(v.data))
	sum := 0.0
	for _, el := range v.data {
		sum += el*el
	}
	for i := range res {
		res[i] = v.data[i] / math.Sqrt(sum)
	}

	return NewVectorFromData(res)
}