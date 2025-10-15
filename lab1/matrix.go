package main

import (
	"fmt"
	"math"
	//"math"
)

// SquareMatrix структура для представления квадратной матрицы
type SquareMatrix struct {
	data [][]float64
}

// NewSquareMatrixFromData создает новую квадратную матрицу из двумерного массива данных
func NewSquareMatrixFromData(data [][]float64) *SquareMatrix {
	n := len(data)

	// Создаем глубокую копию данных для избежания побочных эффектов
	new_data := make([][]float64, n)
	for i := range new_data {
		new_data[i] = make([]float64, n)
		copy(new_data[i], data[i])
	}

	return &SquareMatrix{data: new_data}
}

// NewSquareMatrix создает новую квадратную матрицу заданного размера, инициализированную нулями
func NewSquareMatrix(deg int) *SquareMatrix {
	new_data := make([][]float64, deg)
	for i := range new_data {
		new_data[i] = make([]float64, deg)
	}
	return &SquareMatrix{data: new_data}
}

// GetCol возвращает указанный столбец матрицы как вектор
func (m *SquareMatrix) GetCol(idx int) *Vector {
	res := NewVector(len(m.data))

	// Копируем элементы столбца в вектор
	for i := range m.data {
		res.data[i] = m.data[i][idx]
	}

	return res
}

// Add выполняет поэлементное сложение двух матриц
func (m *SquareMatrix) Add(other *SquareMatrix) *SquareMatrix {
	res := NewSquareMatrix(len(m.data))
	for i := range res.data {
		for j := range res.data[i] {
			res.data[i][j] = m.data[i][j] + other.data[i][j]
		}
	}
	return res
}

// MultMatrix выполняет матричное умножение текущей матрицы на другую матрицу
func (m *SquareMatrix) MultMatrix(other *SquareMatrix) *SquareMatrix {
	res := NewSquareMatrix(len(m.data))

	// Выполняем стандартное матричное умножение
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

// MultAlpha умножает все элементы матрицы на скаляр
func (m *SquareMatrix) MultAlpha(a float64) *SquareMatrix {
	res := NewSquareMatrixFromData(m.data)
	for i := range res.data {
		for j := range res.data[i] {
			res.data[i][j] *= a
		}
	}
	return res
}

// Transpose возвращает транспонированную матрицу
func (m *SquareMatrix) Transpose() *SquareMatrix {
	res := NewSquareMatrix(len(m.data))

	// Меняем местами индексы строк и столбцов
	for i := range res.data {
		for j := range res.data[i] {
			res.data[i][j] = m.data[j][i]
		}
	}
	return res
}

// OperateVector умножает матрицу на вектор
func (m *SquareMatrix) OperateVector(v *Vector) *Vector {
	res := NewVector(len(v.data))

	// Выполняем матрично-векторное умножение
	for i := range m.data {
		sum := 0.0
		for j := range m.data[i] {
			sum += m.data[i][j] * v.data[j]
		}
		res.data[i] = sum
	}
	return res
}

// ExtractMatrix извлекает подматрицу из текущей матрицы по заданным границам
func (m *SquareMatrix) ExtractMatrix(xStart, xEnd, yStart, yEnd int) *SquareMatrix {
	res := NewSquareMatrix(xEnd + 1 - xStart)
	// Копируем указанную область в новую матрицу
	for i := range res.data {
		for j := range res.data[i] {
			res.data[i][j] = m.data[yStart+i][xStart+j]
		}
	}
	return res
}

// CompleteMatrix дополняет текущую матрицу до большего размера, добавляя единичную матрицу
func (m *SquareMatrix) CompleteMatrix(deg int) *SquareMatrix {
	res := NewSquareMatrix(deg)
	start := deg - len(m.data)
	for i := range res.data {
		for j := range res.data[i] {
			if i >= start && j >= start {
				// Копируем элементы из исходной матрицы
				res.data[i][j] = m.data[i-start][j-start]
			} else if i == j {
				// Заполняем диагональ единицами
				res.data[i][j] = 1
			}
		}
	}
	return res
}

// MaxRowNorm вычисляет максимальную строчную норму матрицы (норму бесконечности)
func (m *SquareMatrix) MaxRowNorm() float64 {
	var max float64 = 0.0
	for i := range m.data {
		sum := 0.0
		// Суммируем абсолютные значения элементов строки
		for j := range m.data[i] {
			sum += math.Abs(m.data[i][j])
		}
		// Обновляем максимум если текущая сумма больше
		if sum > max {
			max = sum
		}
	}
	return max
}

// ToString возвращает строковое представление матрицы для вывода
func (m *SquareMatrix) ToString() string {
	res := ""
	for i := range m.data {
		for j := range m.data[i] {
			if m.data[i][j] == 0 {
				// Форматируем нули без десятичных знаков
				res += fmt.Sprintf("%12d", 0)
			} else {
				// Форматируем ненулевые элементы с 5 знаками после запятой
				res += fmt.Sprintf("%12.5f", m.data[i][j])
			}
		}
		res += "\n"
	}
	return res
}