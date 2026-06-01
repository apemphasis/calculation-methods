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

func NewEMatrix(deg int) *SquareMatrix {
	new_data := make([][]float64, deg)
	for i := range new_data {
		new_data[i] = make([]float64, deg)
		new_data[i][i] = 1
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
func (m *SquareMatrix) String() string {
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

// CountW вычисляет вектор отражения W для преобразования Хаусхолдера
func CountW(s *Vector) *Vector {
	// Создаем вектор E того же размера что и s
	e_data := make([]float64, len(s.data))
	// Если первый элемент не нулевой, вычисляем E с учетом знака
	if s.data[0] != 0 {
		e_data[0] = (s.data[0] / math.Abs(s.data[0])) * s.Abs()
	} else {
		// Если первый элемент нулевой, просто используем норму вектора
		e_data[0] = s.Abs()
	}
	// Создаем вектор E из вычисленных данных
	E := NewVectorFromData(e_data)
	// Вычисляем W = s + E
	W := s.Add(E)
	// Нормализуем W (делим на его норму)
	W = W.MultAlpha(1 / W.Abs())

	return W
}

// CountH вычисляет матрицу отражения Хаусхолдера по вектору W
func CountH(w *Vector) *SquareMatrix {
	// Создаем единичную матрицу
	res := NewSquareMatrix(len(w.data))
	for i := range res.data {
		res.data[i][i] = 1
	}

	// Вычисляем внешнее произведение w * w^T
	wwt := NewSquareMatrix(len(w.data))
	for i := range wwt.data {
		for j := range wwt.data[i] {
			wwt.data[i][j] = w.data[i] * w.data[j]
		}
	}

	// Вычисляем H = I - 2 * w * w^T
	res = res.Add(wwt.MultAlpha(-2))
	return res
}

// CountQ вычисляет матрицу преобразования Хаусхолдера для первого столбца
func (m *SquareMatrix) CountQ() *SquareMatrix {
	// Получаем первый столбец матрицы
	W := CountW(m.GetCol(0))
	// Вычисляем матрицу отражения
	return CountH(W)
}

// FindQR выполняет QR разложение методом отражений Хаусхолдера
func (m *SquareMatrix) FindQR() (*SquareMatrix, *SquareMatrix) {
	deg := len(m.data)
	// Инициализируем R как копию исходной матрицы
	R := NewSquareMatrixFromData(m.data)
	// Массив для хранения матриц преобразований
	QArr := make([]*SquareMatrix, 0, deg-1)
	
	// Последовательно обнуляем поддиагональные элементы
	for i := range deg - 1 {
		// Извлекаем подматрицу для текущего шага
		H := R.ExtractMatrix(i, deg-1, i, deg-1).CountQ()
		// Дополняем до полной размерности
		Qi := H.CompleteMatrix(deg)
		
		// Применяем преобразование к R
		R = Qi.MultMatrix(R)
		// Сохраняем матрицу преобразования
		QArr = append(QArr, Qi)
	}

	// Вычисляем полную матрицу Q как произведение преобразований
	Q := NewSquareMatrixFromData(QArr[len(QArr)-1].data)
	for i := len(QArr) - 2; i >= 0; i-- {
		Q = QArr[i].MultMatrix(Q)
	}

	return Q, R
}

// SolveByQR решает систему уравнений методом QR разложения
func (m *SquareMatrix) SolveByQR(b *Vector) *Vector {
	deg := len(m.data)
	R := NewSquareMatrixFromData(m.data)
	bi := NewVectorFromData(b.data)
	
	// Выполняем QR разложение с одновременным преобразованием правой части
	for i := range deg - 1 {
		H := R.ExtractMatrix(i, deg-1, i, deg-1).CountQ()
		Qi := H.CompleteMatrix(deg)
		
		// Применяем преобразование к матрице и правой части
		R = Qi.MultMatrix(R)
		bi = Qi.OperateVector(bi)
	}

	// Решаем верхнюю треугольную систему Rx = bi обратной подстановкой
	solution := GaussReverse(R, bi)

	return solution
}