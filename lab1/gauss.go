package main

// GaussReverse решает верхнюю треугольную систему уравнений обратной подстановкой
func GaussReverse(m *SquareMatrix, b *Vector) *Vector {
	n := len(b.data)

	// Создаем вектор для решения
	solution := NewVector(n)
	// Находим последний элемент решения: x_n = b_n / a_nn
	solution.data[n-1] = b.data[n-1] / m.data[n-1][n-1]

	// Обратный ход: от предпоследнего уравнения к первому
	for i := n - 2; i >= 0; i-- {
		sum := 0.0
		// Суммируем известные члены: a_i,i+1 * x_i+1 + a_i,i+2 * x_i+2 + ... + a_i,n-1 * x_n-1
		for j := 1; i+j < n; j++ {
			sum -= m.data[i][i+j] * solution.data[i+j]
		}
		// Вычисляем x_i = (b_i - сумма) / a_ii
		solution.data[i] = (sum + b.data[i]) / m.data[i][i]
	}

	return solution
}

// Gauss решает нижнюю треугольную систему уравнений прямой подстановкой
func Gauss(m *SquareMatrix, b *Vector) *Vector {
	n := len(b.data)

	// Создаем вектор для решения
	solution := NewVector(n)
	// Находим первый элемент решения: x_0 = b_0 / a_00
	solution.data[0] = b.data[0] / m.data[0][0]

	// Прямой ход: от второго уравнения до последнего
	for i := 1; i < n; i++ {
		sum := 0.0
		// Суммируем известные члены: a_i0 * x_0 + a_i1 * x_1 + ... + a_i,i-1 * x_i-1
		for j := 0; j < i; j++ {
			sum -= m.data[i][j] * solution.data[j]
		}
		// Вычисляем x_i = (b_i - сумма) / a_ii
		solution.data[i] = (sum + b.data[i]) / m.data[i][i]
	}

	return solution
}