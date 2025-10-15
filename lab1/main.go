package main //Козловский Евгений, 13 группа, 5 вариант

import (
	"fmt"
	"math"
)

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

// FindLDLt выполняет LDL^T разложение симметричной матрицы
func (m *SquareMatrix) FindLDLt() (*SquareMatrix, *SquareMatrix, *SquareMatrix) {
	// Создаем копию исходной матрицы
	A_copy := NewSquareMatrixFromData(m.data)
	// Создаем матрицы для диагональной D и нижней треугольной L
	D := NewSquareMatrix(len(m.data))
	L := NewSquareMatrix(len(m.data))

	// Инициализируем L как единичную матрицу
	for i := range L.data {
		L.data[i][i] = 1
	}

	// Выполняем разложение методом Холецкого
	for k := 0; k < len(m.data)-1; k++ {
		D.data[k][k] = A_copy.data[k][k]
		for i := k + 1; i < len(m.data); i++ {
			L.data[i][k] = A_copy.data[i][k] / D.data[k][k]
			for j := k + 1; j <= i; j++ {
				A_copy.data[i][j] -= L.data[i][k] * A_copy.data[j][k]
			}
		}
	}

	// Устанавливаем последний диагональный элемент
	D.data[len(m.data)-1][len(m.data)-1] = A_copy.data[len(m.data)-1][len(m.data)-1]
	// Транспонируем L чтобы получить L^T
	Lt := L.Transpose()

	return L, D, Lt
}

// SolveByLDLt решает систему уравнений используя LDL^T разложение
func (m *SquareMatrix) SolveByLDLt(b *Vector) *Vector {
	// Получаем LDL^T разложение
	L, D, Lt := m.FindLDLt()

	// Решаем Lz = b прямым ходом Гаусса
	z := Gauss(L, b)
	// Решаем Dy = z обратным ходом Гаусса
	y := GaussReverse(D, z)
	// Решаем L^Tx = y обратным ходом Гаусса
	x := GaussReverse(Lt, y)

	return x
}

// FindInverse находит обратную матрицу с помощью LDL^T разложения
func (m *SquareMatrix) FindInverse() *SquareMatrix {
	// Создаем матрицу для результата
	res := NewSquareMatrix(len(m.data))

	// Для каждого столбца единичной матрицы решаем систему
	for i := range res.data {
		// Создаем i-ый базисный вектор
		e := NewVector(len(res.data))
		e.data[i] = 1
		// Решаем систему Ax = e
		col := m.SolveByLDLt(e)
		// Записываем решение в i-ый столбец обратной матрицы
		for j := range col.data {
			res.data[j][i] = col.data[j]
		}
	}

	return res
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

func main() {
	// Создаем матрицу A размером 15x15
	A_data := make([][]float64, 15)
	for i := range A_data {
		A_data[i] = make([]float64, 15)
		for j := range A_data[i] {
			if i == j {
				// Диагональные элементы
				A_data[i][j] = 5 * math.Sqrt(float64(i+1))
				continue
			} else {
				// Недиагональные элементы
				A_data[i][j] = math.Sqrt(float64(i+1)) + math.Sqrt(float64(j+1))
			}
		}
	}

	// Создаем вектор правой части b
	b_data := make([]float64, 15)
	for i := range b_data {
		b_data[i] = 4.5 * math.Sqrt(float64(i+1))
	}

	fmt.Printf("\n==========================  ВХОДНЫЕ ДАННЫЕ  ==========================\n\n")

	A := NewSquareMatrixFromData(A_data)
	b := NewVectorFromData(b_data)

	fmt.Println("A:")
	fmt.Println(A.ToString())

	fmt.Println("b:")
	fmt.Println(b.ToString())

	fmt.Printf("\n\n==========================  LDL^t РАЗЛОЖЕНИЕ  ==========================\n\n")
	
	// Выполняем LDL^T разложение
	L, D, Lt := A.FindLDLt()
	fmt.Println("L:")
	fmt.Println(L.ToString())

	fmt.Println("D:")
	fmt.Println(D.ToString())

	fmt.Println("L^t:")
	fmt.Println(Lt.ToString())

	fmt.Println()

	// Проверяем точность разложения
	fmt.Print("Норма ||A - LDL^t|| (равномерная): ")
	fmt.Println(A.Add(L.MultMatrix(D.MultMatrix(Lt)).MultAlpha(-1)).MaxRowNorm())

	fmt.Printf("\n\n==========================  РЕШЕНИЕ С ПОМОЩЬЮ LDL^t  ==========================\n\n")

	// Решаем систему с помощью LDL^T разложения
	x := A.SolveByLDLt(b)

	fmt.Println("Решение:")
	fmt.Println(x.ToString())
	fmt.Println()
	// Вычисляем норму невязки
	fmt.Print("Норма невязки (евклидова): ")
	fmt.Println(b.Add(A.OperateVector(x).MultAlpha(-1)).Abs())

	fmt.Printf("\n\n==========================  ОБРАТНАЯ С ПОМОЩЬЮ LDL^t  ==========================\n\n")

	// Находим обратную матрицу
	inverseA := A.FindInverse()

	fmt.Println("Обратная:")
	fmt.Println(inverseA.ToString())

	fmt.Println()

	// Проверяем точность обращения
	fmt.Print("Норма ||E - AA^(-1)|| (равномерная): ")

	E := NewSquareMatrix(15)
	for i := range E.data {
		E.data[i][i] = 1
	}

	fmt.Println(E.Add(A.MultMatrix(A.FindInverse()).MultAlpha(-1)).MaxRowNorm())

	fmt.Printf("\n\n==========================  QR РАЗЛОЖЕНИЕ  ==========================\n\n")

	// Выполняем QR разложение
	Q, R := A.FindQR()

	fmt.Println("Q:")
	fmt.Println(Q.ToString())

	fmt.Println("R:")
	fmt.Println(R.ToString())

	fmt.Println()

	// Проверяем точность QR разложения
	fmt.Print("Норма ||A - QR|| (равномерная): ")
	fmt.Println(A.Add(Q.MultMatrix(R).MultAlpha(-1)).MaxRowNorm())

	fmt.Printf("\n\n==========================  РЕШЕНИЕ МЕТОДОМ ОТРАЖЕНИЙ  ==========================\n\n")

	// Решаем систему методом QR разложения
	solution := A.SolveByQR(b)

	fmt.Println("Решение:")
	fmt.Println(solution.ToString())
	fmt.Println()
	// Вычисляем норму невязки
	fmt.Print("Норма невязки (евклидова): ")
	fmt.Println(b.Add(A.OperateVector(solution).MultAlpha(-1)).Abs())
	
	fmt.Println()
}