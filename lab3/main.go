package main

import (
	"fmt"
	//"math/rand"
	"math"
)

func (m *SquareMatrix) findMax() (float64, int, int) {
	max := 0.0
	maxi, maxj := 0, 0
	for i := 0; i < len(m.data); i++ {
		for j := i+1; j < len(m.data); j++ {
			if math.Abs(m.data[i][j]) > max {
				max = math.Abs(m.data[i][j])
				maxi, maxj = i, j
			}
		}
	}

	return max, maxi, maxj
}

func JacobiRotateIteration(A *SquareMatrix, maxi, maxj int) (*SquareMatrix, *SquareMatrix){
	
	p := 2 * A.data[maxi][maxj] / (A.data[maxi][maxi] - A.data[maxj][maxj])
	cos := math.Sqrt((1 + 1/math.Sqrt(1 + p * p)) / 2)
	sin := math.Copysign(math.Sqrt((1 - 1/math.Sqrt(1 + p * p)) / 2), p)

	Q := NewEMatrix(len(A.data))
	Q.data[maxi][maxi] = cos
	Q.data[maxi][maxj] = -sin
	Q.data[maxj][maxi] = sin
	Q.data[maxj][maxj] = cos

	Qt := Q.Transpose()

	Ak := Qt.MultMatrix(A.MultMatrix(Q))

	return Ak, Q
}

func JacobiRotate(A *SquareMatrix, eps float64) (*Vector, *SquareMatrix, int) {
	max, maxi, maxj := A.findMax()
	V := NewEMatrix(len(A.data)) 
	counter := 0
	for {
		counter++
		Ak, Q := JacobiRotateIteration(A, maxi, maxj)
		max, maxi, maxj = Ak.findMax()
		A = Ak
		V = V.MultMatrix(Q)
		if max < eps {
			break 
		}
	}
	lamda := NewVector(len(A.data))
	for i := range lamda.data {
		lamda.data[i] = A.data[i][i]
	}
	return lamda, V, counter
}

func (v *Vector) ScalarProd(other *Vector) float64 {
	res := 0.0
	for i := range v.data {
		res += v.data[i] * other.data[i]
	}
	return res
}

func PowerMethod(A *SquareMatrix, eps float64) (float64, *Vector, int) {
	x0 := NewVector(len(A.data))
	x0.data[0] = 1
	counter := 0
	for{
		yk := A.OperateVector(x0)
		lamda := yk.ScalarProd(x0) / x0.ScalarProd(x0)
		if A.OperateVector(yk).Add(yk.MultAlpha(-1.0 * lamda)).CubicNorm() < eps {
			return lamda, yk.MultAlpha(1.0/yk.CubicNorm()), counter
		}
		counter++
		x0 = yk.MultAlpha(1.0/yk.CubicNorm())
	}
}

func main() {
	// // Создаем матрицу A размером 15x15
	// A_data := make([][]float64, 15)
	// for i := range A_data {
	// 	A_data[i] = make([]float64, 15)
	// 	for j := 0; j <= i; j++ {
	// 		if i == j {
	// 			// Диагональные элементы
	// 			A_data[i][j] = float64(rand.Intn(20))
	// 		} else {
	// 			// Недиагональные элементы
	// 			A_data[i][j] = float64(rand.Intn(20))
	// 			A_data[j][i] = A_data[i][j]
	// 		}
	// 	}
	// }

	A := NewSquareMatrixFromData([][]float64{
		{7, 2, 15, 14, 1, 17, 1, 1, 7, 16, 14, 5, 16, 18, 17}, 
		{2, 3, 13, 1, 5, 2, 16, 0, 8, 7, 14, 3, 18, 7, 15}, 
		{15, 13, 9, 8, 3, 13, 0, 7, 7, 6, 3, 19, 12, 4, 4}, 
		{14, 1, 8, 10, 4, 12, 9, 10, 3, 1, 15, 16, 19, 1, 11}, 
		{1, 5, 3, 4, 4, 6, 15, 2, 17, 19, 18, 14, 6, 18, 6}, 
		{17, 2, 13, 12, 6, 17, 5, 6, 16, 15, 7, 1, 0, 1, 17}, 
		{1, 16, 0, 9, 15, 5, 12, 11, 2, 0, 19, 16, 1, 13, 8}, 
		{1, 0, 7, 10, 2, 6, 11, 12, 19, 7, 1, 6, 5, 16, 4}, 
		{7, 8, 7, 3, 17, 16, 2, 19, 12, 3, 1, 6, 3, 6, 14}, 
		{16, 7, 6, 1, 19, 15, 0, 7, 3, 4, 18, 14, 13, 18, 18}, 
		{14, 14, 3, 15, 18, 7, 19, 1, 1, 18, 4, 15, 13, 6, 9}, 
		{5, 3, 19, 16, 14, 1, 16, 6, 6, 14, 15, 12, 14, 1, 14}, 
		{16, 18, 12, 19, 6, 0, 1, 5, 3, 13, 13, 14, 4, 2, 0}, 
		{18, 7, 4, 1, 18, 1, 13, 16, 6, 18, 6, 1, 2, 19, 13}, 
		{17, 15, 4, 11, 6, 17, 8, 4, 14, 18, 9, 14, 0, 13, 2},
	})

	fmt.Printf("\n==========================  ВХОДНЫЕ ДАННЫЕ  ==========================\n\n")
	fmt.Printf("A:\n%s\n", A)
	
	fmt.Printf("\n\n==========================  МЕТОД ОТРАЖЕНИЙ ЯКОБИ  ==========================\n\n")

	eps := 1e-4
	lamda, V, count := JacobiRotate(NewSquareMatrixFromData(A.data), eps)
	fmt.Printf("Итераций выполнено: %d\n", count)
	for i, l := range lamda.data {
		v := V.GetCol(i)
		errRate := A.OperateVector(v).Add(v.MultAlpha(-1 * l)).CubicNorm()
		fmt.Printf("\n%2d. Собственное значение: %10.5f\n    Собственный вектор: %s\n    Норма невязки: %e\n", i+1, l, v, errRate)
	}

	fmt.Printf("\n\n==========================  СТЕПЕННОЙ МЕТОД  ==========================\n\n")

	maxLamda, maxV, count := PowerMethod(A, eps)
	errRate := A.OperateVector(maxV).Add(maxV.MultAlpha(-1 * maxLamda)).CubicNorm()
	fmt.Printf("Итераций выполнено: %d\nСобственное значение: %10.5f\nСобственный вектор: %s\nНорма невязки: %e\n", count, maxLamda, maxV, errRate)
}