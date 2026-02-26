package main

import (
	"fmt"
	"math"
)

func getBMatrix(A *SquareMatrix) *SquareMatrix {
	B := NewSquareMatrix(len(A.data))
	for i := range A.data {
		for j := range A.data[i] {
			if i == j {
				B.data[i][j] = 0
			} else {
				B.data[i][j] = A.data[i][j] / A.data[i][i]
			}
		}
	}
	return B
}

func MPIIteration(A *SquareMatrix, x0, b *Vector) *Vector {
	xk := NewVector(len(x0.data))

	for i := 0; i < len(x0.data); i++ {
		sum := b.data[i]
		for j := 0; j < len(x0.data); j++ {
			if i != j {
				sum -= A.data[i][j] * x0.data[j]
			}
		}
		xk.data[i] = 1.0 / A.data[i][i] * sum
	}
	return xk
}

func MPI(A *SquareMatrix, x0, b *Vector, eps float64) (*Vector, int, float64) {
	B := getBMatrix(A)
	apost_B := B.MaxRowNorm() / (1.0 - B.MaxRowNorm())
	counter := 0
	missRate := 0.0

	for  {
		counter++
		xk := MPIIteration(A, x0, b)
		missRate = apost_B * xk.Add(x0.MultAlpha(-1)).CubicNorm()
		x0 = NewVectorFromData(xk.data)
		if  missRate < eps{
			break
		}
	}
	return x0, counter, missRate
}

func ZeidelIteration(A *SquareMatrix, x0, b *Vector) *Vector {
	xk := NewVector(len(x0.data))

	for i := 0; i < len(x0.data); i++ {
		sum := b.data[i]
		for j := 0; j < len(x0.data); j++ {
			if j == i { continue }

			if j < i {
				sum -= A.data[i][j] * xk.data[j]
			} else {
				sum -= A.data[i][j] * x0.data[j]
			}
		}
		xk.data[i] = 1.0 / A.data[i][i] * sum
	}
	return xk
}

func Zeidel(A *SquareMatrix, x0, b *Vector, eps float64, limit int) (*Vector, int, float64) {
	B := getBMatrix(A)
	apost_B := B.MaxRowNorm() / (1.0 - B.MaxRowNorm())
	counter := 0
	missRate := 0.0

	for {
		counter++
		xk := ZeidelIteration(A, x0, b)
		missRate = apost_B * xk.Add(x0.MultAlpha(-1)).CubicNorm()
		x0 = NewVectorFromData(xk.data)
		if limit == counter{
			break
		}
	}
	return x0, counter, missRate
}

func RelaxIteration(A *SquareMatrix, x0, b *Vector, w float64) *Vector {
	xk := NewVector(len(x0.data))
	for i := 0; i < len(x0.data); i++ {
		sum := b.data[i]
		for j := 0; j < len(x0.data); j++ {
			if j == i { continue }

			if j < i {
				sum -= A.data[i][j] * xk.data[j]
			} else {
				sum -= A.data[i][j] * x0.data[j]
			}
		}
		xk.data[i] = (1.0 - w) * x0.data[i] + w / A.data[i][i] * sum
	}
	return xk
}

func Relax(A *SquareMatrix, x0, b *Vector, eps, w float64, limit int) (*Vector, int, float64) {
	B := getBMatrix(A)
	apost_B := B.MaxRowNorm() / (1.0 - B.MaxRowNorm())
	counter := 0
	missRate := 0.0

	for {
		counter++
		xk := RelaxIteration(A, x0, b, w)
		missRate = apost_B * xk.Add(x0.MultAlpha(-1)).CubicNorm()
		x0 = NewVectorFromData(xk.data)
		if missRate < eps{
			break
		}
	}
	return x0, counter, missRate
}

func main() {
	// Создаем матрицу A размером 15x15
	A_data := make([][]float64, 15)
	for i := range A_data {
		A_data[i] = make([]float64, 15)
		for j := range A_data[i] {
			if i == j {
				// Диагональные элементы
				A_data[i][j] = float64(8 * (i + 1))
				continue
			} else {
				// Недиагональные элементы
				A_data[i][j] = -0.01 / float64((i + 1) * (i + 1) * (i + 1) + (j + 1))
			}
		}
	}

	// Создаем вектор правой части b
	tx_data := make([]float64, 15)
	
	for i := range tx_data {
		tx_data[i] = float64(i + 4)
	}

	fmt.Printf("\n==========================  ВХОДНЫЕ ДАННЫЕ  ==========================\n\n")

	A := NewSquareMatrixFromData(A_data)
	tx := NewVectorFromData(tx_data)

	b := A.OperateVector(tx)

	fmt.Println("Искомое решение:")
	fmt.Println(tx)
	fmt.Println("\nA:")
	fmt.Println(A)

	fmt.Println("\nb:")
	fmt.Println(b)

	
	eps := 1e-5
	fmt.Printf("\nПогрешность: %e\n", eps)

	B := getBMatrix(A)
	
	expected_k := int(math.Ceil(math.Log(eps * (1.0 - B.MaxRowNorm()) / b.CubicNorm()) / math.Log(B.MaxRowNorm()))) - 1
	fmt.Printf("\nОжидаемое кол-во итераций: %d\n", expected_k)

	fmt.Printf("\n===== Метод МПИ: =====\n")

	mpiRes, mpiCount, mpiErr := MPI(A, NewVector(15), b, eps)

	fmt.Println("Решение:")
	fmt.Println(mpiRes)
	fmt.Printf("\nПонадобилось итераций: %d\nАпостериорная погрешность: %e\n\n", mpiCount, mpiErr)
	
	fmt.Printf("\n===== Метод Зейделя: =====\n")
	
	ziedelRes, ziedelCount, ziedelErr := Zeidel(A, NewVector(15), b, eps, mpiCount)

	fmt.Println("Решение:")
	fmt.Println(ziedelRes)
	fmt.Printf("\nПонадобилось итераций: %d\nАпостериорная погрешность: %e\n\n", ziedelCount, ziedelErr)

	fmt.Printf("\n===== Метод Релаксации W = 0.5: =====\n")
	
	relaxRes, relaxCount, relaxErr := Relax(A, NewVector(15), b, eps, 0.5, mpiCount)

	fmt.Println("Решение:")
	fmt.Println(relaxRes)
	fmt.Printf("\nПонадобилось итераций: %d\nАпостериорная погрешность: %e\n\n", relaxCount, relaxErr)

	fmt.Printf("\n===== Метод Релаксации W = 1.5: =====\n")
	
	relaxRes, relaxCount, relaxErr = Relax(A, NewVector(15), b, eps, 1.5, mpiCount)

	fmt.Println("Решение:")
	fmt.Println(relaxRes)
	fmt.Printf("\nПонадобилось итераций: %d\nАпостериорная погрешность: %e\n\n", relaxCount, relaxErr)
}