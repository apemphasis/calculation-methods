package main //Козловский Евгений, 13 группа, 5 вариант

import (
	"fmt"
	"math"
)

func CountW(s *Vector) *Vector {

	e_data := make([]float64, len(s.data))
	if s.data[0] != 0{

		e_data[0] = (s.data[0] / math.Abs(s.data[0])) * s.Abs()
	}else{
		e_data[0] = s.Abs()
	}
	E := NewVectorFromData(e_data)
	W := s.Add(E)
	W = W.MultAlpha(1 / W.Abs())

	return W
}

func CountH(w *Vector) *SquareMatrix {
	res := NewSquareMatrix(len(w.data))

	for i := range(res.data){
		res.data[i][i] = 1
	}

	wwt := NewSquareMatrix(len(w.data))
	for i := range(wwt.data){
		for j := range(wwt.data[i]){
			wwt.data[i][j] = w.data[i] * w.data[j]
		}
	}
	
	res = res.Add(wwt.MultAlpha(-2))
	return res
}

func (m *SquareMatrix) CountQ() *SquareMatrix {
	W := CountW(m.GetCol(0))
	return CountH(W)
}

func (m *SquareMatrix) FindLDLt() (*SquareMatrix, *SquareMatrix, *SquareMatrix) {
	A_copy := NewSquareMatrixFromData(m.data)
	D := NewSquareMatrix(len(m.data))
	L := NewSquareMatrix(len(m.data))

	for i := range(L.data) {
		L.data[i][i] = 1
	}

	for k := 0; k < len(m.data) - 1; k++{
		D.data[k][k] = A_copy.data[k][k]
		for i := k + 1; i < len(m.data) ; i++{
			L.data[i][k] = A_copy.data[i][k] / D.data[k][k]
			for j:= k + 1; j <= i; j++{
				A_copy.data[i][j] -= L.data[i][k] * A_copy.data[j][k]
			}
		}
	}

	D.data[len(m.data) - 1][len(m.data) - 1] = A_copy.data[len(m.data) - 1][len(m.data) - 1]
	Lt := L.Transpose()

	return L, D, Lt
}

func (m *SquareMatrix) SolveByLDLt(b *Vector) *Vector {
	L, D, Lt := m.FindLDLt()

	z := Gauss(L, b)
	y := GaussReverse(D, z)
	x := GaussReverse(Lt, y)

	return x
}

func (m *SquareMatrix) FindInverse() *SquareMatrix {
	res := NewSquareMatrix(len(m.data))

	for i := range(res.data) {
		e := NewVector(len(res.data))
		e.data[i] = 1
		col := m.SolveByLDLt(e)
		for j := range(col.data){
			res.data[j][i] = col.data[j]
		}
	}

	return res
}

func (m *SquareMatrix) FindQR() (*SquareMatrix, *SquareMatrix) {
	deg := len(m.data)
	R := NewSquareMatrixFromData(m.data)
	QArr := make([]*SquareMatrix, 0, deg - 1)
	for i := range(deg - 1){
		H := R.ExtractMatrix(i, deg - 1, i, deg - 1).CountQ()
		Qi := H.CompleteMatrix(deg)
		
		R = Qi.MultMatrix(R)
		QArr = append(QArr, Qi)
	}

	Q := NewSquareMatrixFromData(QArr[len(QArr) - 1].data)
	for i := len(QArr) - 2; i >= 0; i--{
		Q = QArr[i].MultMatrix(Q)
	}

	return Q, R
}

func (m *SquareMatrix) SolveByQR(b *Vector) *Vector{

	deg := len(m.data)
	R := NewSquareMatrixFromData(m.data)
	bi := NewVectorFromData(b.data)
	for i := range(deg - 1){
		H := R.ExtractMatrix(i, deg - 1, i, deg - 1).CountQ()
		Qi := H.CompleteMatrix(deg)
		
		R = Qi.MultMatrix(R)
		bi = Qi.OperateVector(bi)
	}

	solution := GaussReverse(R, bi)

	return solution
}

func main() {

	A_data := make([][]float64, 15)
	for i := range A_data {
		A_data[i] = make([]float64, 15)
		for j := range A_data[i] {
			if i == j {
				A_data[i][j] = 5 * math.Sqrt(float64(i + 1))
				continue
			}else{
				A_data[i][j] = math.Sqrt(float64(i + 1)) + math.Sqrt(float64(j + 1))
			}
		}
	}

	b_data := make([]float64, 15)

	for i := range(b_data){
		b_data[i] = 4.5 * math.Sqrt(float64(i + 1))
	}

	fmt.Printf("\n==========================  ВХОДНЫЕ ДАННЫЕ  ==========================\n\n")

	A := NewSquareMatrixFromData(A_data)
	b := NewVectorFromData(b_data)

	fmt.Println("A:")
	fmt.Println(A.ToString())

	fmt.Println("b:")
	fmt.Println(b.ToString())

	fmt.Printf("\n\n==========================  LDL^t РАЗЛОЖЕНИЕ  ==========================\n\n")
	
	L, D, Lt := A.FindLDLt()
	fmt.Println("L:")
	fmt.Println(L.ToString())

	fmt.Println("D:")
	fmt.Println(D.ToString())

	fmt.Println("L^t:")
	fmt.Println(Lt.ToString())

	fmt.Println()

	fmt.Print("Норма ||A - LDL^t|| (максимальная по строкам): ")
	fmt.Println(A.Add(L.MultMatrix(D.MultMatrix(Lt)).MultAlpha(-1)).MaxRowNorm())

	fmt.Printf("\n\n==========================  РЕШЕНИЕ С ПОМОЩЬЮ LDL^t  ==========================\n\n")

	x := A.SolveByLDLt(b)

	fmt.Println("Решение:")
	fmt.Println(x.ToString())
	fmt.Println()
	fmt.Print("Норма невязки (евклидова): ")
	fmt.Println(b.Add(A.OperateVector(x).MultAlpha(-1)).Abs())

	fmt.Printf("\n\n==========================  ОБРАТНАЯ С ПОМОЩЬЮ LDL^t  ==========================\n\n")

	inverseA := A.FindInverse()

	fmt.Println("Обратная:")
	fmt.Println(inverseA.ToString())

	fmt.Println()

	fmt.Print("Норма ||E - AA^(-1)|| (максимальная по строкам): ")

	E := NewSquareMatrix(15)
	for i := range(E.data){
		E.data[i][i] = 1
	}

	fmt.Println(E.Add(A.MultMatrix(A.FindInverse()).MultAlpha(-1)).MaxRowNorm())

	fmt.Printf("\n\n==========================  QR РАЗЛОЖЕНИЕ  ==========================\n\n")

	Q, R := A.FindQR()

	fmt.Println("Q:")
	fmt.Println(Q.ToString())

	fmt.Println("R:")
	fmt.Println(R.ToString())

	fmt.Println()

	fmt.Print("Норма ||A - QR|| (максимальная по строкам): ")
	fmt.Println(A.Add(Q.MultMatrix(R).MultAlpha(-1)).MaxRowNorm())

	fmt.Printf("\n\n==========================  РЕШЕНИЕ МЕТОДОМ ОТРАЖЕНИЙ  ==========================\n\n")

	solution := A.SolveByQR(b)

	fmt.Println("Решение:")
	fmt.Println(solution.ToString())
	fmt.Println()
	fmt.Print("Норма невязки (евклидова): ")
	fmt.Println(b.Add(A.OperateVector(solution).MultAlpha(-1)).Abs())
	
	fmt.Println()
}
