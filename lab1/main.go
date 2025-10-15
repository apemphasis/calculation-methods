package main //Kozlovsky_13gr_5v

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

func (v *Vector) ClearZeros() {
	for i := range(v.data){
		if math.Abs(v.data[i]) < 1e-10{
			v.data[i] = 0
		}
	}
}

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

func (m *SquareMatrix) CountQ(col int) *SquareMatrix {
	W := CountW(m.GetCol(col))
	return CountH(W)
}

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

func (v *Vector) ToStringF() string {
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
	res.ClearZeros()
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

func (m *SquareMatrix) ClearZeros() {
	for i := range(m.data){
		for j := range(m.data[i]){
			if math.Abs(m.data[i][j]) < 1e-10{
				m.data[i][j] = 0
			}
		}
	}
} //! CLEAR ZEROS

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

func GaussReverse(m *SquareMatrix, b *Vector) *Vector {
	n := len(b.data)

	solution := NewVector(n)
	solution.data[n - 1] = b.data[n-1] / m.data[n-1][n-1]

	for i := n - 2; i >= 0; i--{
		sum := 0.0
		for j := 1; i + j < n; j++ {
			sum -= m.data[i][i+j] * solution.data[i + j]
		}
		solution.data[i] = (sum + b.data[i]) / m.data[i][i]
	}

	return solution
}

func Gauss(m *SquareMatrix, b *Vector) *Vector {
	n := len(b.data)

	solution := NewVector(n)
	solution.data[0] = b.data[0] / m.data[0][0]

	for i := 1; i < n; i++{
		sum := 0.0
		for j := 0; j < i; j++ {
			sum -= m.data[i][j] * solution.data[j]
		}
		solution.data[i] = (sum + b.data[i]) / m.data[i][i]
	}

	return solution
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
		H := R.ExtractMatrix(i, deg - 1, i, deg - 1).CountQ(0)
		Qi := H.CompleteMatrix(deg)
		Qi.ClearZeros()
		
		R = Qi.MultMatrix(R)
		R.ClearZeros()
		QArr = append(QArr, Qi)
	}

	Q := NewSquareMatrixFromData(QArr[len(QArr) - 1].data)
	for i := len(QArr) - 2; i >= 0; i--{
		Q = QArr[i].MultMatrix(Q)
	}

	return Q, R
}

func (m *SquareMatrix) SolveByQR(b *Vector) *Vector{
	Q, R := m.FindQR()

	Qtb := Q.Transpose().OperateVector(b)
	solution := GaussReverse(R, Qtb)

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

	fmt.Printf("\n==========================  INPUT DATA  ==========================\n\n")

	A := NewSquareMatrixFromData(A_data)
	b := NewVectorFromData(b_data)

	fmt.Println("A:")
	fmt.Println(A.ToString())

	fmt.Println("b:")
	fmt.Println(b.ToString())

	fmt.Printf("\n\n==========================  LDL^t DECOMPOSITION  ==========================\n\n")

	L, D, Lt := A.FindLDLt()

	fmt.Println("L:")
	fmt.Println(L.ToString())

	fmt.Println("D:")
	fmt.Println(D.ToString())

	fmt.Println("L^t:")
	fmt.Println(Lt.ToString())

	fmt.Printf("\n\n==========================  SOLVE BY LDL^t  ==========================\n\n")

	x := A.SolveByLDLt(b)

	fmt.Println("SOLUTION:")
	fmt.Println(x.ToStringF())

	fmt.Printf("\n\n==========================  INVERSE BY LDL^t  ==========================\n\n")

	inverseA := A.FindInverse()

	fmt.Println("INVERSED:")
	fmt.Println(inverseA.ToString())

	fmt.Printf("\n\n==========================  QR DECOMPOSITION  ==========================\n\n")

	Q, R := A.FindQR()

	fmt.Println("Q:")
	fmt.Println(Q.ToString())

	fmt.Println("R:")
	fmt.Println(R.ToString())

	fmt.Printf("\n\n==========================  SOLVE BY QR  ==========================\n\n")

	solution := A.SolveByQR(b)

	fmt.Println("SOLUTION:")
	fmt.Println(solution.ToStringF())
	
	fmt.Printf("\n\n==========================  TEST  ==========================\n\n")

	fmt.Println("L * D * L^t:")
	fmt.Println(L.MultMatrix(D.MultMatrix(Lt)).ToString())

	fmt.Println("Q * R:")
	QR := Q.MultMatrix(R)
	fmt.Println(QR.ToString())

	// fmt.Println("A:")
	// fmt.Println(A.ToString())
}
