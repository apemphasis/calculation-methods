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

func MultArrQ(arr []*SquareMatrix) *SquareMatrix {
	res := NewSquareMatrixFromData(arr[len(arr) - 1].data)
	for i := len(arr) - 2; i >= 0; i--{
		res = arr[i].MultMatrix(res)
	}
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
			res += fmt.Sprintf("%12.3e", v.data[i])
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
			res += fmt.Sprintf("%12.3f", v.data[i])
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
				res += fmt.Sprintf("%12.3e", m.data[i][j])
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

func main() {

	matrix := SquareMatrix{data: make([][]float64, 15)}
	for i := range matrix.data {
		matrix.data[i] = make([]float64, 15)
		for j := range matrix.data[i] {
			if i == j {
				matrix.data[i][j] = 5 * math.Sqrt(float64(i))
				continue
			}else{
				matrix.data[i][j] = math.Sqrt(float64(i)) + math.Sqrt(float64(j))
			}
		}
	}

	b_data := make([]float64, 15)

	for i := range(b_data){
		b_data[i] = 4.5 * math.Sqrt(float64(i))
	}

	//A := NewSquareMatrixFromData([][]float64{{-2, -3, 3}, {4, 3, -3}, {4, 0, 9}})
	//b := NewVectorFromData([]float64{1, 1, 4})
	A := &matrix
	b := NewVectorFromData(b_data)

	fmt.Println(b.ToStringF())

	// метод отражений
	deg := len(A.data)
	R := NewSquareMatrixFromData(A.data)
	QArr := make([]*SquareMatrix, 0, deg - 1)
	for i := range(deg - 1){
		H := R.ExtractMatrix(i, deg - 1, i, deg - 1).CountQ(0)
		Qi := H.CompleteMatrix(deg)
		Qi.ClearZeros()
		
		R = Qi.MultMatrix(R)
		R.ClearZeros()
		QArr = append(QArr, Qi)
		fmt.Println(Qi.ToString())
	}

	fmt.Printf("\n==============================================\n\n")

	fmt.Printf("QR - РАЗЛОЖЕНИЕ ДЛЯ A\n\n")

	fmt.Println("Q:")
	Q := MultArrQ(QArr)
	fmt.Println(Q.ToString())

	fmt.Println("R:")
	fmt.Println(R.ToString())

	fmt.Printf("\n==============================================\n\n")

	fmt.Println("РЕШЕНИЕ:")

	Qtb := Q.Transpose().OperateVector(b)
	solution := GaussReverse(R, Qtb)
	fmt.Println(solution.ToStringF())
	fmt.Println(b.ToStringF())
	fmt.Println(A.OperateVector(solution).ToStringF())

	fmt.Printf("\n==============================================\n\n")

	fmt.Println("ПРОВЕРКА:")

	fmt.Println("Q*R:")
	QR := Q.MultMatrix(R)
	QR.ClearZeros()
	fmt.Println(QR.ToString())

	fmt.Println("A:")
	fmt.Println(A.ToString())
}
