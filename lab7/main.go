package main

import (
	"fmt"
	"math"
)

var a, b float64 = 1, 2

var eps float64 = 1e-6

func f(x float64) float64 {
	return x * math.Log(x)
}

var ACTUAl float64 = 2* math.Log(2) - 0.75
var ACTUAl_2 float64 = math.Sqrt(math.Pi) * math.Exp(0.25)

func evalSimpson(N, h float64) float64 {
	sum := 0.0

	for i := 0.0; i < N; i++ {
		temp := f(a + h*i) + 4 * f(a + h*i + h/2) + f(a + h*i + h)
		sum += temp
	}

	sum *= h/6

	// fmt.Printf("N: %.f; h: %.4f\n", N, h)
	// fmt.Println(sum)
	// fmt.Println(ACTUAl)
	// fmt.Printf("Погрешность: %e\n\n", math.Abs(sum - ACTUAl))
	return sum
}

func evalSimpsonOdd(N, h float64) float64 {
	M := N/2 
	sum := h/3 *(f(a) + f(b))

	coef := 4*h / 3
	temp := 0.0
	for i := 0.0; i < M; i++ {
		temp += f(a + h*(2*i+1))
	}
	sum += temp * coef

	coef = 2*h / 3
	temp = 0.0
	for i := 1.0; i < M; i++ {
		temp += f(a + h*(2*i))
	}	
	sum += temp * coef
	// fmt.Printf("N: %.f; M: %.f; h: %.4f\n", N, M, h)
	// fmt.Println(sum)
	// fmt.Println(ACTUAl)
	// fmt.Printf("Погрешность: %e\n\n", math.Abs(sum - ACTUAl))
	return sum
}

func evalRightRectangles(N, h float64) float64 {
	sum := 0.0
	for i := 1.0; i <= N; i++ {
		sum += f(a + h*i)
	}
	sum *= h
	// fmt.Printf("N: %.f;  h: %e\n", N, h)
	// fmt.Println(sum)
	// fmt.Println(ACTUAl)
	// fmt.Printf("Погрешность: %e\n\n", math.Abs(sum - ACTUAl))
	return sum
}

func calcHByRunge(method func(float64, float64) float64, m float64) float64 {
	N := 2.0
	
	e := 1.0
	k := 1
	for math.Abs(e) > eps && k > 1 {
		Ih := method(N, (b-a)/N)
		Ih2 := method(2*N,(b-a)/N/2)
		e = (Ih2-Ih)/(math.Pow(2, m)-1)
		fmt.Println("e: ", e)
		N *= 2
		k++
	}
	return N
}

func f2(t float64) float64 {
	return t * math.Exp(t) / 4
}

func main() {
	
	var N1, N2 float64 = 5, 6
	_ = N2

	var h1 float64 = (b-a)/N1
	
	//evalSimpson(N1, h1)


	var h2 float64 = 0.125
	// var M float64 = 4

	// evalSimpsonOdd(2*M, h2)

	// h2  = 0.1
	// M = 5

	// evalSimpsonOdd(2*M, h2)

	// evalRightRectangles(696572.0, 1.0/696572)


	N1 = calcHByRunge(evalSimpson, 4)
	h1 = (b-a)/N1

	sum := evalSimpson(N1, h1)
	fmt.Printf("N: %.f; h: %.4f\n", N1, h1)
	fmt.Println(sum)
	fmt.Println(ACTUAl)
	fmt.Printf("Погрешность: %e\n\n", math.Abs(sum - ACTUAl))

	N2 = calcHByRunge(evalSimpsonOdd, 4)
	h2 = (b-a)/N2

	sum = evalSimpsonOdd(N2, h2)
	fmt.Printf("N: %.f; h: %.4f\n", N2, h2)
	fmt.Println(sum)
	fmt.Println(ACTUAl)
	fmt.Printf("Погрешность: %e\n\n", math.Abs(sum - ACTUAl))


	N2 = calcHByRunge(evalRightRectangles, 1)
	h2 = (b-a)/N2

	sum = evalSimpsonOdd(N2, h2)
	fmt.Printf("N: %.f; h: %e\n", N2, h2)
	fmt.Println(sum)
	fmt.Println(ACTUAl)
	fmt.Printf("Погрешность: %e\n\n", math.Abs(sum - ACTUAl))

	sum = KF_NAST(32)

	fmt.Printf("\nЗначение при n: %d       = %f\n", 32, sum)
	fmt.Printf("Погрешность: %e\n\n", math.Abs(sum - ACTUAl_2))
}

func KF_NAST(n int) float64 {
	matrix := make([][]float64, n)
	for i := range matrix {
		matrix[i] = make([]float64, n)
	}

	for i := 0; i < len(matrix)-1; i++ {
		matrix[i][i+1] = math.Sqrt(float64(i+1)/2)
		matrix[i+1][i] = math.Sqrt(float64(i+1)/2)
	}

	//A := NewSquareMatrixFromData(matrix)
	//fmt.Println(A)
	lamda_eps := 1e-6
	lamda, V, _ := JacobiRotate(NewSquareMatrixFromData(matrix), lamda_eps)
	//fmt.Printf("Итераций МВЯ выполнено: %d\n", count)

	// lamda - вектор собственных значений - узлы xk
	// V - матрица собственных векторов: первая координата используется в формуле Ak
	sum := 0.0
	for i, l := range lamda.data {
		v := V.GetCol(i).Normalize()
		Ak := v.data[0] * v.data[0] * math.Sqrt(math.Pi)
		//fmt.Printf("\nX%d = %f    A%d = %f\n", i, l, i, Ak)
		sum += Ak * f2(l)
		//errRate := A.OperateVector(v).Add(v.MultAlpha(-1 * l)).CubicNorm()
		//fmt.Printf("\n%2d. Собственное значение: %10.5f\n    Собственный вектор: %s\n    Норма невязки: %e\n", i+1, l, v, errRate)
	}
	return sum
}

func KF_NAST_LOOP() (float64, int) {
	n := 8
	count := 0
	
	Sn := KF_NAST(n)
	for {
		Sn2 := KF_NAST(2*n)
		if math.Abs(Sn-Sn2) < 1e-4 * (1 + math.Abs(Sn2)) && count > 3 {
			fmt.Printf("\nРазница: %e\n", math.Abs(Sn-Sn2))
			return Sn2, 2*n
		}
		fmt.Printf("\n%d complete\n", count)
		Sn = Sn2
		n *= 2
		count++
	}
}