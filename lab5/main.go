package main

import (
	"fmt"
	"image/color"
	"math"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

var (
	A float64 = -2
	B float64 = 2
)

func test(x float64) float64 {
	return math.Pow(2, x) + math.Pow(2, -x)
}

func f1(x float64) float64 {
	return x * math.Sinh(x)
}

func f2(x float64) float64 {
	return math.Abs(5*x + 3)
}

func GetNewtonCoef(x, f []float64) ([]float64, error) {
	if len(x) != len(f) {
		return nil, fmt.Errorf("arr length must be equal")
	}
	N := len(x) + 1

	matrix := make([][]float64, N)
	for i := 0; i < len(matrix); i++ {
		matrix[i] = make([]float64, N-1)
	}

	matrix[0] = x
	matrix[1] = f

	for i := 2; i < len(matrix); i++ {
		gap := i - 1
		for j := 0; j < len(matrix[0])-i+1; j++ {
			el := (matrix[i-1][j+1] - matrix[i-1][j]) / (matrix[0][j+gap] - matrix[0][j])
			matrix[i][j] = el
		}
	}

	result := make([]float64, N-1)
	for i := 0; i < N-1; i++ {
		result[i] = matrix[i+1][0]
	}

	// fmt.Println("coef matrix")

	// for i := range matrix[0] {
	// 	for j := range matrix {
	// 		fmt.Printf(" %f ", matrix[j][i])
	// 	}
	// 	fmt.Println()
	// }

	// fmt.Println()

	// for i := range matrix {
	// 	for j := range matrix[0] {
	// 		fmt.Printf(" %f ", matrix[i][j])
	// 	}
	// 	fmt.Println()
	// }

	// fmt.Println()

	return result, nil
}

func GetUniformVertexes(n int) []float64 {
	length := B - A
	res := make([]float64, n+1)

	for i := range res {
		res[i] = A + length/float64(n)*float64(i)
	}

	return res
}

func GetChebyshovskyVertexes(n int) []float64 {
	res := make([]float64, n+1)
	div := 2 * float64(n+1)
	for i := range res {
		res[i] = (A+B)/2 + (B-A)/2*math.Cos(math.Pi*(2*float64(i)+1)/div)
	}
	return res
}

func BuildNewton(n int, function func(float64) float64, getVertex func(int) []float64) func(float64) float64 {
	x := getVertex(n)
	fmt.Printf("x: %#v\n", x)

	f := make([]float64, n+1)
	for i := range f {
		f[i] = function(x[i])
	}

	fmt.Printf("f: %#v\n", f)

	coefs, _ := GetNewtonCoef(x, f)
	fmt.Printf("coefs: %#v\n", coefs)

	return func(v float64) float64 {
		res := 0.0
		//fmt.Print("polynom: ")
		for i := range coefs {
			el := coefs[i]
			//fmt.Print(coefs[i])
			for j := 0; j < i; j++ {
				el *= v - x[j]
				//fmt.Printf("*(x-%f)", x[j])
			}
			res += el
			//fmt.Print(" + ")
		}
		//fmt.Println()
		return res
	}
}

func BuildSpline(n int, function func(float64) float64) ([]func(float64) float64, []float64) {
	x := GetUniformVertexes(n)

	f := make([]float64, n+1)
	for i := range f {
		f[i] = function(x[i])
	}

	h := make([]float64, n+1)
	for i := 1; i < n+1; i++ {
		h[i] = x[i] - x[i-1]
	}

	MMatrix := make([][]float64, n+1)
	for i := range MMatrix {
		MMatrix[i] = make([]float64, n+1)
	}
	for j := 1; j < n; j++ {
		MMatrix[j][j-1] = h[j] / (h[j] + h[j+1])
		MMatrix[j][j] = 2
		MMatrix[j][j+1] = h[j+1] / (h[j] + h[j+1])
	}
	MMatrix[0][0] = 1
	MMatrix[n][n] = 1

	vector := make([]float64, n+1)
	for i := range vector {
		if i == 0 {
			vector[i] = 2*math.Cosh(A) + A*math.Sinh(A)
		} else if i == n {
			vector[i] = 2*math.Cosh(B) + B*math.Sinh(B)
		} else {
			xes := []float64{x[i-1], x[i], x[i+1]}
			fes := make([]float64, 3)
			for j := range fes {
				fes[j] = function(xes[j])
			}
			coefs, _ := GetNewtonCoef(xes, fes)
			vector[i] = coefs[2] * 6
		}
	}

	// for i := range MMatrix {
	// 	for j := range MMatrix[0] {
	// 		fmt.Printf("%10f", MMatrix[i][j])
	// 	}
	// 	fmt.Println()
	// }

	// fmt.Println()

	// for i := range MMatrix {
	// 	fmt.Printf("%10f", vector[i])
	// }
	// fmt.Println()

	mat := NewSquareMatrixFromData(MMatrix)
	M := mat.SolveByQR(NewVectorFromData(vector))

	res := make([]func(float64) float64, n)

	for i := range res {
		res[i] = func(v float64) float64 {
			s := M.data[i] * math.Pow((x[i+1]-v), 3) / (6 * h[i+1])
			s += M.data[i+1] * math.Pow((v-x[i]), 3) / (6 * h[i+1])
			s += (f[i+1] - M.data[i+1]*h[i+1]*h[i+1]/6) * (v - x[i]) / h[i+1]
			s += (f[i] - M.data[i]*h[i+1]*h[i+1]/6) * (x[i+1] - v) / h[i+1]
			return s
		}
	}

	return res, x
}

type GraphFunc struct {
	f        func(float64) float64
	name     string
	a        float64
	b        float64
	color    color.RGBA
	isLegend bool
}

func BuildGraph(title string, funcs ...GraphFunc) {

	p := plot.New()

	p.Title.Text = title
	p.X.Label.Text = "x"
	p.Y.Label.Text = "y"

	points := make([]plotter.XYs, len(funcs))

	for i := range points {
		points[i] = make(plotter.XYs, 100)
		for j := range points[i] {
			x := funcs[i].a + (funcs[i].b-funcs[i].a)/99*float64(j)
			points[i][j].X = x
			points[i][j].Y = funcs[i].f(x)
		}

		line, _ := plotter.NewLine(points[i])
		line.LineStyle.Width = vg.Points(2)
		line.LineStyle.Color = funcs[i].color

		p.Add(line)
		if funcs[i].isLegend {
			p.Legend.Add(funcs[i].name, line)
		}
	}

	if !funcs[0].isLegend {
		line, _ := plotter.NewLine(points[0])
		line.LineStyle.Width = vg.Points(2)
		line.LineStyle.Color = funcs[len(funcs)-2].color
		p.Legend.Add(funcs[0].name, line)
	}
	p.Legend.Top = false

	p.Save(6*vg.Inch, 4*vg.Inch, title+".png")
}

func main() {
	// x := []float64{3, 0, 1, 4, 5}
	// f := []float64{3, 5, 3, 2, 1}

	// res, _ := GetNewtonCoef(x, f);

	// for _, el := range res {
	// 	fmt.Printf(" %f ", el)
	// }

	for i := 2; i <= 16; i *= 2 {
		newtonf1 := BuildNewton(i, f1, GetUniformVertexes)
		newtonf2 := BuildNewton(i, f2, GetUniformVertexes)
		BuildGraph(fmt.Sprintf("Задание 1/Интерполяция многочленом Ньтона ф-и f1, n = %d", i),
			GraphFunc{newtonf1, "Многочлен Ньютона", A, B, color.RGBA{255, 0, 0, 127}, true},
			GraphFunc{f1, "f1", A, B, color.RGBA{0, 0, 255, 127}, true},
		)
		BuildGraph(fmt.Sprintf("Задание 1/Интерполяция многочленом Ньтона ф-и f2, n = %d", i),
			GraphFunc{newtonf2, "Многочлен Ньютона", A, B, color.RGBA{255, 0, 0, 127}, true},
			GraphFunc{f2, "f2", A, B, color.RGBA{0, 0, 255, 127}, true},
		)
		newtonf1_ch := BuildNewton(i, f1, GetChebyshovskyVertexes)
		newtonf2_ch := BuildNewton(i, f2, GetChebyshovskyVertexes)
		BuildGraph(fmt.Sprintf("Задание 2/Интерполяция многочленом Ньтона ф-и f1, n = %d", i),
			GraphFunc{newtonf1_ch, "Многочлен Ньютона", A, B, color.RGBA{255, 0, 0, 127}, true},
			GraphFunc{f1, "f1", A, B, color.RGBA{0, 0, 255, 127}, true},
		)
		BuildGraph(fmt.Sprintf("Задание 2/Интерполяция многочленом Ньтона ф-и f2, n = %d", i),
			GraphFunc{newtonf2_ch, "Многочлен Ньютона", A, B, color.RGBA{255, 0, 0, 127}, true},
			GraphFunc{f2, "f2", A, B, color.RGBA{0, 0, 255, 127}, true},
		)
		spline_f1, vertexes_f1 := BuildSpline(i, f1)
		splineGraph_f1 := make([]GraphFunc, len(spline_f1))
		for i := range splineGraph_f1 {
			var red int = 255 / len(spline_f1) * i
			splineGraph_f1[i] = GraphFunc{
				spline_f1[i],
				"spline",
				vertexes_f1[i],
				vertexes_f1[i+1],
				color.RGBA{uint8(red), 0, 0, 255},
				false,
			}
		}
		splineGraph_f1 = append(splineGraph_f1, GraphFunc{
			f1,
			"f1",
			A,
			B,
			color.RGBA{0, 0, 255, 127},
			true,
		})
		BuildGraph(fmt.Sprintf("Задание 3/Интерполяция сплайном ф-и f1, n = %d", i), splineGraph_f1...)

		spline_f2, vertexes_f2 := BuildSpline(i, f2)
		splineGraph_f2 := make([]GraphFunc, len(spline_f2))
		for i := range splineGraph_f2 {
			var red int = 255 / len(spline_f2) * i
			splineGraph_f2[i] = GraphFunc{
				spline_f2[i],
				"spline",
				vertexes_f2[i],
				vertexes_f2[i+1],
				color.RGBA{uint8(red), 0, 0, 255},
				false,
			}
		}
		splineGraph_f2 = append(splineGraph_f2, GraphFunc{
			f2,
			"f2",
			A,
			B,
			color.RGBA{0, 0, 255, 127},
			true,
		})
		BuildGraph(fmt.Sprintf("Задание 3/Интерполяция сплайном ф-и f2, n = %d", i), splineGraph_f2...)
	}

}
