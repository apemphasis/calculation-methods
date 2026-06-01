package main

import (
	"fmt"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
	"math"
)

type Point struct {
	X, Y float64
}

func BezierCurve5(points []Point, numPoints int) []Point {
	if len(points) != 6 {
		panic("Для кривой Безье 5-й степени необходимо 6 контрольных точек")
	}

	curve := make([]Point, numPoints)
	
	for i := 0; i < numPoints; i++ {
		t := float64(i) / float64(numPoints-1)
		
		var x, y float64
		
		for k := 0; k <= 5; k++ {
			bernstein := bernsteinPoly(5, k, t)
			x += bernstein * points[k].X
			y += bernstein * points[k].Y
		}
		
		curve[i] = Point{x, y}
	}
	
	return curve
}

func bernsteinPoly(n, k int, t float64) float64 {
	if k < 0 || k > n {
		return 0
	}
	
	// Рекуррентная формула:
	// B(n,k,t) = (1-t)*B(n-1,k,t) + t*B(n-1,k-1,t)
	if n == 0 {
		return 1
	}
	if k == 0 {
		return math.Pow(1-t, float64(n))
	}
	if k == n {
		return math.Pow(t, float64(n))
	}
	
	return (1-t)*bernsteinPoly(n-1, k, t) + t*bernsteinPoly(n-1, k-1, t)
}

// plotBezierCurve строит график кривой Безье с контрольными точками
func plotBezierCurve(points []Point, title string, filename string) error {
	p := plot.New()
	p.Title.Text = title
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"
	
	numPoints := 100
	curve := BezierCurve5(points, numPoints)
	
	// Преобразуем точки кривой в формат для построения
	curvePoints := make(plotter.XYs, len(curve))
	for i, pt := range curve {
		curvePoints[i].X = pt.X
		curvePoints[i].Y = pt.Y
	}
	
	// Преобразуем контрольные точки в формат для построения
	controlPoints := make(plotter.XYs, len(points))
	for i, pt := range points {
		controlPoints[i].X = pt.X
		controlPoints[i].Y = pt.Y
	}
	
	// Создаем линию кривой
	line, err := plotter.NewLine(curvePoints)
	if err != nil {
		return err
	}
	line.Color = plotutil.Color(0)
	line.Width = vg.Points(2)
	
	// Создаем точки для контрольных точек
	scatter, err := plotter.NewScatter(controlPoints)
	if err != nil {
		return err
	}
	scatter.Color = plotutil.Color(1)
	scatter.Shape = plotutil.Shape(1) 
	scatter.Radius = vg.Points(5)
	
	p.Add(line, scatter)
	
	p.Legend.Add("Кривая Безье 5-й степени", line)
	p.Legend.Add("Контрольные точки", scatter)
	
	// Сохраняем график
	if err := p.Save(8*vg.Inch, 6*vg.Inch, filename); err != nil {
		return err
	}
	
	fmt.Printf("График сохранен в файл: %s\n", filename)
	return nil
}

func main() {
	// Первая кривая
	points1 := []Point{
		{4, 3},    // первая точка
		{3, 5},    // промежуточная
		{1, 6},    // промежуточная
		{-1, 4},   // промежуточная
		{-2, 1},   // промежуточная
		{-2, -1},  // последняя точка
	}
	
	// Вторая кривая
	points2 := []Point{
		{4, 3},     // первая точка
		{2, 4},     // промежуточная
		{0, 5},     // промежуточная
		{-2, 2},    // промежуточная
		{-1, 0},    // промежуточная
		{-2, -1},   // последняя точка
	}
	
	// Третья кривая
	points3 := []Point{
		{4, 3},     // первая точка
		{5, 1},     // промежуточная
		{3, -2},    // промежуточная
		{0, -1},    // промежуточная
		{-1, -2},   // промежуточная
		{-2, -1},   // последняя точка
	}
	
	// Четвертая кривая
	points4 := []Point{
		{4, 3},     // первая точка
		{2, 1},     // промежуточная
		{0, 2},     // промежуточная
		{-1, 0},    // промежуточная
		{-1.5, -1}, // промежуточная
		{-2, -1},   // последняя точка
	}
	
	// Построение графиков
	fmt.Println("Начинаем построение кривых Безье 5-й степени...")
	
	if err := plotBezierCurve(points1, "Кривая Безье 5-й степени - Набор 1", "bezier_curve1.png"); err != nil {
		fmt.Printf("Ошибка при построении графика 1: %v\n", err)
	}
	
	if err := plotBezierCurve(points2, "Кривая Безье 5-й степени - Набор 2", "bezier_curve2.png"); err != nil {
		fmt.Printf("Ошибка при построении графика 2: %v\n", err)
	}
	
	if err := plotBezierCurve(points3, "Кривая Безье 5-й степени - Набор 3", "bezier_curve3.png"); err != nil {
		fmt.Printf("Ошибка при построении графика 3: %v\n", err)
	}
	
	if err := plotBezierCurve(points4, "Кривая Безье 5-й степени - Набор 4", "bezier_curve4.png"); err != nil {
		fmt.Printf("Ошибка при построении графика 4: %v\n", err)
	}
	
	fmt.Println("Построение завершено!")
}