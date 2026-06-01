package main

import (
	"fmt"
	"image/color"
	"math"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

// Исходная функция
func f(x, u float64) float64 {
	return (x*x - 1) / (u*u + 1)
}

// Функция F(z) для метода Ньютона (здесь t - это t_{i+1})
func F(z, t, h, cst float64) float64 {
	return 3*h*(t*t-1)/(8*(z*z+1)) - z + cst
}

// Производная F'(z) для метода Ньютона
func F_1(z, t, h float64) float64 {
	return -3*h*(t*t-1)*z/(4*(z*z+1)*(z*z+1)) - 1
}

const (
	a   float64 = 0
	b   float64 = 1
	eps float64 = 1e-3
)

// Норма по правилу Рунге. Учитываем, что массивы имеют размер N+1 и 2N+1
func CubicNorm(coarse, fine []float64) float64 {
	m := 0.0
	for i := 0; i < len(coarse); i++ {
		if 2*i < len(fine) {
			m = max(m, math.Abs(coarse[i]-fine[2*i]))
		}
	}
	return m
}

func SaveChart(y []float64, title string) {
	p := plot.New()
	p.Title.Text = title
	p.X.Label.Text = "x"
	p.Y.Label.Text = "y"

	points := make(plotter.XYs, len(y))
	tau := (b - a) / float64(len(y)-1) // len(y)-1 это N отрезков

	for i := range points {
		points[i].X = a + tau*float64(i)
		points[i].Y = y[i]
	}

	line, _ := plotter.NewLine(points)
	line.LineStyle.Width = vg.Points(2)
	line.LineStyle.Color = color.RGBA{255, 0, 0, 255}

	p.Add(line)
	p.Save(6*vg.Inch, 4*vg.Inch, title+".png")
}

// Вспомогательная функция: 1 полный прогон RK2 для фиксированного N
// Возвращает массив размером N+1
func rungeKuttFixedN(N int) []float64 {
	y := make([]float64, N+1)
	y[0] = 0 // Начальное условие
	h := (b - a) / float64(N)

	for i := 0; i < N; i++ {
		t := a + h*float64(i)
		
		// Явные формулы по таблице Бутчера
		k1 := f(t, y[i])
		k2 := f(t+2*h, y[i]+2*h*k1)
		
		y[i+1] = y[i] + h*(0.75*k1 + 0.25*k2)
	}
	return y
}

func RungeKutt() []float64 {
	N := 4
	y_prev := rungeKuttFixedN(N)
	
	for {
		N *= 2
		y := rungeKuttFixedN(N)
		
		if CubicNorm(y_prev, y)/3.0 <= eps { // Погрешность O(h^2)
			return y
		}
		y_prev = y
	}
}

func Adams() []float64 {
	N := 4
	y_prev := make([]float64, 0)

	for count := 0; count < 15; count++ { // Ограничитель на случай бесконечного цикла
		h := (b - a) / float64(N)
		
		// 1. ПОЛУЧАЕМ ПРАВИЛЬНЫЙ РАЗГОН ДЛЯ ТЕКУЩЕГО ШАГА H
		startPts := rungeKuttFixedN(N)
		y := make([]float64, N+1)
		y[0] = startPts[0]
		y[1] = startPts[1]
		y[2] = startPts[2]

		// 2. ОСНОВНОЙ ЦИКЛ АДАМСА
		for i := 2; i < N; i++ { // Идем до N, чтобы посчитать y[N]
			t_i := a + h*float64(i)
			t_next := a + h*float64(i+1) // <--- ИСПРАВЛЕНО: время следующего шага
			
			// Считаем константу из известных предыдущих узлов
			cst := y[i] + (h/24.0)*(19*f(t_i, y[i]) - 5*f(t_i-h, y[i-1]) + f(t_i-2*h, y[i-2]))
			
			curr := y[i] // Начальное приближение
			
			for counter := 0; counter < 1000; counter++ {
				// В F и F_1 обязательно передаем t_next, так как неявность зависит от f_{i+1}
				next := curr - F(curr, t_next, h, cst)/F_1(curr, t_next, h) 
				
				if math.Abs(next-curr) < 1e-6 { // Точность Ньютона лучше брать с запасом
					curr = next
					break
				}
				curr = next
			}
			y[i+1] = curr
		}

		// 3. ПРОВЕРКА ПО РУНГЕ
		if len(y_prev) > 0 {
			if CubicNorm(y_prev, y)/15.0 <= eps { // Погрешность O(h^4)
				return y
			}
		}
		
		y_prev = y
		N *= 2
	}

	return y_prev
}

func main() {
	fmt.Println("Метод Рунге-Кутта 2-го порядка: ")
	yRK := RungeKutt()
	fmt.Printf("Кол-во узлов: %5d\n", len(yRK)) // Выведет кол-во точек (N+1)
	SaveChart(yRK, "Runge-Kutt")

	fmt.Println()

	fmt.Println("Метод Адамса 4-го порядка: ")
	yAd := Adams() // Теперь не зависит от результата RK напрямую!
	fmt.Printf("Кол-во узлов: %5d\n", len(yAd))
	SaveChart(yAd, "Adams")

	sub := make([]float64, len(yAd))
	for i := range sub {
		sub[i] = math.Abs(yRK[2*i] - yAd[i])
	}
	SaveChart(sub, "sub")
}