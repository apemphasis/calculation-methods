package main

import (
	"fmt"
	"math"
)

const q = 0.52

const eps = 1e-6

func phi(x float64) float64 {
    return math.Sqrt(math.Sqrt(18 * x * x - 5 * x + 8))
}

func f(x float64) float64 {
	return x*x*x*x - 18 * x*x + 5 * x - 8
}

func f_1(x float64) float64 {
	return 4*x*x*x - 36*x + 5 
}

func mpiMiss(curr, prev float64) float64 {
	return q / (1-q) * math.Abs(curr - prev)
}

func main() {
	fmt.Println("=== Метод простой итерации ===")

	x := 4.5

	for counter := 0; ; counter++{
		t := phi(x)
		if miss := mpiMiss(t, x); miss <= eps {
			fmt.Printf("x: %.10f\n", t)
			fmt.Printf("Miss: %e\n", miss)
			fmt.Printf("Done: %d iterations\n", counter)
			break
		}
		x = t
	}

	fmt.Println()
	fmt.Println("=== Метод Ньютона ===")

	x = 5.0

	for counter := 0; ; counter++{
		t := x - f(x)/f_1(x)
		if math.Abs(t-x) < eps {
			fmt.Printf("x: %.10f\n", t)
			fmt.Printf("Miss: %e\n", math.Abs(t-x))
			fmt.Printf("Done: %d iterations\n", counter)
			break
		}
		x=t
	}
}