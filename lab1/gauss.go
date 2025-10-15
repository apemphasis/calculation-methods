package main

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