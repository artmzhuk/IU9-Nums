package gauss

import "math"

func Evaluate(xVals, a, b, c, d []float64, xArg float64) float64 {
	n := len(xVals)
	// Найти индекс i такого, что x_i <= xArg < x_{i+1}
	i := 0
	for j := 0; j < n-1; j++ {
		if xVals[j] <= xArg && xArg < xVals[j+1] {
			i = j
			break
		}
	}
	if xArg == xVals[n-1] {
		i = n - 2 // Берем последний интервал
	}

	dx := xArg - xVals[i]
	return a[i] + b[i]*dx + c[i]*dx*dx + d[i]*dx*dx*dx
}

func CubicSplineCoefficients(x, y []float64) ([]float64, []float64, []float64, []float64) {
	n := len(x)
	h := (x[n-1] - x[0]) / float64(n-1)

	c := SolveTridiagonalSystem(n, h, y)
	a := make([]float64, n-1)
	copy(a, y[:n-1])

	b := make([]float64, n-1)
	for i := 0; i < n-2; i++ {
		b[i] = (y[i+1]-y[i])/h - h/3.0*(c[i+1]+2*c[i])
	}
	b[n-2] = (y[n-1]-y[n-2])/h - 2.0/3.0*h*c[n-2]

	d := make([]float64, n-1)
	for i := 0; i < n-2; i++ {
		d[i] = (c[i+1] - c[i]) / (3 * h)
	}
	d[n-2] = -c[n-2] / (3 * h)

	return a, b, c[:n-1], d
}

func SolveTridiagonalSystem(n int, h float64, y []float64) []float64 {
	a := make([][]float64, n)
	for i := range a {
		a[i] = make([]float64, n+1)
	}

	for i := 1; i < n-1; i++ {
		a[i][i-1] = 1
		a[i][i] = 4
		a[i][i+1] = 1
		a[i][n] = (y[i+1] - 2*y[i] + y[i-1]) / (h * h)
	}

	a[0][0] = 1

	a[n-1][n-1] = 1

	return GaussElimination(a)
}

func GaussElimination(a [][]float64) []float64 {
	n := len(a)
	x := make([]float64, n)

	for i := 0; i < n; i++ {
		maxEl := math.Abs(a[i][i])
		maxRow := i
		for k := i + 1; k < n; k++ {
			if math.Abs(a[k][i]) > maxEl {
				maxEl = math.Abs(a[k][i])
				maxRow = k
			}
		}

		for k := i; k <= n; k++ {
			a[maxRow][k], a[i][k] = a[i][k], a[maxRow][k]
		}

		for k := i + 1; k < n; k++ {
			c := -a[k][i] / a[i][i]
			for j := i; j <= n; j++ {
				if i == j {
					a[k][j] = 0
				} else {
					a[k][j] += c * a[i][j]
				}
			}
		}
	}

	x[n-1] = a[n-1][n] / a[n-1][n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = a[i][n]
		for j := i + 1; j < n; j++ {
			x[i] -= a[i][j] * x[j]
		}
		x[i] /= a[i][i]
	}

	return x
}
