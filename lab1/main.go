package main

import (
	"fmt"
)

type Spline struct {
	Coefficients [][]float64
}

func cubicSplineInterpolation(x, y []float64) Spline {
	n := len(x)
	h := make([]float64, n-1)
	alpha := make([]float64, n-1)
	l := make([]float64, n)
	mu := make([]float64, n)
	z := make([]float64, n)

	for i := 0; i < n-1; i++ {
		h[i] = x[i+1] - x[i]
	}

	for i := 1; i < n-1; i++ {
		alpha[i] = (3.0/h[i])*(y[i+1]-y[i]) - (3.0/h[i-1])*(y[i]-y[i-1])
	}

	l[0], mu[0], z[0] = 1.0, 0.0, 0.0
	for i := 1; i < n-1; i++ {
		l[i] = 2.0*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1]
		mu[i] = h[i] / l[i]
		z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]
	}

	l[n-1], z[n-1] = 1.0, 0.0
	c := make([]float64, n)
	b := make([]float64, n-1)
	d := make([]float64, n-1)

	cSolution := Gauss(h, l, mu, z)

	for j := n - 2; j >= 0; j-- {
		c[j] = cSolution[j]
		b[j] = (y[j+1]-y[j])/h[j] - h[j]*(c[j+1]+2.0*c[j])/3.0
		d[j] = (c[j+1] - c[j]) / (3.0 * h[j])
	}

	a := append([]float64{}, y...)
	coefficients := make([][]float64, n-1)
	for i := 0; i < n-1; i++ {
		coefficients[i] = []float64{a[i], b[i], c[i], d[i]}
	}

	return Spline{Coefficients: coefficients}
}

func Gauss(a, b, c, d []float64) []float64 {
	n := len(a)
	x := make([]float64, n)
	l := make([]float64, n)
	u := make([]float64, n)

	l[0] = b[0]
	u[0] = c[0] / l[0]
	d[0] = d[0] / l[0]

	for i := 1; i < n; i++ {
		l[i] = b[i] - a[i]*u[i-1]
		if i < n-1 {
			u[i] = c[i] / l[i]
		}
		d[i] = (d[i] - a[i]*d[i-1]) / l[i]
	}

	x[n-1] = d[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = d[i] - u[i]*x[i+1]
	}

	return x
}

func main() {
	x := []float64{1, 2, 3, 4, 5, 6}
	y := []float64{1.0002, 1.0341, 0.6, 0.40105, 0.1, 0.23975}

	spline := cubicSplineInterpolation(x, y)

	for i, coeffs := range spline.Coefficients {
		fmt.Printf(" %d: a = %f, b = %f, c = %f, d = %f\n", i+1, coeffs[0], coeffs[1], coeffs[2], coeffs[3])
	}
}
