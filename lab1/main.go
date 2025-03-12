package main

import (
	"fmt"
	"lab1/gauss"
	"math"
)

// 4arctg(x)/(x^2+1)
func f(x float64) float64 {
	return 4 * math.Atan(x) / (x*x + 1)
}

func main() {
	a, b := 0.0, 4.0
	n := 32
	h := (b - a) / float64(n-1)

	x := make([]float64, n)
	y := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = a + float64(i)*h
		y[i] = f(x[i])
	}

	aCoeff, bCoeff, cCoeff, dCoeff := gauss.CubicSplineCoefficients(x, y)
	fmt.Println("Коэффициенты сплайна:")
	for i := range aCoeff {
		fmt.Printf("Шаг %d: a:%.5f, b:%.5f, c:%.5f, d:%.5f\n", i, aCoeff[i], bCoeff[i], cCoeff[i], dCoeff[i])
	}

	/*	fmt.Println("a:", aCoeff)
		fmt.Println("b:", bCoeff)
		fmt.Println("c:", cCoeff)
		fmt.Println("d:", dCoeff)*/

	fmt.Println("Значения в промежуточных точках:")
	for i := 1; i < n; i++ {
		k := a + (float64(i)-0.5)*h
		splineValue := gauss.Evaluate(x, aCoeff, bCoeff, cCoeff, dCoeff, k)
		fmt.Printf("S[%d] = %f\n", i-1, splineValue)
	}

	var xQuery float64
	fmt.Print("Введите произвольную точку: ")
	fmt.Scan(&xQuery)

	fmt.Println("Значение сплайна:", gauss.Evaluate(x, aCoeff, bCoeff, cCoeff, dCoeff, xQuery))
	fmt.Println("Значение функции:", f(xQuery))
}
