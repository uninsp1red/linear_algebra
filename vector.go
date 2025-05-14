package vector

import (
	"errors"
	"log"
	"math"
)

type Vector struct {
	data []float64
}

func (v Vector) Data() []float64 {
	return v.data
}

func (v *Vector) Add(other Vector) {
	if len(v.data) != len(other.data) {
		log.Printf("vectors dims: (%d != %d)",
			len(v.data), len(other.data))
		return
	}

	for i := range v.data {
		v.data[i] += other.data[i]
	}
}

func (v *Vector) Minus(other Vector) {
	if len(v.data) != len(other.data) {
		log.Printf("vectors dims: (%d != %d)",
			len(v.data), len(other.data))
		return
	}

	for i := range v.data {
		v.data[i] -= other.data[i]
	}
}
func (v *Vector) Scale(s float64) {
	for i := range v.data {
		v.data[i] *= s
	}
}

func (v Vector) Copy() Vector {
	data := make([]float64, len(v.data))
	copy(data, v.data)
	return Vector{data: data}
}

func (v Vector) ScaleCopy(s float64) Vector {
	newVector := v.Copy()
	(&newVector).Scale(s)
	return newVector
}

func (v Vector) Dot(other Vector) (float64, error) {
	if len(v.data) != len(other.data) {
		return 0, errors.New("vectors len are different")
	}

	var result float64

	for i := range v.data {
		result += v.data[i] * other.data[i]
	}

	return result, nil
}

func (v Vector) Len() int {
	return len(v.data)
}

func NewVector(elem ...float64) Vector {
	return Vector{data: elem}
}

func (v Vector) Project(other Vector) (Vector, error) {
	dot, err := v.Dot(other)
	if err != nil {
		return Vector{}, err
	}

	otherNorm, err := other.Dot(other)
	if err != nil {
		return Vector{}, err
	}

	scale := dot / otherNorm

	return other.ScaleCopy(scale), nil
}

func (v Vector) Norm() float64 {
	var sum float64
	for _, val := range v.data {
		sum += val * val
	}
	return math.Sqrt(sum)
}

func (v Vector) isEqual(other Vector) (bool, error) {
	if v.Len() != other.Len() {
		return false, errors.New("diff dims")
	}
	dim := v.Len()
	for i := 0; i < dim; i++ {
		if v.data[i] != other.data[i] {
			return false, nil
		}
	}
	return true, nil
}
