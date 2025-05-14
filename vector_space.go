package vector_space

import (
	"errors"
	"linear_shell/vector"
	"math"
	"sort"
)

type Basis struct {
	vectors []vector.Vector
	dim     int
}

func NewBasis(vectors ...vector.Vector) (Basis, error) {
	if len(vectors) == 0 {
		return Basis{}, errors.New("empty basis")
	}

	b := Basis{
		vectors: vectors,
		dim:     len(vectors),
	}

	if ind, err := b.AreLinearlyIndependent(); err != nil || len(ind) != b.dim {
		return Basis{}, errors.New("vectors are not independent")
	}

	return b, nil
}

func (b Basis) Dimension() int {
	return b.dim
}

func (b Basis) Vectors() []vector.Vector {
	return b.vectors
}

func (b Basis) AreLinearlyIndependent() ([]int, error) {
	n := b.dim
	vectorLen := b.vectors[0].Len()
	matrix := make([][]float64, n)
	for i, v := range b.vectors {
		matrix[i] = make([]float64, vectorLen)
		copy(matrix[i], v.Data())
	}

	indexes := make([]int, n)
	for i := range indexes {
		indexes[i] = i
	}

	independentIdx := []int{}
	rank := 0

	for col := 0; col < vectorLen && rank < n; col++ {
		maxRow := rank
		for i := rank + 1; i < n; i++ {
			if math.Abs(matrix[i][col]) > math.Abs(matrix[maxRow][col]) {
				maxRow = i
			}
		}

		if math.Abs(matrix[maxRow][col]) < 1e-10 {
			continue
		}

		matrix[rank], matrix[maxRow] = matrix[maxRow], matrix[rank]
		indexes[rank], indexes[maxRow] = indexes[maxRow], indexes[rank]

		independentIdx = append(independentIdx, indexes[rank])
		pivot := matrix[rank][col]
		for j := col; j < vectorLen; j++ {
			matrix[rank][j] /= pivot
		}

		for i := 0; i < n; i++ {
			if i != rank {
				factor := matrix[i][col]
				for j := col; j < vectorLen; j++ {
					matrix[i][j] -= factor * matrix[rank][j]
				}
			}
		}

		rank++
	}

	sort.Ints(independentIdx)
	return independentIdx, nil
}

func (b Basis) Orthogonalize() (Basis, error) {
	if b.dim == 0 {
		return Basis{}, errors.New("empty basis")
	}

	orthogonalVectors := make([]vector.Vector, b.dim)
	orthogonalVectors[0] = b.vectors[0].Copy()

	for i := 1; i < b.dim; i++ {
		orthogonalVectors[i] = b.vectors[i].Copy()

		for j := 0; j < i; j++ {
			projection, err := orthogonalVectors[i].Project(orthogonalVectors[j])
			if err != nil {
				return Basis{}, err
			}
			orthogonalVectors[i].Minus(projection)
		}
	}

	return NewBasis(orthogonalVectors...)
}
