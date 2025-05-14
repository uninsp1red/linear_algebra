package linear_span

import (
	"errors"
	"fmt"
	"linear_shell/vector"
	"linear_shell/vector_space"
	"math"
	"math/rand"
)

type LinearSpan struct {
	vectors []vector.Vector
	basis   vector_space.Basis
}

func NewLinearSpan(vectors ...vector.Vector) (*LinearSpan, error) {
	if len(vectors) == 0 {
		return nil, errors.New("empty set of vectors")
	}

	vectorsCopy := make([]vector.Vector, len(vectors))
	for i, v := range vectors {
		vectorsCopy[i] = v.Copy()
	}

	span := &LinearSpan{
		vectors: vectorsCopy,
	}

	basis, err := span.computeBasis()
	if err != nil {
		return nil, err
	}

	span.basis = basis
	return span, nil
}

func (ls *LinearSpan) computeBasis() (vector_space.Basis, error) {
	tempBasis, err := vector_space.NewBasis(ls.vectors...)
	if err != nil {
		return vector_space.Basis{}, err
	}

	independentIndexes, err := tempBasis.AreLinearlyIndependent()
	if err != nil {
		return vector_space.Basis{}, err
	}

	independentVectors := make([]vector.Vector, len(independentIndexes))
	for i, idx := range independentIndexes {
		independentVectors[i] = ls.vectors[idx]
	}

	return vector_space.NewBasis(independentVectors...)
}

func (ls LinearSpan) Vectors() []vector.Vector {
	return ls.vectors
}

func (ls LinearSpan) Basis() vector_space.Basis {
	return ls.basis
}

func (ls LinearSpan) ContainsVector(v vector.Vector) (bool, error) {
	if v.Len() != ls.vectors[0].Len() {
		return false, errors.New("vector dimensions don't match")
	}

	if ls.basis.Dimension() == v.Len() {
		return true, nil
	}

	dim := ls.basis.Dimension()
	vectorLen := v.Len()
	otherMatrix := make([][]float64, dim)

	basisVectors := ls.basis.Vectors()
	for i := 0; i < dim; i++ {
		otherMatrix[i] = make([]float64, vectorLen+1)
		copy(otherMatrix[i], basisVectors[i].Data())
		otherMatrix[i][vectorLen] = v.Data()[i]
	}

	for i := 0; i < dim; i++ {
		maxRow := i
		for j := i + 1; j < dim; j++ {
			if math.Abs(otherMatrix[j][i]) > math.Abs(otherMatrix[maxRow][i]) {
				maxRow = j
			}
		}

		if math.Abs(otherMatrix[maxRow][i]) < 1e-10 {
			continue
		}

		otherMatrix[i], otherMatrix[maxRow] = otherMatrix[maxRow], otherMatrix[i]

		p := otherMatrix[i][i]
		for j := i; j <= vectorLen; j++ {
			otherMatrix[i][j] /= p
		}

		for j := 0; j < dim; j++ {
			if j != i {
				factor := otherMatrix[j][i]
				for k := i; k <= vectorLen; k++ {
					otherMatrix[j][k] -= factor * otherMatrix[i][k]
				}
			}
		}
	}

	for i := 0; i < dim; i++ {
		isZeroRow := true
		for j := 0; j < vectorLen; j++ {
			if math.Abs(otherMatrix[i][j]) > 1e-10 {
				isZeroRow = false
				break
			}
		}

		if isZeroRow && math.Abs(otherMatrix[i][vectorLen]) > 1e-10 {
			return false, nil
		}
	}

	return true, nil
}

// для ввода сообственных коэф
func (ls LinearSpan) GenVecFromCoeffs(coeffs ...float64) (vector.Vector, error) {
	dim := ls.basis.Dimension()

	if dim != len(coeffs) {
		return vector.Vector{}, errors.New("too many coeffs")
	}
	zeroData := make([]float64, ls.vectors[0].Len())

	basisVectors := ls.basis.Vectors()
	resultVector := vector.NewVector(zeroData...)

	for i := 0; i < dim; i++ {
		tempVector := basisVectors[i].ScaleCopy(coeffs[i])
		resultVector.Add(tempVector)
	}
	return resultVector, nil
}
func (ls LinearSpan) GenerateVector() vector.Vector {
	dim := ls.basis.Dimension()
	coeffs := make([]float64, dim)

	for i := 0; i < dim; i++ {
		coeffs[i] = rand.Float64()*2 - 1
	}

	basisVectors := ls.basis.Vectors()
	resultData := make([]float64, basisVectors[0].Len())

	for i := 0; i < dim; i++ {
		for j := 0; j < len(resultData); j++ {
			resultData[j] += coeffs[i] * basisVectors[i].Data()[j]
		}
	}

	return vector.NewVector(resultData...)
}

func (ls LinearSpan) GetOrthogonalBasis() (vector_space.Basis, error) {
	return ls.basis.Orthogonalize()
}

func (ls LinearSpan) PrintLSBasis() {
	for _, v := range ls.basis.Vectors() {
		fmt.Println(v)
	}
}
