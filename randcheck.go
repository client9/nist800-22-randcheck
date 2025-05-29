package randcheck

import (
	"fmt"
	"math"

	"golang.org/x/exp/constraints"
	"gonum.org/v1/gonum/mathext"
)

func RunAll[T constraints.Integer](b []T) error {
	var pval float64

	pval = MonoBit(b)
	if pval < 0.01 {
		return fmt.Errorf("failed monobit with pval=%f", pval)
	}
	pval = Frequency(b, 3)
	if pval < 0.01 {
		return fmt.Errorf("failed frequency with pval=%f\n", pval)
	}
	pval = Runs(b)
	if pval < 0.01 {
		return fmt.Errorf("failed Runs with pval=%f", pval)
	}
	pval = LongestRunOfOnes8(b)
	if pval < 0.01 {
		return fmt.Errorf("failed LongestRunOfOnes 8 with pval=%f", pval)
	}
	pval1, pval2 := Serial(b, 3)
	if pval1 < 0.01 || pval2 < 0.01 {
		return fmt.Errorf("failed Serial with pval= %f %f", pval1, pval2)
	}
	pval = AppEntropy(b, 3)
	if pval < 0.01 {
		return fmt.Errorf("failed App Entropy with pval=%f", pval)
	}
	return nil

}

func BitString(s string) []byte {
	bits := []byte(s)
	for i, b := range bits {
		switch b {
		case '0':
			bits[i] = 0
		case '1':
			bits[i] = 1
		case 0, 1:
			// NOP
		default:
			panic("bad char in bit string")
		}
	}
	return bits
}

func MonoBit[T constraints.Integer](bits []T) float64 {

	n := float64(len(bits))

	if n < 100 {
		panic("input too small")
	}

	sum := 0
	for _, val := range bits {
		if val == 0 {
			sum--
		} else if val == 1 {
			sum++
		} else {
			panic("unnknwn val")
		}
	}

	sobs := math.Abs(float64(sum)) / math.Sqrt(n)

	pval := math.Erfc(sobs / math.Sqrt(2))

	return pval

}

func Frequency[T constraints.Integer](bits []T, blocksize int) float64 {
	n := len(bits)

	if n < 100 {
		panic("input too small")
	}

	blocks := len(bits) / blocksize

	obs := make([]float64, blocks)

	for i := 0; i < blocks; i++ {
		ones := 0
		for j := 0; j < blocksize; j++ {
			ones += int(bits[i*blocksize+j])
		}
		obs[i] = float64(ones) / float64(blocksize)
	}

	sum := 0.0
	for _, val := range obs {
		sum += math.Pow(val-0.5, 2)
	}

	chisqr := 4.0 * float64(blocksize) * sum

	pval := mathext.GammaIncRegComp(0.5*float64(blocks), 0.5*chisqr)
	return pval

}

func Runs[T constraints.Integer](bits []T) float64 {

	n := len(bits)
	nf := float64(n)

	if n < 100 {
		panic("input too small")
	}

	pi := float64(0)
	for _, val := range bits {
		if val == 1 {
			pi++
		}
	}
	pi = pi / nf

	runs := 1.0
	last := bits[0]
	for i := 1; i < n; i++ {
		current := bits[i]
		if current == last {
			continue
		}
		last = current
		runs++
	}

	pic := pi * (1.0 - pi)
	return math.Erfc((runs - (2.0 * nf * pic)) / (2.0 * math.Sqrt(2.0*nf) * pic))
}

func onesRun[T constraints.Integer](bits []T) int {
	maxRun := 0
	currentRun := 0
	for _, b := range bits {
		if b == 1 {
			currentRun++
			continue
		}
		if currentRun > maxRun {
			maxRun = currentRun
		}
		currentRun = 0
	}
	if currentRun > maxRun {
		maxRun = currentRun
	}
	return maxRun
}

// 2.4 Longest Run of Ones
func LongestRunOfOnes8[T constraints.Integer](bits []T) float64 {
	const M = 8
	n := len(bits)
	N := float64(n) / float64(M)

	ones := make([]int, 0, n/M)
	for i := 0; i < n/M; i += 1 {
		ones = append(ones, onesRun(bits[i*M:i*M+M]))
	}
	// remap
	v0 := 0.0
	v1 := 0.0
	v2 := 0.0
	v3 := 0.0
	for _, val := range ones {
		switch val {
		case 0, 1:
			v0++
		case 2:
			v1++
		case 3:
			v2++
		default:
			v3++
		}
	}
	const K float64 = 3.0

	const pi0 float64 = 0.2148
	const pi1 float64 = 0.3672
	const pi2 float64 = 0.2305
	const pi3 float64 = 0.1875

	x0 := math.Pow(v0-N*pi0, 2) / (pi0)
	x1 := math.Pow(v1-N*pi1, 2) / (pi1)
	x2 := math.Pow(v2-N*pi2, 2) / (pi2)
	x3 := math.Pow(v3-N*pi3, 2) / (pi3)

	chisqr := (x0 + x1 + x2 + x3) / N

	pval := mathext.GammaIncRegComp(0.5*K, 0.5*chisqr)
	return pval
}

// 2.5 Binary Matrix Rank Test
// 40,000 bits min

// 2.6 Discrete Fourier Transform (Spectral) Test
// 1000 bits min

// 2.7 Non-overlapping Template Matching Test

// 2.8 Overlapping Template Matching Test
// 1,000,000 bits

// 2.9 Mauer's Universal Statistic Test
// 1,000,000 bits

// 2.10 Linear Complexity Test
// 1,000,000 bits

// 2.11 Serial Test
// Choose m and n such that m <  log2 n-2.

func serialPsiSqr[T constraints.Integer](bits []T, m int) float64 {

	n := len(bits)
	nf := float64(n)

	bits = append(bits, bits[:m-1]...)
	bins := make([]int, (1 << m))

	for i := 0; i < n; i++ {
		idx := sliceToInt(bits[i : i+m])
		bins[idx] = bins[idx] + 1
	}

	// compute sum of squares
	sumsqr := 0.0
	for _, count := range bins {
		sumsqr += float64(count * count)
	}
	psisqr := (math.Pow(2, float64(m))/nf)*sumsqr - nf
	return psisqr
}

func Serial[T constraints.Integer](bits []T, m int) (float64, float64) {
	psisqr_m := serialPsiSqr(bits, m)
	psisqr_m1 := serialPsiSqr(bits, m-1)
	psisqr_m2 := serialPsiSqr(bits, m-2)

	del_psisqr := psisqr_m - psisqr_m1
	delsqr_psisqr := psisqr_m - 2.0*psisqr_m1 + psisqr_m2

	K := math.Pow(2, float64(m-1))
	pval1 := mathext.GammaIncRegComp(0.5*K, 0.5*del_psisqr)

	K = math.Pow(2, float64(m-2))
	pval2 := mathext.GammaIncRegComp(0.5*K, 0.5*delsqr_psisqr)

	return pval1, pval2
}

// 2.12 Approximate Entropy Test

func slicePrint[T constraints.Integer](bits []T) string {
	buf := make([]byte, len(bits), len(bits))
	for i := 0; i < len(bits); i++ {
		if bits[i] == 0 {
			buf[i] = '0'
		} else {
			buf[i] = '1'
		}
	}
	return string(buf)
}

func sliceToInt[T constraints.Integer](bits []T) int {
	val := bits[0]
	for i := 1; i < len(bits); i++ {
		val = val<<1 | bits[i]
	}
	return int(val)
}

// bits, bits to test
// m, blocksize to test
func appEntropyPhi[T constraints.Integer](bits []T, blocksize int) float64 {
	m := blocksize
	n := len(bits)

	bits = append(bits, bits[:m-1]...)
	bins := make([]int, (1 << blocksize))

	for i := 0; i < n; i++ {
		idx := sliceToInt(bits[i : i+m])
		bins[idx] = bins[idx] + 1
	}

	C := make([]float64, len(bins))
	for i := 0; i < len(bins); i++ {
		C[i] = float64(bins[i]) / float64(n)
	}
	phi := float64(0.0)
	for _, c := range C {

		// not part of the spec or example
		// document's example uses "0 Log 0"
		if c == 0 {
			continue
		}

		phi += c * math.Log(c)
	}
	//fmt.Printf("Phi = %f\n", phi)
	return phi
}

func appEntropyM[T constraints.Integer](bits []T, m int) float64 {
	ae := appEntropyPhi(bits, m) - appEntropyPhi(bits, m+1)

	// this can dip below zero, maybe due to rounding errors
	if ae <= 0 {
		return 0.0
	}
	return ae
}

func AppEntropy[T constraints.Integer](bits []T, m int) float64 {

	// TODO constraints
	// m= 3 ==> len(b) >= 512

	n := float64(len(bits))

	ApEn := appEntropyM(bits, m)

	chisqr := 2.0 * n * (math.Log(2.0) - ApEn)
	K := float64(int(1) << m) // i.e. 2^blocksize
	pval := mathext.GammaIncRegComp(0.5*K, 0.5*chisqr)
	return pval
}

// 2.13 Cumulative Sums Test

// 2.14 Random Excursions Test
// 1,000,000 bits

// 2.15 Random Excursions Variant Test
// 1,000,000 bits
