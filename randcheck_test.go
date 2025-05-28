package randcheck

import (
	"gonum.org/v1/gonum/mathext"
	"math"
	"testing"
)

func equalFloat64(a, b float64) bool {
	tolerance := 1e-6 // Define a small tolerance value
	return math.Abs(a-b) < tolerance
}
func TestMonobit(t *testing.T) {

	bits := BitString("1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000")

	pval := MonoBit(bits)

	want := 0.109599

	if !equalFloat64(pval, want) {
		t.Fail()
	}
}

func TestFrequency(t *testing.T) {

	bits := BitString("1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000")

	pval := Frequency(bits, 10)

	want := 0.706438

	if !equalFloat64(pval, want) {
		t.Fail()
	}
}

func TestRuns(t *testing.T) {

	bits := BitString("1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000")

	pval := Runs(bits)

	want := 0.500798

	if !equalFloat64(pval, want) {
		t.Fail()
	}
}

func TestLongestRuns(t *testing.T) {

	bits := BitString("11001100000101010110110001001100111000000000001001001101010100010001001111010110100000001101011111001100111001101101100010110010")
	pval := LongestRunOfOnes8(bits)

	want := 0.180598
	//want := 0.180609

	if !equalFloat64(pval, want) {
		t.Fail()
	}
}

func TestSerial1(t *testing.T) {
	bits := BitString("0011011101")
	pval1, pval2 := Serial(bits, 3)

	want1 := 0.8087921354109989
	want2 := 0.6703200460356394

	// example is incorrect
	//want1 := 0.9057
	//want2 := 0.8805

	if !equalFloat64(pval1, want1) {
		t.Fail()
	}

	if !equalFloat64(pval2, want2) {
		t.Fail()
	}
}

func TestSerial2(t *testing.T) {
	const m = 2

	del_psisqr := 0.339764
	delsqr_psisqr := 0.336400

	K := math.Pow(2, float64(m-1))
	pval1 := mathext.GammaIncRegComp(0.5*K, 0.5*del_psisqr)

	K = math.Pow(2, float64(m-2))
	pval2 := mathext.GammaIncRegComp(0.5*K, 0.5*delsqr_psisqr)

	want1 := 0.843764
	want2 := 0.561915

	if !equalFloat64(pval1, want1) {
		t.Fail()
	}

	if !equalFloat64(pval2, want2) {
		t.Fail()
	}
}

func TestAppEntropy1(t *testing.T) {
	bits := BitString("0100110101")
	pval := AppEntropy(bits, 3)

	want := 0.261961

	if !equalFloat64(pval, want) {
		t.Fail()
	}
}

func TestAppEntropy2(t *testing.T) {
	bits := BitString("1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000")

	pval := AppEntropy(bits, 2)

	want := 0.235301

	if !equalFloat64(pval, want) {
		t.Fail()
	}
}
