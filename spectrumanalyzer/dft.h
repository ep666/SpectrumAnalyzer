#pragma once
#include <vector>
#include <complex>

namespace dsp
{
	std::vector<std::complex<double>> dft(const std::vector<double>& signal);

	std::vector<std::complex<double>> dft(const std::vector<std::complex<double>>& signal);

	std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& spectrum);

	void fftRecursive(std::vector<std::complex<double>>& signal);

	void fft(std::vector<std::complex<double>>& signal);

	void ifft(std::vector<std::complex<double>>& spectrum);
}