#define _USE_MATH_DEFINES

#include "dft.h"
#include <cmath>
#include <algorithm>

std::vector<std::complex<double>> dsp::dft(const std::vector<double>& signal)
{
	using namespace std::complex_literals;

	int N = signal.size();
	std::vector<std::complex<double>> out;
	if (N == 0)
		return out;

	for (int k = 0; k < N; ++k)
	{
		std::complex<double> temp = 0;
		for (int n = 0; n < N; ++n)
		{
			temp += signal[n] * std::exp(-1i * ((2.0 * M_PI) / N * n * k));
		}
		out.push_back(temp);
	}
	return out;
}

std::vector<std::complex<double>> dsp::dft(const std::vector<std::complex<double>>& signal)
{
	using namespace std::complex_literals;

	int N = signal.size();
	std::vector<std::complex<double>> out;
	if (N == 0)
		return out;

	for (int k = 0; k < N; ++k)
	{
		std::complex<double> temp = 0;
		for (int n = 0; n < N; ++n)
		{
			temp += signal[n] * std::exp(-1i * ((2.0 * M_PI) / N * n * k));
		}
		out.push_back(temp);
	}
	return out;
}

std::vector<std::complex<double>> dsp::idft(const std::vector<std::complex<double>>& spectrum)
{
	using namespace std::complex_literals;

	std::vector<std::complex<double>> out;
	int N = spectrum.size();
	if (N == 0)
		return out;

	for (int k = 0; k < N; ++k)
	{
		std::complex<double> temp;
		for (int n = 0; n < N; ++n)
		{
			temp += spectrum[n] * std::exp(1i * ((2.0 * M_PI) / N * n * k));
		}
		out.push_back(temp * (1.0 / N));
	}

	return out;
}

void dsp::fftRecursive(std::vector<std::complex<double>>& signal)
{
	if (signal.size() <= 1)
		return;

	int N = signal.size();

	using namespace std::complex_literals;
	using complex = std::vector<std::complex<double>>;

	complex odd;
	complex even;

	for (int i = 0; i < N / 2; ++i)
	{
		even.push_back(signal[2 * i]);
		odd.push_back(signal[2 * i + 1]);
	}

	fftRecursive(odd);
	fftRecursive(even);

	for (int k = 0; k < N / 2; ++k)
	{
		std::complex<double> temp = odd[k] * std::exp(-1i * (k * 2.0 * M_PI / N));

		signal[k] = even[k] + temp;
		signal[k + N / 2] = even[k] - temp;
	}
}

void dsp::fft(std::vector<std::complex<double>>& signal)
{
	int N = signal.size();
	if (N < 2)
		return;
	using namespace std::complex_literals;
	using complex = std::complex<double>;

	//Bit-reversal permutation
	for (int i = 1, j = 0; i < N; ++i)
	{
		int bit = N >> 1;
		for (; j >= bit; bit >>= 1)
			j -= bit;

		j += bit;
		if (i < j)
			std::swap(signal[i], signal[j]);
	}

	for (int len = 2; len <= N; len <<= 1)
	{
		for (int i = 0; i < N; i += len)
		{
			for (int j = 0; j < len / 2; ++j)
			{
				complex wn = std::exp(-1i * (j * 2.0 * M_PI / len));
				complex even = signal[i + j];
				complex odd = signal[i + j + len / 2];
				signal[i + j] = even + odd * wn;
				signal[i + j + len / 2] = even - odd * wn;
			}
		}
	}
}

void dsp::ifft(std::vector<std::complex<double>>& spectrum)
{
	int N = spectrum.size();
	if (N < 2)
		return;
	using namespace std::complex_literals;
	using complex = std::complex<double>;

	//Bit-reversal permutation
	for (int i = 1, j = 0; i < N; ++i)
	{
		int bit = N >> 1;
		for (; j >= bit; bit >>= 1)
			j -= bit;

		j += bit;
		if (i < j)
			std::swap(spectrum[i], spectrum[j]);
	}

	for (int len = 2; len <= N; len <<= 1)
	{
		for (int i = 0; i < N; i += len)
		{
			for (int j = 0; j < len / 2; ++j)
			{
				complex wn = std::exp(1i * (j * 2.0 * M_PI / len));
				complex even = spectrum[i + j];
				complex odd = spectrum[i + j + len / 2];
				spectrum[i + j] = (even + odd * wn) / 2.0;
				spectrum[i + j + len / 2] = (even - odd * wn) / 2.0;
			}
		}
	}
}
