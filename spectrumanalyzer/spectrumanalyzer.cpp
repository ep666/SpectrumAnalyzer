#include "spectrumanalyzer.h"

#include <QFile>
#include <QFileDialog>
#include <QGridLayout>
#include <QString>
#include <QVector>

#include "dft.h"

SpectrumAnalyzer::SpectrumAnalyzer(QWidget* parent)
	: QWidget(parent), m_isComplex(false)
{
	m_pleFileName = new QLineEdit("File name");
	m_pleFileName->setReadOnly(true);

	m_pcmdAddFile = new QPushButton("Add file");
	QObject::connect(m_pcmdAddFile, SIGNAL(clicked()), SLOT(addFile()));

	m_pcmdHilbertTransform = new QPushButton("Hilbert transform");
	QObject::connect(m_pcmdHilbertTransform, SIGNAL(clicked()), SLOT(hilbertTransform()));

	m_pcbSelectWindow = new QComboBox();
	m_pcbSelectWindow->addItem("None");
	m_pcbSelectWindow->addItem("Triangle");
	m_pcbSelectWindow->addItem("Hamming");
	m_pcbSelectWindow->addItem("Blackman");
	m_pcbSelectWindow->addItem("Hann");

	m_pcmdApplyWindow = new QPushButton("Apply window");
	QObject::connect(m_pcmdApplyWindow, SIGNAL(clicked()), SLOT(addWindow()));

	m_pcpBaseSignal = new QCustomPlot();
	m_pcpBaseSignal->addGraph();
	m_pcpBaseSignal->plotLayout()->insertRow(0);
	m_pcpBaseSignal->plotLayout()->addElement(0, 0,
		new QCPTextElement(m_pcpBaseSignal, "Base signal", QFont("arial", 8, QFont::Courier)));
	m_pcpBaseSignal->replot();

	m_pcpCurrentSignal = new QCustomPlot();
	m_pcpCurrentSignal->addGraph();
	m_pcpCurrentSignal->addGraph();
	m_pcpCurrentSignal->plotLayout()->insertRow(0);
	m_pcpCurrentSignal->plotLayout()->addElement(0, 0,
		new QCPTextElement(m_pcpCurrentSignal, "Current signal", QFont("arial", 8, QFont::Courier)));
	m_pcpCurrentSignal->replot();

	m_pcpSpectrum = new QCustomPlot();
	m_pcpSpectrum->addGraph();
	m_pcpSpectrum->plotLayout()->insertRow(0);
	m_pcpSpectrum->plotLayout()->addElement(0, 0,
		new QCPTextElement(m_pcpSpectrum, "Spectrum", QFont("arial", 8, QFont::Courier)));
	m_pcpCurrentSignal->replot();


	m_psbFrequency = new QDoubleSpinBox();
	m_psbFrequency->setMinimum(1);
	m_psbFrequency->setMaximum(100000);
	QLabel* freq = new QLabel("f, Hz");

	m_psbSampleRate = new QDoubleSpinBox();
	m_psbSampleRate->setMinimum(1);
	m_psbSampleRate->setMaximum(1000000);
	QLabel* sampleRate = new QLabel("Sample rate, Hz");

	QLabel* heightWindow = new QLabel("Window:");

	QGridLayout* layout = new QGridLayout;
	QVBoxLayout* vlay = new QVBoxLayout;
	vlay->addWidget(sampleRate);
	vlay->addWidget(m_psbSampleRate);
	vlay->addWidget(freq);
	vlay->addWidget(m_psbFrequency);
	vlay->addWidget(heightWindow);
	vlay->addWidget(m_pcbSelectWindow);
	vlay->addWidget(m_pcmdApplyWindow);
	vlay->addWidget(m_pcmdHilbertTransform);

	layout->addWidget(m_pleFileName, 0, 0);
	layout->addWidget(m_pcmdAddFile, 0, 1);
	layout->addWidget(m_pcpBaseSignal, 1, 0);
	layout->addWidget(m_pcpCurrentSignal, 2, 0);
	layout->addLayout(vlay, 1, 1, 3, 1, Qt::AlignVCenter);
	layout->addWidget(m_pcpSpectrum, 3, 0);

	setLayout(layout);

	setMinimumSize(640, 480);
}

void SpectrumAnalyzer::addFile()
{
	m_baseSignal.clear();
	m_currentSignal.clear();
	m_spectrum.clear();
	m_pcmdHilbertTransform->setEnabled(true);
	m_isComplex = false;

	QString fileName = QFileDialog::getOpenFileName(this, "Open the file");
	m_pleFileName->setText(fileName);
	QFile file(fileName);
	if (!file.open(QIODevice::ReadOnly | QFile::Text)) {
		QMessageBox::warning(this, "Warning", "Cannot open file : " + file.errorString());
		file.close();
		return;
	}

	while (!file.atEnd())
	{
		m_baseSignal.push_back(file.readLine().toDouble());
	}

	file.close();

	printBaseSignal();

	m_currentSignal = m_baseSignal;

	prepareCurrentSignal();
	printCurrentSignal();

	m_spectrum = m_currentSignal;

	fftAndPrint(m_isComplex);
}

void SpectrumAnalyzer::printBaseSignal()
{
	QVector<double> x, yReal, yImag;

	for (int i = 0; i < m_baseSignal.size(); ++i)
	{
		x.push_back(i);
		yReal.push_back(m_baseSignal[i].real());
		yImag.push_back(m_baseSignal[i].imag());
	}

	m_pcpBaseSignal->graph(0)->setData(x, yReal);
	m_pcpBaseSignal->graph(0)->setPen(QPen(Qt::blue));
	m_pcpBaseSignal->graph(0)->rescaleAxes();
	m_pcpBaseSignal->replot();
}

void SpectrumAnalyzer::prepareCurrentSignal()
{
	while (!isPowerOfTwo(m_currentSignal.size()))
	{
		m_currentSignal.push_back(0);
	}
}

void SpectrumAnalyzer::printCurrentSignal()
{
	QVector<double> x, yReal, yImag;

	for (int i = 0; i < m_currentSignal.size(); ++i)
	{
		double factor = 1 / (m_psbSampleRate->value());
		x.push_back(i * factor);
		yReal.push_back(m_currentSignal[i].real());
		yImag.push_back(m_currentSignal[i].imag());
	}

	m_pcpCurrentSignal->graph(0)->setData(x, yReal);
	m_pcpCurrentSignal->graph(0)->setPen(QPen(Qt::blue));
	m_pcpCurrentSignal->graph(0)->rescaleAxes();
	m_pcpCurrentSignal->graph(1)->setData(x, yImag);
	m_pcpCurrentSignal->graph(1)->setPen(QPen(Qt::red));
	m_pcpCurrentSignal->graph(1)->rescaleAxes(true);
	m_pcpCurrentSignal->replot();
}

void SpectrumAnalyzer::fftAndPrint(bool isComplex, bool isTransformed)
{
	if (!isTransformed)
		dsp::fft(m_spectrum);

	int n = m_spectrum.size();
	QVector<double> x, y;
	double bps = m_psbSampleRate->value() / m_psbFrequency->value();
	double period = m_currentSignal.size() / bps;
	double xValue = m_psbFrequency->value() / period;

	if (!isComplex)
	{
		for (int i = n / 2; i < n; ++i)
		{
			x.push_back((i - n) * xValue);
			y.push_back(std::abs(m_spectrum[i]) / n);
		}
		for (int i = 0; i < n / 2; ++i)
		{
			x.push_back(i * xValue);
			y.push_back(std::abs(m_spectrum[i]) / n);
		}
	}
	else
		for (int i = 0; i < n; ++i)
		{
			x.push_back(i * xValue);
			y.push_back(std::abs(m_spectrum[i]) / n);
		}

	m_pcpSpectrum->graph(0)->setData(x, y);
	m_pcpSpectrum->graph(0)->rescaleAxes();
	m_pcpSpectrum->replot();
}

bool SpectrumAnalyzer::isPowerOfTwo(int n)
{
	if (n <= 0)
		return false;

	return (ceil(log2(n)) == floor(log2(n)));
}

void SpectrumAnalyzer::addWindow()
{
	m_currentSignal = m_baseSignal;
	prepareCurrentSignal();
	int N = m_currentSignal.size();
	switch (m_pcbSelectWindow->currentIndex())
	{
	case 0:
		printCurrentSignal();
		m_spectrum = m_currentSignal;
		fftAndPrint(m_isComplex);
		break;

	case 1:
		for (int n = 0; n < N; ++n)
		{
			m_currentSignal[n] *= (1.0 * N - 2.0 * std::abs(1.0 * n - N / 2.0)) / N;
		}
		printCurrentSignal();
		m_spectrum = m_currentSignal;
		fftAndPrint(m_isComplex);
		break;

	case 2:
		for (int n = 0; n < N; ++n)
		{
			m_currentSignal[n] *= 0.54 - 0.46 * std::cos(2.0 * M_PI * n / N);
		}
		printCurrentSignal();
		m_spectrum = m_currentSignal;
		fftAndPrint(m_isComplex);
		break;

	case 3:
		for (int n = 0; n < N; ++n)
		{
			m_currentSignal[n] *= 0.42 - 0.5 * std::cos(2.0 * M_PI * n / N) + 0.08 * std::cos(4.0 * M_PI * n / N);
		}
		printCurrentSignal();
		m_spectrum = m_currentSignal;
		fftAndPrint(m_isComplex);
		break;

	case 4:
		for (int n = 0; n < N; ++n)
		{
			m_currentSignal[n] *= 0.5 - 0.5 * cos(2.0 * M_PI * n / N);
		}
		printCurrentSignal();
		m_spectrum = m_currentSignal;
		fftAndPrint(m_isComplex);
		break;
	}
}

void SpectrumAnalyzer::hilbertTransform()
{
	using namespace std::complex_literals;
	int filterOrder = 128;

	m_pcmdHilbertTransform->setEnabled(false);
	m_isComplex = true;

	std::vector<std::complex<double>> hilbertCore;

	for (int k = -filterOrder; k <= filterOrder; ++k)
	{
		if (k != 0)
			hilbertCore.push_back(1.0 / (M_PI * k) * (1 - std::cos(M_PI * k)));
		else
			hilbertCore.push_back(0);
	}

	int N = m_baseSignal.size();
	int M = hilbertCore.size();

	std::vector<std::complex<double>> tempSignal = m_baseSignal;
	std::vector<std::complex<double>> out;

	if (N < 2)
		return;

	for (int i = 0; i < M; ++i)
		tempSignal.push_back(0);

	for (int i = 0; i < N; ++i)
		hilbertCore.push_back(0);

	//Convolution
	out = m_baseSignal;

	for (int m = 0; m < out.size(); ++m)
	{
		for (int n = m; n < out.size(); ++n)
		{
			out[n] += tempSignal[m] * 1i * hilbertCore[n - m];
		}
	}

	m_baseSignal = out;
	m_currentSignal = out;

	prepareCurrentSignal();
	printCurrentSignal();

	m_spectrum = m_currentSignal;
	fftAndPrint(m_isComplex);
}
