#pragma once

#include <QPushButton>
#include <QComboBox>
#include <QLineEdit>
#include <QTextEdit>
#include <QSpinBox>
#include "qcustomplot.h"
#include <vector>
#include <complex>

class SpectrumAnalyzer : public QWidget
{
    Q_OBJECT

public:
    SpectrumAnalyzer(QWidget *parent = Q_NULLPTR);
    using complex = std::complex<double>;

private:
    QLineEdit*              m_pleFileName;
    QPushButton*            m_pcmdAddFile;
    QPushButton*            m_pcmdApplyWindow;
    QPushButton*            m_pcmdHilbertTransform;
    QComboBox*              m_pcbSelectWindow;
    QCustomPlot*            m_pcpBaseSignal;
    QCustomPlot*            m_pcpCurrentSignal;
    QCustomPlot*            m_pcpSpectrum;
    QDoubleSpinBox*         m_psbFrequency;
    QDoubleSpinBox*         m_psbSampleRate;
    std::vector<complex>    m_baseSignal;
    std::vector<complex>    m_currentSignal;
    std::vector<complex>    m_spectrum;
    bool                    m_isComplex;

    void printBaseSignal();
    void prepareCurrentSignal();
    void printCurrentSignal();
    void fftAndPrint(bool isComplex = false, bool isTransformed = false);
    bool isPowerOfTwo(int n);
    

private slots:
    void addFile();
    void addWindow();
    void hilbertTransform();
};
