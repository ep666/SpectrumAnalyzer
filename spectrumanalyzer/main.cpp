#include <QtWidgets/QApplication>
#include "spectrumanalyzer.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    SpectrumAnalyzer sa;
    sa.resize(1024, 768);
    sa.show();
    return a.exec();
}
