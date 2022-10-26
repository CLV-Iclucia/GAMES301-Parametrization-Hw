#include "MeshParamWidget.h"
#include <iostream>

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateLayout();
}

MeshParamWidget::~MeshParamWidget()
{
}

void MeshParamWidget::CreateTabWidget(void)
{
	pbPrintInfo = new QPushButton(tr("Print Mesh Info"));
	tutteParamCirBoundaryUniform = new QPushButton(tr("Tutte Parametrization-Circular Boundary-Uniform Weight"));
	tutteParamCirBoundaryFloater = new QPushButton(tr("Tutte Parametrization-Circular Boundary-Floater Weight"));
	tutteParamSqrBoundaryUniform = new QPushButton(tr("Tutte Parametrization-Square Boundary-Uniform Weight"));
	tutteParamSqrBoundaryFloater = new QPushButton(tr("Tutte Parametrization-Square Boundary-Floater Weight"));
	connect(pbPrintInfo, SIGNAL(clicked()), SIGNAL(PrintInfoSignal()));
	connect(tutteParamCirBoundaryUniform, SIGNAL(clicked()), SIGNAL(PerformTutteSignalCirUniform()));
	connect(tutteParamCirBoundaryFloater, SIGNAL(clicked()), SIGNAL(PerformTutteSignalCirFloater()));
	connect(tutteParamSqrBoundaryUniform, SIGNAL(clicked()), SIGNAL(PerformTutteSignalSqrUniform()));
	connect(tutteParamSqrBoundaryFloater, SIGNAL(clicked()), SIGNAL(PerformTutteSignalSqrFloater()));
	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(pbPrintInfo);
	layout->addWidget(tutteParamCirBoundaryUniform);
	layout->addWidget(tutteParamCirBoundaryFloater);
	layout->addWidget(tutteParamSqrBoundaryUniform);
	layout->addWidget(tutteParamSqrBoundaryFloater);
	layout->addStretch();
	wParam = new QWidget();
	wParam->setLayout(layout);
	saParam = new QScrollArea();
	saParam->setFocusPolicy(Qt::NoFocus);
	saParam->setFrameStyle(QFrame::NoFrame);
	saParam->setWidget(wParam);
	saParam->setWidgetResizable(true);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	this->setLayout(layout);
}
