#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "../PolyMesh/include/PolyMesh/IOManager.h"
#include "Numerical.h"
using VertPtr = acamcad::polymesh::MVert*;
using HalfEdgePtr = acamcad::polymesh::MHalfedge*;
using acamcad::MPoint3;
using acamcad::MVector3;

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
public:
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	void Clear(void);
	void UpdateMesh(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
	void PerformTutte(bool, bool);
	void TutteParam(void);
signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
protected:
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);

private:
	void FixBoundaryCircle(std::vector<MPoint3>&) const;
	void FixBoundarySquare(std::vector<MPoint3>&) const;
	std::vector<Real> CalcUniformWeight(VertPtr) const;
	std::vector<Real> CalcFloaterWeight(VertPtr) const;
	void PolarMap(std::vector<VertPtr>&, VertPtr) const;
	void DrawPoints(void) const;
	void DrawWireframe(void) const;
	void DrawHiddenLines(void) const;
	void DrawFlatLines(void) const;
	void DrawFlat(void) const;
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;
protected:
	acamcad::polymesh::PolyMesh* polyMesh = new acamcad::polymesh::PolyMesh();
	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	acamcad::MPoint3 ptMin;
	acamcad::MPoint3 ptMax;
	
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;
	bool isBoundaryCircle;
	bool isFloaterWeight;
};
