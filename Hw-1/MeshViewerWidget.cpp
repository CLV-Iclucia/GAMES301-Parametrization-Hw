#include <QtCore>
#include "MeshViewerWidget.h"

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0,0,0),
	ptMax(0,0,0),
	isEnableLighting(true),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false)
{
}

MeshViewerWidget::~MeshViewerWidget(void)
{
}

bool MeshViewerWidget::LoadMesh(const std::string & filename)
{
	Clear();

	bool read_OK = acamcad::polymesh::loadMesh(filename, polyMesh);
	std::cout << "Load mesh from file " << filename << std::endl;
	if (read_OK)
	{
		strMeshFileName = QString::fromStdString(filename);
		QFileInfo fi(strMeshFileName);
		strMeshPath = fi.path();
		strMeshBaseName = fi.baseName();
		UpdateMesh();
		update();
		return true;
	}
	return false;
}

void MeshViewerWidget::Clear(void)
{
	polyMesh->clear();
}

void MeshViewerWidget::UpdateMesh(void)
{
	polyMesh->updateFacesNormal();
	polyMesh->updateMeshNormal();
	polyMesh->updateVerticesNormal();
	if (polyMesh->numVertices() == 0)
	{
		std::cerr << "ERROR: UpdateMesh() No vertices!" << std::endl;
		return;
	}
	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;

	for (const auto& vh : polyMesh->vertices())
	{
		auto p = vh->position();
		for (size_t i = 0; i < 3; i++)
		{
			ptMin[i] = ptMin[i] < p[i] ? ptMin[i] : p[i];
			ptMax[i] = ptMax[i] > p[i] ? ptMax[i] : p[i];
		}
	}

	double avelen = 0.0;
	double maxlen = 0.0;
	double minlen = DBL_MAX;
	for (const auto& eh : polyMesh->edges()) {
		double len = eh->length();
		maxlen = len > maxlen ? len : maxlen;
		minlen = len < minlen ? len : minlen;
		avelen += len;
	}

	SetScenePosition((ptMin + ptMax)*0.5, (ptMin - ptMax).norm()*0.5);
	std::cout << "Information of the input mesh:" << std::endl;
	std::cout << "  [V, E, F] = [" << polyMesh->numVertices()<< ", " << polyMesh->numEdges() << ", " << polyMesh->numPolygons() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	std::cout << "  Edge Length: [" << minlen << ", " << maxlen << "]; AVG: " << avelen / polyMesh->numEdges() << std::endl;
}

bool MeshViewerWidget::SaveMesh(const std::string & filename)
{
	return acamcad::polymesh::writeMesh(filename, polyMesh);
}

bool MeshViewerWidget::ScreenShot()
{
	update();
	QString filename = strMeshPath + "/" + QDateTime::currentDateTime().toString("yyyyMMddHHmmsszzz") + QString(".png");
	QImage image = grabFramebuffer();
	image.save(filename);
	std::cout << "Save screen shot to " << filename.toStdString() << std::endl;
	return true;
}

void MeshViewerWidget::SetDrawBoundingBox(bool b)
{
	isDrawBoundingBox = b;
	update();
}
void MeshViewerWidget::SetDrawBoundary(bool b)
{
	isDrawBoundary = b;
	update();
}
void MeshViewerWidget::EnableLighting(bool b)
{
	isEnableLighting = b;
	update();
}
void MeshViewerWidget::EnableDoubleSide(bool b)
{
	isTwoSideLighting = b;
	update();
}

void MeshViewerWidget::ResetView(void)
{
	ResetModelviewMatrix();
	ViewCenter();
	update();
}

void MeshViewerWidget::ViewCenter(void)
{
	if (polyMesh->numVertices()!=0)
	{
		UpdateMesh();
	}
	update();
}

void MeshViewerWidget::CopyRotation(void)
{
	CopyModelViewMatrix();
}

void MeshViewerWidget::LoadRotation(void)
{
	LoadCopyModelViewMatrix();
	update();
}

void MeshViewerWidget::PerformTutte(bool boundaryShape, bool weightStyle)
{
	isBoundaryCircle = boundaryShape;
	isFloaterWeight = weightStyle;
	TutteParam();
}
void MeshViewerWidget::FixBoundarySquare(std::vector<MPoint3>& boundaryCoord) const
{
	int boundaryLength = polyMesh->boundaryVertices().size();
	int cnt = boundaryLength >> 2;
	for (int i = 0; i < boundaryLength; i++) // fix the boundary to a circle
	{
		if (i < cnt)
			boundaryCoord.push_back({ i * 1.0 / cnt, 0.0, 0.0 });
		else if (i < (cnt << 1))
			boundaryCoord.push_back({ 1.0, (i - cnt) * 1.0 / cnt, 0.0 });
		else if (i < (cnt * 3))
			boundaryCoord.push_back({ 1.0 - (i - (cnt << 1)) * 1.0 / cnt, 1.0, 0.0 });
		else boundaryCoord.push_back({ 0.0, 1.0 - (i - (cnt * 3.0)) / cnt, 0.0 });
	}
}
void MeshViewerWidget::TutteParam(void)
{

	if (polyMesh == nullptr) return;
	auto vertices = polyMesh->vertices();
	std::vector<VertPtr> internalVertices;
	for (auto it = polyMesh->vertices_begin(); it != polyMesh->vertices_end(); it++)
		if (!polyMesh->isBoundary(*it)) internalVertices.push_back(*it);
	std::map<VertPtr, int> mp;
	std::map<VertPtr, int> bd_mp;
	for (int i = 0; i < internalVertices.size(); i++)
		mp[internalVertices[i]] = i;
	auto boundaryVertices = polyMesh->boundaryVertices();
	std::vector<MPoint3> boundaryCoord;
	if (isBoundaryCircle)FixBoundaryCircle(boundaryCoord);
	else FixBoundarySquare(boundaryCoord);
	for (int i = 0; i < boundaryVertices.size(); i++)
		bd_mp[boundaryVertices[i]] = i;
	Matrix<Real> A(internalVertices.size(), internalVertices.size());
	Vec b[2] = {Vec(internalVertices.size()), Vec(internalVertices.size())};
	for (int i = 0; i < internalVertices.size(); i++)
	{
		auto lambda = isFloaterWeight ? CalcFloaterWeight(internalVertices[i])
										: CalcUniformWeight(internalVertices[i]);
		auto adj = polyMesh->vertAdjacentVertices(internalVertices[i]);
		for (int k = 0; k < adj.size(); k++)
		{
			if (!polyMesh->isBoundary(adj[k]))
				A[i][mp[adj[k]]] = -lambda[k];
			else
			{
				b[0][i] += lambda[k] * boundaryCoord[bd_mp[adj[k]]].x();
				b[1][i] += lambda[k] * boundaryCoord[bd_mp[adj[k]]].y();
			}
		}
		A[i][i] = 1.0;
	}
	for (auto vert : boundaryVertices)
		vert->setPosition(boundaryCoord[bd_mp[vert]]);
	Matrix<Real> inv(internalVertices.size(), internalVertices.size());
	inverse(A, inv);
	Vec u = inv * b[0], v = inv * b[1];
	for (int i = 0; i < internalVertices.size(); i++)
		internalVertices[i]->setPosition(u[i], v[i], 0.0);
	update();
}

void MeshViewerWidget::PrintMeshInfo(void)
{
	std::cout << "Mesh Info:\n";
	std::cout << "  [V, E, F] = [" << polyMesh->numVertices() << ", " << polyMesh->numEdges() << ", " << polyMesh->numPolygons() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	
}

void MeshViewerWidget::DrawScene(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&projectionmatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&modelviewmatrix[0]);
	//DrawAxis();
	if (isDrawBoundingBox) DrawBoundingBox();
	if (isDrawBoundary) DrawBoundary();
	if (isEnableLighting) glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);
	DrawSceneMesh();
	if (isEnableLighting) glDisable(GL_LIGHTING);
}

void MeshViewerWidget::DrawSceneMesh(void)
{
	if (polyMesh->numVertices() == 0) { return; }
	SetMaterial();
	switch (drawmode)
	{
	case POINTS:
		DrawPoints();
		break;
	case WIREFRAME:
		DrawWireframe();
		break;
	case HIDDENLINES:
		DrawHiddenLines();
		break;
	case FLATLINES:
		DrawFlatLines();
		break;
	case FLAT:
		glColor3d(0.8, 0.8, 0.8);
		DrawFlat();
		break;
	default:
		break;
	}
}

void MeshViewerWidget::FixBoundaryCircle(std::vector<MPoint3>& boundaryCoord) const
{
	int boundaryLength = polyMesh->boundaryVertices().size();
	for (int i = 0; i < boundaryLength; i++) // fix the boundary to a circle
	{
		boundaryCoord.push_back({cos(i * 2.0 * M_PI / boundaryLength),
			sin(i * 2.0 * M_PI / boundaryLength), 0});
	}
}

std::vector<Real> MeshViewerWidget::CalcUniformWeight(VertPtr vertex) const
{
	auto adj = polyMesh->vertAdjacentVertices(vertex);
	std::vector<Real> lambda;
	for (int i = 0; i < adj.size(); i++)
		lambda.push_back(1.0 / adj.size());
	return lambda;
}

void MeshViewerWidget::PolarMap(std::vector<VertPtr>& adj, VertPtr p) const
{
	Real theta = 0.0;
	for (auto vert : adj)
	{
		MVector3 dif(vert->position(), p->position());
		vert->setPosition(dif.x(), dif.y(), dif.z());
	}
	p->setPosition(0.0, 0.0, 0.0);
	std::vector<Real> Ang(adj.size() - 1);
	for (int i = 1; i < adj.size(); i++)
	{
		MVector3 dif(adj[i]->position());
		Ang[i - 1] = vectorAngle(MVector3(adj[i - 1]->position()), dif);
		theta += Ang[i - 1];
	}
	theta += vectorAngle(MVector3(adj[adj.size() - 1]->position()), MVector3(adj[0]->position()));
	Real dist = MVector3(adj[0]->position()).norm();
	adj[0]->setPosition(dist, 0.0, 0.0);
	MVector3 last = {1.0, 0.0, 0.0};//get a unit vector
	for (int i = 1; i < adj.size(); i++)
	{
		MVector3 dif(adj[i]->position());
		dist = dif.norm();
		Real ang = 2 * M_PI * Ang[i - 1] / theta;
		MVector3 dir(last.x() * cos(ang) - last.y() * sin(ang), 
					last.x() * sin(ang) + last.y() * cos(ang), 0.0);
		adj[i]->setPosition(dir.x() * dist, dir.y() * dist, 0.0);
		last = dir;
	}
}


std::vector<Real> MeshViewerWidget::CalcFloaterWeight(VertPtr vertex) const
{
	auto adj = polyMesh->vertAdjacentVertices(vertex);
	std::vector<MPoint3> coords;
	for (auto vert : adj)
		coords.push_back(vert->position());
	coords.push_back(vertex->position());
	PolarMap(adj, vertex);
	std::vector<int> r(adj.size());//r[i] : adj[i] intersects with segment r[i]-r[i+1]
	int ptr = 1, last_sgn = sign(cross(MVector3(adj[1]->position()), 
									MVector3(adj[0]->position())).z());
	for (int k = 2; k < adj.size(); k++)
	{
		int sgn = sign(cross(MVector3(adj[k]->position()), 
						MVector3(adj[0]->position())).z());
		if (last_sgn != sgn)
		{
			ptr = k - 1;
			break;
		}
		else last_sgn = sgn;
	}//now we have r[0] = ptr
	r[0] = ptr;
	int nxt = ptr == adj.size() - 1 ? 0 : ptr + 1;
	for (int k = 1; k < adj.size(); k++)
	{
		while (k == ptr || nxt == ptr)
		{
			ptr = nxt;
			nxt++;
			if (nxt == adj.size()) nxt = 0;
		}
		MVector3 p(adj[k]->position());
		last_sgn = sign(cross(MVector3(adj[ptr]->position()), p).z());
		int sgn = sign(cross(MVector3(adj[nxt]->position()), p).z());
		while (sgn == last_sgn || k == ptr || nxt == k)
		{
			ptr = nxt;
			nxt++;
			if (nxt == adj.size()) nxt = 0;
			last_sgn = sgn;
			sgn = sign(cross(MVector3(adj[nxt]->position()), p).z());
		}
		r[k] = ptr;
	} 
	std::vector<Real> lambda(adj.size());
	auto calcArea = [](const MPoint3& pA, const MPoint3& pB, const MPoint3& pC) -> Real
	{
		return cross(MVector3(pA, pB), MVector3(pA, pC)).norm();
	};
	for (int k = 0; k < adj.size(); k++)
	{
		int nxt = r[k] == adj.size() - 1 ? 0 : r[k] + 1;
		Real area = calcArea(adj[k]->position(), adj[r[k]]->position(), adj[nxt]->position());
		Real deltaA = calcArea({0.0, 0.0, 0.0}, adj[r[k]]->position(), adj[nxt]->position());
		Real deltaB = calcArea({0.0, 0.0, 0.0}, adj[k]->position(), adj[nxt]->position());
		Real deltaC = calcArea({0.0, 0.0, 0.0}, adj[k]->position(), adj[r[k]]->position());
		lambda[k] += deltaA / area;
		lambda[r[k]] += deltaB / area;
		lambda[nxt] += deltaC / area;
	}
	for (Real& x : lambda)
		x /= adj.size();
	for (int k = 0; k < adj.size(); k++)//restore the boundary coordinates
		adj[k]->setPosition(coords[k]);
	vertex->setPosition(coords[adj.size()]);
	return lambda;
}

void MeshViewerWidget::DrawPoints(void) const
{
	glColor3d(1.0, 0.5, 0.5);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : polyMesh->vertices()) {
		glNormal3dv(vh->normal().data());
		glVertex3dv(vh->position().data());
	}
	glEnd();
}

void MeshViewerWidget::DrawWireframe(void) const
{
	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	for (const auto& eh : polyMesh->edges()) {
		auto heh = eh->halfEdge();
		auto v0 = heh->fromVertex();
		auto v1 = heh->toVertex();
		glNormal3dv(v0->normal().data());
		glVertex3dv(v0->position().data());
		glNormal3dv(v1->normal().data());
		glVertex3dv(v1->position().data());
	}
	glEnd();
}

void MeshViewerWidget::DrawHiddenLines() const
{
	glLineWidth(1.0);
	float backcolor[4];
	glGetFloatv(GL_COLOR_CLEAR_VALUE, backcolor);
	glColor4fv(backcolor);
	glDepthRange(0.01, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawFlat();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawFlat();
	}
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3d(.3, .3, .3);
	DrawFlat();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::DrawFlatLines(void) const
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.5f, 2.0f);
	glShadeModel(GL_FLAT);
	//glColor3d(0.8, 0.8, 0.8);
	glColor3d(1.0, 1.0, 1.0);
	DrawFlat();
	glDisable(GL_POLYGON_OFFSET_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawWireframe();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawWireframe();
	}
}

void MeshViewerWidget::DrawFlat(void) const
{
	glBegin(GL_TRIANGLES);
	for (const auto& fh : polyMesh->polyfaces())
	{
		glNormal3dv(fh->normal().data());
		for (const auto& fvh :polyMesh->polygonVertices(fh))
		{
			glVertex3dv(fvh->position().data());
		}
	}
	glEnd();
}


void MeshViewerWidget::DrawBoundingBox(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(.3, .7, .3);
	glBegin(GL_LINES);
	for (const auto& i : { 0, 1 })
	{
		for (const auto& j : { 0, 1 })
		{
			for (const auto& k : { 0, 1 })
			{
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(~i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], ~j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], ~k ? ptMin[2] : ptMax[2]);
			}
		}
	}
	glEnd();
	glLineWidth(linewidth);
}

void MeshViewerWidget::DrawBoundary(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(0.1, 0.1, 0.1);
	glBegin(GL_LINES);

	for (const auto& eh : polyMesh->edges()) {
		if (polyMesh->isBoundary(eh)) {
			auto heh = eh->halfEdge();
			auto v0 = heh->fromVertex();
			auto v1 = heh->toVertex();
			glNormal3dv(v0->normal().data());
			glVertex3dv(v0->position().data());
			glNormal3dv(v1->normal().data());
			glVertex3dv(v1->position().data());
		}
	}
	
	glEnd();
	glLineWidth(linewidth);
}
