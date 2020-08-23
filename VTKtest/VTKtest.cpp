#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingContextOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

#include <vtkSmartPointer.h>
#include <vtkParametricFunctionSource.h>
#include <vtkParametricSpline.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkProperty.h>

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkGlyph3DMapper.h>

#include <vtkSphereSource.h>
#include <vtkNamedColors.h>

#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

int main(int, char *[])
{
	vtkSmartPointer<vtkNamedColors> colors =
		vtkSmartPointer<vtkNamedColors>::New();

	// Create three points. We will join (Origin and P0) with a red line and (Origin and P1) with a green line
	double origin[3] = { 0.0, 0.0, 0.0 };
	double p0[3] = { 1.0, 0.0, 0.0 };
	double p1[3] = { 0.0, 1.0, 0.0 };
	double p2[3] = { 0.0, 1.0, 2.0 };
	double p3[3] = { 1.0, 2.0, 3.0 };

	// Create a vtkPoints object and store the points in it
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->InsertNextPoint(origin);
	points->InsertNextPoint(p0);
	points->InsertNextPoint(p1);
	points->InsertNextPoint(p2);
	points->InsertNextPoint(p3);

	vtkSmartPointer<vtkParametricSpline> spline =
		vtkSmartPointer<vtkParametricSpline>::New();
	spline->SetPoints(points);

	vtkSmartPointer<vtkParametricFunctionSource> functionSource =
		vtkSmartPointer<vtkParametricFunctionSource>::New();
	functionSource->SetParametricFunction(spline);
	functionSource->Update();


	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	//transform->RotateWXYZ(double angle, double x, double y, double z);
	transform->RotateWXYZ(110, 1, 1, 0);

	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
		vtkSmartPointer<vtkTransformPolyDataFilter>::New();

	transformFilter->SetTransform(transform);
	transformFilter->SetInputConnection(functionSource->GetOutputPort());
	transformFilter->Update();

	// Setup actor and mapper-------------------------------
	vtkSmartPointer<vtkPolyDataMapper> mapper_orig =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper_orig->SetInputConnection(functionSource->GetOutputPort());

	vtkPolyData* polydata_orig = functionSource->GetOutput();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(transformFilter->GetOutputPort());

	vtkPolyData* polydata = functionSource->GetOutput();

	vtkSmartPointer<vtkDoubleArray> ds = vtkSmartPointer<vtkDoubleArray>::New();
	ds->SetName("Distances");
	ds->SetNumberOfComponents(1);
	ds->SetNumberOfValues(polydata->GetNumberOfPoints() - 1);

	vtkSmartPointer<vtkDoubleArray> kv = vtkSmartPointer<vtkDoubleArray>::New();
	kv->SetName("Curvatures");
	kv->SetNumberOfComponents(1);
	kv->SetNumberOfValues(polydata->GetNumberOfPoints() - 1);

	vtkSmartPointer<vtkDoubleArray> tr = vtkSmartPointer<vtkDoubleArray>::New();
	tr->SetName("Torsions");
	tr->SetNumberOfComponents(1);
	tr->SetNumberOfValues(polydata->GetNumberOfPoints() - 1);

	vtkSmartPointer<vtkDoubleArray> k2d = vtkSmartPointer<vtkDoubleArray>::New();
	k2d->SetName("Curvatures2d");
	k2d->SetNumberOfComponents(1);
	k2d->SetNumberOfValues(polydata->GetNumberOfPoints() - 1);

	//-----------------------------------------
	vtkSmartPointer<vtkDoubleArray> dxi = vtkSmartPointer<vtkDoubleArray>::New();
	dxi->SetName("dxi");
	dxi->SetNumberOfComponents(3);
	dxi->SetNumberOfTuples(polydata->GetNumberOfPoints());

	vtkSmartPointer<vtkDoubleArray> dti = vtkSmartPointer<vtkDoubleArray>::New();
	dti->SetName("dti");
	dti->SetNumberOfComponents(3);
	dti->SetNumberOfTuples(polydata->GetNumberOfPoints());

	vtkSmartPointer<vtkDoubleArray> dni = vtkSmartPointer<vtkDoubleArray>::New();
	dni->SetName("dni");
	dni->SetNumberOfComponents(3);
	dni->SetNumberOfTuples(polydata->GetNumberOfPoints());

	vtkSmartPointer<vtkDoubleArray> dbi = vtkSmartPointer<vtkDoubleArray>::New();
	dbi->SetName("dbi");
	dbi->SetNumberOfComponents(3);
	dbi->SetNumberOfTuples(polydata->GetNumberOfPoints());


	vtkSmartPointer<vtkDoubleArray> ti = vtkSmartPointer<vtkDoubleArray>::New();
	ti->SetName("ti");
	ti->SetNumberOfComponents(3);
	ti->SetNumberOfTuples(polydata->GetNumberOfPoints());

	vtkSmartPointer<vtkDoubleArray> ni = vtkSmartPointer<vtkDoubleArray>::New();
	ni->SetName("ni");
	ni->SetNumberOfComponents(3);
	ni->SetNumberOfTuples(polydata->GetNumberOfPoints());

	vtkSmartPointer<vtkDoubleArray> bi = vtkSmartPointer<vtkDoubleArray>::New();
	bi->SetName("bi");
	bi->SetNumberOfComponents(3);
	bi->SetNumberOfTuples(polydata->GetNumberOfPoints());

	//dxi,ti
	for (vtkIdType i = 1; i < polydata->GetNumberOfPoints() - 1; i++) {
		double pnow[3], pnext[3], plast[3];
		polydata->GetPoint(i, pnow);
		polydata->GetPoint(i + 1, pnext);
		polydata->GetPoint(i - 1, plast);

		if (i == 10) std::cout << pnow[1] << std::endl;

		double squaredDistance1 = vtkMath::Distance2BetweenPoints(pnow, plast);
		double distance1 = sqrt(squaredDistance1);
		double squaredDistance2 = vtkMath::Distance2BetweenPoints(pnext, pnow);
		double distance2 = sqrt(squaredDistance2);

		double tmpSub1[3], tmpSub2[3];
		vtkMath::Subtract(pnow, plast, tmpSub1);
		vtkMath::Subtract(pnext, pnow, tmpSub2);

		vtkMath::MultiplyScalar(tmpSub1, (double)1 / (distance1*2.0));
		vtkMath::MultiplyScalar(tmpSub2, (double)1 / (distance2*2.0));

		double tmpAdd[3];
		vtkMath::Add(tmpSub1, tmpSub2, tmpAdd);
		dxi->SetTuple(i, tmpAdd);

		double tmpNorm;
		tmpNorm = vtkMath::Norm(tmpAdd);

		double tmpT[3];
		dxi->GetTuple(i, tmpT);

		vtkMath::MultiplyScalar(tmpT, 1.0 / tmpNorm);
		ti->SetTuple(i, tmpT);
	}

	//dti,ni,kv,bi
	for (vtkIdType i = 2; i < polydata->GetNumberOfPoints() - 2; i++) {
		double pnow[3], pnext[3], plast[3];
		polydata->GetPoint(i, pnow);
		polydata->GetPoint(i + 1, pnext);
		polydata->GetPoint(i - 1, plast);

		double squaredDistance1 = vtkMath::Distance2BetweenPoints(pnow, plast);
		double distance1 = sqrt(squaredDistance1);
		double squaredDistance2 = vtkMath::Distance2BetweenPoints(pnext, pnow);
		double distance2 = sqrt(squaredDistance2);

		double tnow[3], tnext[3], tlast[3];
		ti->GetTuple(i, tnow);
		ti->GetTuple(i + 1, tnext);
		ti->GetTuple(i - 1, tlast);

		double tmpSub1[3], tmpSub2[3];
		vtkMath::Subtract(tnow, tlast, tmpSub1);
		vtkMath::Subtract(tnext, tnow, tmpSub2);

		vtkMath::MultiplyScalar(tmpSub1, (double)1 / (distance1*2.0));
		vtkMath::MultiplyScalar(tmpSub2, (double)1 / (distance2*2.0));

		double tmpAdd[3];
		vtkMath::Add(tmpSub1, tmpSub2, tmpAdd);
		dti->SetTuple(i, tmpAdd);

		double tmpNorm;
		tmpNorm = vtkMath::Norm(tmpAdd);

		double tmpN[3];
		dti->GetTuple(i, tmpN);

		vtkMath::MultiplyScalar(tmpN, 1.0 / tmpNorm);
		ni->SetTuple(i, tmpN);

		double tmpNorm2;
		tmpNorm2 = vtkMath::Norm(tmpAdd);
		kv->SetValue(i, tmpNorm2);

		

		double tmpni[3], tmpti[3], tmpbi[3];
		ti->GetTuple(i, tmpti);
		ni->GetTuple(i, tmpni);

		vtkMath::Cross(tmpti, tmpni, tmpbi);
		bi->SetTuple(i, tmpbi);
	}

	//dni,dbi
	for (vtkIdType i = 3; i < polydata->GetNumberOfPoints() - 3; i++) {
		double pnow[3], pnext[3], plast[3];
		polydata->GetPoint(i, pnow);
		polydata->GetPoint(i + 1, pnext);
		polydata->GetPoint(i - 1, plast);

		double squaredDistance1 = vtkMath::Distance2BetweenPoints(pnow, plast);
		double distance1 = sqrt(squaredDistance1);
		double squaredDistance2 = vtkMath::Distance2BetweenPoints(pnext, pnow);
		double distance2 = sqrt(squaredDistance2);

		double nnow[3], nnext[3], nlast[3];
		ni->GetTuple(i, nnow);
		ni->GetTuple(i + 1, nnext);
		ni->GetTuple(i - 1, nlast);

		double bnow[3], bnext[3], blast[3];
		bi->GetTuple(i, bnow);
		bi->GetTuple(i + 1, bnext);
		bi->GetTuple(i - 1, blast);

		double tmpSub1[3], tmpSub2[3];
		vtkMath::Subtract(nnow, nlast, tmpSub1);
		vtkMath::Subtract(nnext, nnow, tmpSub2);

		double tmpSub3[3], tmpSub4[3];
		vtkMath::Subtract(bnow, blast, tmpSub3);
		vtkMath::Subtract(bnext, bnow, tmpSub4);

		vtkMath::MultiplyScalar(tmpSub1, (double)1 / (distance1*2.0));
		vtkMath::MultiplyScalar(tmpSub2, (double)1 / (distance2*2.0));

		vtkMath::MultiplyScalar(tmpSub3, (double)1 / (distance1*2.0));
		vtkMath::MultiplyScalar(tmpSub4, (double)1 / (distance2*2.0));

		double tmpAdd[3];
		vtkMath::Add(tmpSub1, tmpSub2, tmpAdd);
		dni->SetTuple(i, tmpAdd);

		double tmpAdd1[3];
		vtkMath::Add(tmpSub3, tmpSub4, tmpAdd1);
		dbi->SetTuple(i, tmpAdd1);

		double tmp1[3], tmp2, tmp3[3], tmp4[3];
		dni->GetTuple(i, tmp1);
		tmp2 = kv->GetValue(i);
		ti->GetTuple(i, tmp3);

		vtkMath::MultiplyScalar(tmp3, tmp2);
		vtkMath::Add(tmp1, tmp3, tmp4);

		double tmptr;
		tmptr = vtkMath::Norm(tmp4);

		tr->SetValue(i, tmptr);
		std::cout << "point" << i << ": " << tmptr << std::endl;
		//std::cout << tmptr << std::endl;tr
	}

	vtkSmartPointer<vtkActor> actor_orig =
		vtkSmartPointer<vtkActor>::New();
	actor_orig->SetMapper(mapper_orig);
	actor_orig->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());
	actor_orig->GetProperty()->SetLineWidth(1.0);

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
	actor->GetProperty()->SetLineWidth(1.0);

	// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> polyData =
		vtkSmartPointer<vtkPolyData>::New();
	polyData->SetPoints(points);


	// Setup render window, renderer, and interactor
	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderer->AddActor(actor_orig);
	renderer->AddActor(actor);
	renderer->SetBackground(colors->GetColor3d("Silver").GetData());

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
		vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); //like paraview

	renderWindowInteractor->SetInteractorStyle(style);

	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}