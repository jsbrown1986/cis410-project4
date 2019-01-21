/*=========================================================================

			Project 4 for CIS 410 (W18)
			EvaluateVectorFieldAtLocation, AdvectWithEulerStep,
			CalculateArcLength, and Interpolation functions 
			implemented by Jacob Brown 2/1/2018

===========================================================================*/
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
	return idx[1] * dims[0] + idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// Interpolation Function
float Interpolation(float value, float a, float b, float fa, float fb) {
	return (b - a) != 0 ? fa + ((value - a) / (b - a)) * (fb - fa) : 0;
}

// ****************************************************************************
//  Function: EvaluateVectorFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

void EvaluateVectorFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F, float *rv)
{
	int bottomLeftPointGet[2];
	int bottomRightPointGet[2];
	int topLeftPointGet[2];
	int topRightPointGet[2];

	// Checks to see if values are outside arrays before searching
	if (pt[0] < X[0] || pt[0] >= X[dims[0] - 1] || pt[1] < Y[0] || pt[1] >= Y[dims[1] - 1]) {
		rv[0] = 0;
		rv[1] = 0;
		return;
	}

	// Run through X array, seeing if its between two values
	for (int i = 0; i < dims[0] - 1; i++) {
		if (X[i] <= pt[0] && X[i + 1] > pt[0]) {
			bottomLeftPointGet[0] = i;
			break;
		}
	}

	// Run through Y-array, seeing if its between two values
	for (int i = 0; i < dims[1] - 1; i++) {
		if (Y[i] <= pt[1] && Y[i + 1] > pt[1]) {
			bottomLeftPointGet[1] = i;
			break;
		}
	}

	// Uses logical point index to find other indices of the cell
	bottomRightPointGet[0] = bottomLeftPointGet[0] + 1;
	bottomRightPointGet[1] = bottomLeftPointGet[1];

	topLeftPointGet[0] = bottomLeftPointGet[0];
	topLeftPointGet[1] = bottomLeftPointGet[1] + 1;

	topRightPointGet[0] = bottomLeftPointGet[0] + 1;
	topRightPointGet[1] = bottomLeftPointGet[1] + 1;

	// Uses indices to find all EIGHT count 'em EIGHT F-values of the cell
	float bottom_left_f_x = F[2 * GetPointIndex(bottomLeftPointGet, dims)];
	float bottom_left_f_y = F[2 * GetPointIndex(bottomLeftPointGet, dims) + 1];

	float bottom_right_f_x = F[2 * GetPointIndex(bottomRightPointGet, dims)];
	float bottom_right_f_y = F[2 * GetPointIndex(bottomRightPointGet, dims) + 1];

	float top_left_f_x = F[2 * GetPointIndex(topLeftPointGet, dims)];
	float top_left_f_y = F[2 * GetPointIndex(topLeftPointGet, dims) + 1];

	float top_right_f_x = F[2 * GetPointIndex(topRightPointGet, dims)];
	float top_right_f_y = F[2 * GetPointIndex(topRightPointGet, dims) + 1];

	// Performs interpolation to find the x-component of velocity
	float bottom_interp = Interpolation(pt[0], X[bottomLeftPointGet[0]], X[bottomRightPointGet[0]], bottom_left_f_x, bottom_right_f_x);
	float top_interp = Interpolation(pt[0], X[topLeftPointGet[0]], X[topRightPointGet[0]], top_left_f_x, top_right_f_x);

	// Performs interpolation to find the y-component of velocity
	float left_interp = Interpolation(pt[1], Y[bottomLeftPointGet[1]], Y[topLeftPointGet[1]], bottom_left_f_y, top_left_f_y);
	float right_interp = Interpolation(pt[1], Y[bottomRightPointGet[1]], Y[topRightPointGet[1]], bottom_right_f_y, top_right_f_y);

	// Returns x, y velocity
	rv[0] = Interpolation(pt[1], Y[bottomLeftPointGet[1]], Y[topLeftPointGet[1]], bottom_interp, top_interp);
	rv[1] = Interpolation(pt[0], X[bottomLeftPointGet[0]], X[bottomRightPointGet[0]], left_interp, right_interp);
}

// ****************************************************************************
//  Function: AdvectWithEulerStep
//
//  Arguments:
//     pt: the seed location (two-dimensions)
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//     h: The size of the Euler step
//     nsteps: The number of Euler steps to take
//     output_locations (output): An array of size 2*(nsteps+1).  It's first entry
//        should be the seed location.  The second entry should be the result
//        of the first advection step.  The final entry should be the result
//        of the final advection step.
//     speeds (output): An array of size (nsteps+1).  It's first entry should be the
//        speed at the seed location.  It's final entry should be the speed
//        at the location of the result of the final advection step.
//        Recall that speed is the magnitude of the velocity.
//
// ****************************************************************************

void
AdvectWithEulerStep(const float *pt, const int *dims, const float *X, 
                    const float *Y, const float *F, 
                    float h, int nsteps, float *output_locations, float *speeds)
{

	float pointGive[2], speedGet[2];
	float vector1 = 0.0f;
	float vector2 = 0.0f;

	// Gets the first seed location
	output_locations[0] = pt[0];
	output_locations[1] = pt[1];

	pointGive[0] = pt[0];
	pointGive[1] = pt[1];

	// Advects from current step to new step
	for (int i = 0; i <= nsteps; i++) {

		// Gets speed at current location
		EvaluateVectorFieldAtLocation(pointGive, dims, X, Y, F, speedGet);

		// Calculates new location (if applicable)
		if (i < nsteps) {
			output_locations[(2 * i) + 2] = pointGive[0] + (h * speedGet[0]);
			output_locations[(2 * i) + 3] = pointGive[1] + (h * speedGet[1]);

			// Gets this new location for next evaluation of speed
			pointGive[0] = output_locations[(2 * i) + 2];
			pointGive[1] = output_locations[(2 * i) + 3];
		}

		// Finds magnitude of current velocities via Pythagorean theorem
		vector1 = pow(speedGet[0], 2);
		vector2 = pow(speedGet[1], 2);

		speeds[i] = sqrt(vector1 + vector2);
	}
}

// ****************************************************************************
//  Function: CalculateArcLength
//
//  Arguments:
//     locations: an array of 2D locations.
//     nlocations: the number of locations in the array "locations".
//
//  Returns: the arc length, meaning the distance between each successive
//           pair of points
//
// ****************************************************************************

float
CalculateArcLength(const float *output_locations, int nlocations)
{

	float arc_length = 0.0f;
	float x1, x2, y1, y2, x_sqr, y_sqr;

	// Calculates the arc length between every two locations
	for (int i = 0; i < nlocations - 1; i++) {
		x1 = output_locations[(2 * i) + 0];
		y1 = output_locations[(2 * i) + 1];

		x2 = output_locations[(2 * i) + 2];
		y2 = output_locations[(2 * i) + 3];

		x_sqr = pow(x2 - x1, 2);
		y_sqr = pow(y2 - y1, 2);

		arc_length += sqrt(x_sqr + y_sqr);
	}

	return arc_length;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// VTK files are only 3D, so the vector data is all of the form (X,Y,0).
// Remove the 0's since it is counter-intuitive for students who are 
// thinking of this as 2D data.
float *
Convert3DVectorDataTo2DVectorData(const int *dims, const float *F)
{
    float *rv = new float[dims[0]*dims[1]*2];
    int index3D = 0;
    int index2D = 0;
    for (int i = 0 ; i < dims[0] ; i++)
       for (int j = 0 ; j < dims[1] ; j++)
       {
           rv[index2D]   = F[index3D];
           rv[index2D+1] = F[index3D+1];
           index2D += 2;
           index3D += 3;
       }

    return rv;
}

vtkPolyData *
CreateVTKPolyData(int nseeds, int nsteps, float **output_locations, float **speeds)
{
    int numPoints = nseeds*(nsteps+1);
    vtkPoints *pts = vtkPoints::New();
    pts->SetNumberOfPoints(numPoints);
    vtkFloatArray *var = vtkFloatArray::New();
    var->SetName("speed");
    var->SetNumberOfTuples(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nseeds ; i++)
    {
        for (int j = 0 ; j < nsteps+1 ; j++)
        {
            double pt[3];
            pt[0] = output_locations[i][2*j];
            pt[1] = output_locations[i][2*j+1];
            pt[2] = 0.;
            pts->SetPoint(ptIdx, pt);
            var->SetTuple1(ptIdx, speeds[i][j]);
            if (j > 0)
            {
                vtkIdType ids[2] = { ptIdx-1, ptIdx };
                lines->InsertNextCell(2, ids);
            }
            ptIdx++;
        }
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(pts);
    pd->GetPointData()->AddArray(var);
    pd->GetPointData()->SetActiveScalars("speed");
    pd->SetLines(lines);
    lines->Delete();
    var->Delete();
    pts->Delete();

    return pd;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj4_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F_3D = (float *) rgrid->GetPointData()->GetVectors()->GetVoidPointer(0);
    float *F = Convert3DVectorDataTo2DVectorData(dims, F_3D);
    
    
    const int npts = 10;
    float pt[npts][3] =
         {
            {10.1119, 1.22062, 0},
            {8.62376, 13.3839, 0},
            {1.55026, 1.26123, 0},
            {6.9736, 0.653565, 0},
            {2, 2.74117, 0},
            {8.93699, 10.4111, 0},
            {6.08791, -0.533753, 0},
            {10.0543, 1.38024, 0},
            {3.84128, -0.768977, 0},
            {6.66757, 6.0259, 0},
         };


    for (i = 0 ; i < npts ; i++)
    {
       float vec[2];
       EvaluateVectorFieldAtLocation(pt[i], dims, X, Y, F, vec);
       cerr << "Velocity at (" << pt[i][0] <<", "<<pt[i][1] << ") is (" << vec[0] << ", " << vec[1] << ")" << endl;
    }

    float h = 0.01;
    const int nsteps = 5000;
    float **output_locations = new float*[2*(npts+1)];
    float **speeds = new float*[npts+1];
    for (i = 0 ; i < npts ; i++)
    {
       output_locations[i] = new float[(nsteps+1)*2];
       speeds[i] = new float[nsteps];
       AdvectWithEulerStep(pt[i], dims, X, Y, F, h, nsteps, output_locations[i], speeds[i]);
       float length = CalculateArcLength(output_locations[i], nsteps+1);
       cerr << "Arc length for (" << pt[i][0] << ", " << pt[i][1] << ") is " << length << endl;
    }

    vtkPolyData *pd = CreateVTKPolyData(npts, nsteps, output_locations, speeds);

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInput(pd);
    writer->Write();
*/

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);
  
    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);
  
    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();
  
    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);
  
    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);
  
    ren1->GetActiveCamera()->SetFocalPoint(5,5,0);
    ren1->GetActiveCamera()->SetPosition(5,5,30);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    delete [] F;
    pd->Delete();
}
