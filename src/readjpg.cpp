#include <vtkImageViewer2.h>
#include <vtkJPEGReader.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkImageToStructuredGrid.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkStructuredGrid.h"
#include "vtkUniformGrid.h"
#include "vtkPolyData.h"
#include "vtkDelaunay2D.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkSmartPointer.h"

vtkSmartPointer<vtkDelaunay2D> triangulateImage ( int* dims, 
                                                double* origin,
                                                unsigned int nSamp )
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  unsigned int numPoints = dims[0]*dims[1];
  points->SetNumberOfPoints(nSamp * nSamp);
  double xstride =  static_cast<double>((dims[0]-1-origin[0])) / 
                    static_cast<double>((nSamp - 1));
  double ystride =  static_cast<double>((dims[1]-1-origin[0])) / 
                    static_cast<double>((nSamp - 1));
   
  int ptId = 0; double pt[3]; pt[2] = 0;
  for (unsigned int iY = 0; iY < nSamp; ++iY)
  {
    for (unsigned int iX = 0; iX < nSamp; ++iX)
    {
      pt[0] = iX * xstride; pt[1] = iY * ystride;
      points->SetPoint(ptId, pt);
      ptId += 1;  
    }
  }
 
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New(); 
  polydata->SetPoints(points);
  vtkSmartPointer<vtkDelaunay2D> triangulator = vtkSmartPointer<vtkDelaunay2D>::New();
  triangulator->SetInputData(polydata);
  triangulator->Update();
  return triangulator;

}


int main(int argc, char* argv[])
{
  vtkNew<vtkNamedColors> colors;

  // Verify input arguments
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0]
              << " Filename(.jpeg/jpg) e.g. Pileated.jpg " << std::endl;
    return EXIT_FAILURE;
  }

  // Read the image
  vtkNew<vtkJPEGReader> jpegReader;
  jpegReader->SetFileName(argv[1]);
  jpegReader->Update();
  // Visualize
  vtkNew<vtkImageViewer2> imageViewer;
  vtkImageData* imagedata = jpegReader->GetOutput(); 

  std::cout << imagedata->GetNumberOfPoints() << std::endl;
  std::cout << imagedata->GetNumberOfCells() << std::endl;
  double x[3];
  imagedata->GetPoint(24160255, x);
  std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
  double origin[3];
  imagedata->GetOrigin(origin);
  std::cout << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
  std::cout << imagedata->GetActualMemorySize() * 1024 << std::endl;
  int dims[3];
  imagedata->GetDimensions(dims);
  std::cout << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;

  vtkSmartPointer<vtkDelaunay2D> triangulator = triangulateImage(dims, origin, 10);
  vtkSmartPointer<vtkPolyData> polytri = triangulator->GetOutput();
  std::cout << "HERE" << std::endl;
  std::cout << polytri->GetNumberOfPoints() << std::endl;
  std::cout << polytri->GetNumberOfCells() << std::endl;
  double p[3];
  polytri->GetPoint(25,p);
  std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
  vtkNew<vtkXMLPolyDataWriter> writer;
  writer->SetFileName("test.vtp");
  writer->SetInputData(polytri);  
  writer->SetDataModeToBinary();
  writer->Write();
  

  imageViewer->SetInputData(imagedata);
  
  
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  imageViewer->SetupInteractor(renderWindowInteractor);
  imageViewer->Render();
  imageViewer->GetRenderer()->ResetCamera();
  imageViewer->GetRenderer()->SetBackground(
      colors->GetColor3d("DarkSlateGray").GetData());
  imageViewer->GetRenderWindow()->SetWindowName("JPEGReader");
  imageViewer->Render();

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
