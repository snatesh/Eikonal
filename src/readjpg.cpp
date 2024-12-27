#include<jcodec.hh>
#include <vtkImageViewer2.h>
#include <vtkJPEGReader.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkPointData.h>
#include <map>
#include <vector>

void writeVTP(vtkSmartPointer<vtkPolyData> polytri, const char* ofname)
{
  // write mesh
  vtkNew<vtkXMLPolyDataWriter> writer;
  writer->SetFileName(ofname);
  writer->SetInputData(polytri);  
  writer->SetDataModeToBinary();
  writer->Write();
}

int usage(char* argv[])
{
  std::cout << "Usage: " << argv[0]
            << " mode args\n"
            << " if (mode=0), args = filename.jpg nSamp nRuns ctol \n"
            << " if (mode=1), args = channel1.vtp channel2.vtp channel3.vtp\n";
  return EXIT_FAILURE;

}

int main(int argc, char* argv[])
{

  unsigned int mode = static_cast<unsigned int>(atoi(argv[1]));
  if (mode != 0 && mode != 1)
  {
    return usage(argv);
  }
  
  if (mode == 0)
  {
    if (argc != 6) { return usage(argv); }
    // jacobi poly
    double a = 0.5, b = 0.5, c = 0.5;
    unsigned int N = 55;
    unsigned int m = 6;
    unsigned int nthreads = 6;
    double ctol = atof(argv[5]);
    // Read the image
    vtkNew<vtkJPEGReader> jpegReader;
    jpegReader->SetFileName(argv[2]);
    jpegReader->Update();
    vtkSmartPointer<vtkImageData> imagedata = jpegReader->GetOutput(); 

    // create cubic interpolator on image
     vtkSmartPointer<vtkImageInterpolator> interpolator = 
      vtkSmartPointer<vtkImageInterpolator>::New();
    interpolator->Initialize(imagedata);
    interpolator->SetInterpolationModeToCubic();
    // get triangulation settings
    unsigned int nSamp = static_cast<unsigned int>(atoi(argv[3]));
    unsigned int nRuns = static_cast<unsigned int>(atoi(argv[4]));
    Triangulator* T = new Triangulator  ( N, m, a, b, c, nSamp, nRuns, 
                                          "xtri_N55_n9_M91_m12.txt",
                                          "ytri_N55_n9_M91_m12.txt",            
                                          "wtri_N55_n9_M91_m12.txt",
                                           imagedata, interpolator, nthreads);
    T->run();
    Compressor* C = new Compressor(T);
    C->run(ctol);
    // write grids with compressed data
    writeVTP(T->polytri1, "test1.vtp");
    writeVTP(T->polytri2, "test2.vtp");
    writeVTP(T->polytri3, "test3.vtp");
  
    // Visualize original image
    vtkNew<vtkNamedColors> colors;
    vtkNew<vtkImageViewer2> imageViewer;
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

    // clean
    delete T;
    delete C;
  }
  else if (mode == 1)
  {
    if (argc != 5) { return usage(argv); }
    const char *channel1, *channel2, *channel3;
    channel1 = argv[2];
    channel2 = argv[3];
    channel3 = argv[4];
    // test decompressor reading
    Decompressor* D = new Decompressor ( "test1.vtp", "test2.vtp", "test3.vtp" );
    D->run();
    D->writeImage("test.png");
    delete D;
  }
  else
  {
    std::cout << "Usage: " << argv[0]
              << " mode args\n "
              << "if (mode=0), args = filename.jpg nSamp nRuns\n"
              << "if (mode=1), args = channel1.vtp channel2.vtp channel3.vtp\n";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}




/* Random test codes */

  //for (int icell = 0; icell < polytri->GetNumberOfCells(); ++icell)
  //{

  //  polytri->GetCellPoints(icell, ptids);
  //  polytri->GetPoints()->GetPoint(ptids->GetId(0), v1);
  //  polytri->GetPoints()->GetPoint(ptids->GetId(1), v2);
  //  polytri->GetPoints()->GetPoint(ptids->GetId(2), v3);
  //
  //} 



  // testing coef expansions and derivatives
  //double* F = (double*) calloc(N, sizeof(double));
  //double* cF = (double*) calloc(M, sizeof(double));
  //double* cdxF = (double*) calloc(M, sizeof(double));
  //double* cdyF = (double*) calloc(M, sizeof(double));
  //for (unsigned int i = 0; i < N; ++i) { F[i] = std::pow(X[i] + Y[i], 6); }
  //Pm->computeCoeffs(F, X, Y, W, cF);
  //cblas_dgemv( CblasColMajor, CblasNoTrans,  M, M, 1.0, Dx, M, cF, 1, 0.0, cdxF, 1 );
  //cblas_dgemv( CblasColMajor, CblasNoTrans,  M, M, 1.0, Dy, M, cF, 1, 0.0, cdyF, 1 );
  //for (unsigned int i = 0; i < M; ++i) 
  //{
  //  printf("%.17g\n", cF[i]);
  //}
  //printf("\n");
  //for (unsigned int i = 0; i < M; ++i) 
  //{
  //  printf("%.17g\n", cdxF[i]);
  //}
  //printf("\n");
  //for (unsigned int i = 0; i < M; ++i) 
  //{
  //  printf("%.17g\n", cdyF[i]);
  //}
  //printf("\n");
  //free(F); free(cF);
  //free(cdxF); free(cdyF);

  //std::cout << imagedata->GetNumberOfPoints() << std::endl;
  //std::cout << imagedata->GetNumberOfCells() << std::endl;
  //double x[3];
  //imagedata->GetPoint(24160255, x);
  //std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
  //double origin[3];
  //imagedata->GetOrigin(origin);
  //std::cout << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
  //std::cout << imagedata->GetActualMemorySize() * 1024 << std::endl;
  //int dims[3];
  //imagedata->GetDimensions(dims);
  //std::cout << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;



  //double xqtr[3], xqt[3]; 
  //xqtr[2] = 0; xqt[2] = 0;
  //double interpwts[3];

  //for (unsigned int j = 0; j < N; ++j)
  //{
  //  // get quad point in ref
  //  xqtr[0] = X[j]; xqtr[1] = Y[j]; 
  //  // map quad point into each cell
  //  for (int i = 0; i < polytri->GetNumberOfCells(); ++i)
  //  {
  //    int subid;
  //    polytri->GetCell(i)->EvaluateLocation(subid, xqtr, xqt, interpwts);
  //    points->InsertNextPoint(xqt); 
  //  } 
  //}



  //vtkDoubleArray* gnint = vtkDoubleArray::New();
  //gnint->SetName("Integral(GradNorm)");
  //gnint->SetNumberOfComponents(3);
  //gnint->SetNumberOfTuples(polytri->GetNumberOfCells());
  //double sums[3];
  //for (unsigned int i = 0; i < polytri->GetNumberOfCells(); ++i)
  //{
  //  sums[0] = intgn1[i]; sums[1] = intgn2[i]; sums[2] = intgn3[i];
  //  gnint->SetTuple(i, sums);
  //}

  //polytri->GetCellData()->AddArray(gnint); 
  //vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  //points->DeepCopy(polytri->GetPoints());


  //std::cout << "NEW num points " << polytri->GetNumberOfPoints() << std::endl;
  //vtkNew<vtkXMLPolyDataWriter> writer;
  //writer->SetFileName("test.vtp");
  //writer->SetInputData(polytri);  
  ////writer->SetDataModeToBinary();
  //writer->SetDataModeToAscii();
  //writer->Write();


    //std::cout << imagedata->GetPointData()->GetNumberOfTuples() << std::endl;
    //std::cout << imagedata->GetPointData()->GetNumberOfComponents() << std::endl;
    //std::cout << imagedata->GetScalarTypeAsString() << std::endl;
    //std::cout << imagedata->GetNumberOfScalarComponents() << std::endl;
