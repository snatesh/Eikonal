#include <vtkImageViewer2.h>
#include <vtkImageReader2.h>
#include <vtkJPEGReader.h>
#include <vtkPNGReader.h>
#include <vtkTIFFReader.h>
#include <vtkNamedColors.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageSSIM.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <cmath>


vtkSmartPointer<vtkImageData> readImage ( const std::string& pref,
                                          const std::string& ext )
{
  std::string fname = pref + ext;
  vtkSmartPointer<vtkImageReader2> reader;
  if (ext == ".jpg" || ext == ".jpeg")
  {
    vtkSmartPointer<vtkJPEGReader> jpgreader = 
      vtkSmartPointer<vtkJPEGReader>::New();
    reader = jpgreader; 
  }
  else if (ext == ".png")
  {
    vtkSmartPointer<vtkPNGReader> pngreader = 
      vtkSmartPointer<vtkPNGReader>::New();
    reader = pngreader;
  }
  else if (ext == ".tiff")
  {
    vtkSmartPointer<vtkTIFFReader> tiffreader = 
      vtkSmartPointer<vtkTIFFReader>::New();
    reader = tiffreader;
  }
  reader->SetFileName(fname.c_str());
  reader->Update();
  return reader->GetOutput();
}


int main(int argc, char* argv[])
{

  vtkSmartPointer<vtkImageSSIM> ssim = vtkSmartPointer<vtkImageSSIM>::New();
  ssim->SetInputToRGB();
  //ssim->SetInputData(readImage("test",".png")); 
  ssim->SetInputData(readImage("../matlab/raw_images/00046_00_30s",".jpg"));  

  //ssim->SetInputData(readImage("../matlab/raw_images/00049_01_10s",".jpg"));  
  //ssim->SetImageData(readImage("../matlab/raw_images/00049_01_10s",".png"));  
  ssim->SetImageData(readImage("../matlab/raw_images/00046_00_30s",".png"));  
  ssim->Update();
  
  vtkSmartPointer<vtkImageData> ssim_output =  ssim->GetOutput();
  double tup[3], ave = 0;
  int numtuples = ssim_output->GetPointData()->GetNumberOfTuples();
  vtkSmartPointer<vtkDataArray> scalars = ssim_output->GetPointData()->GetScalars();
  for (unsigned int i = 0; i < numtuples; ++i)
  {
    scalars->GetTuple(i, tup);
    ave += (tup[0] + tup[1] + tup[2]) / 3;
  }
  ave /= numtuples; 
 
  std::cout << "ssim average: " << ave << std::endl;
 
  //vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  //writer->SetFileName("ssim.vti");
  //writer->SetInputData(ssim_output);
  //writer->SetDataModeToBinary();
  //writer->Write();

  return 0;

}
