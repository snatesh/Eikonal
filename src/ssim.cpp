#include <vtkImageViewer2.h>
#include <vtkImageReader2.h>
#include <vtkJPEGReader.h>
#include <vtkPNGReader.h>
#include <vtkPNMReader.h>
#include <vtkTIFFReader.h>
#include <vtkNamedColors.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageSSIM.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <filesystem>
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
  else if (ext == ".ppm")
  {
    vtkSmartPointer<vtkPNMReader> ppmreader = 
      vtkSmartPointer<vtkPNMReader>::New();
    reader = ppmreader; 
  }
  reader->SetFileName(fname.c_str());
  reader->Update();
  return reader->GetOutput();
}

int main(int argc, char* argv[])
{
  vtkSmartPointer<vtkImageSSIM> ssim = vtkSmartPointer<vtkImageSSIM>::New();
  ssim->SetInputToRGB();
  
  std::string f1WithExt(argv[1]);
  std::filesystem::path p1(f1WithExt);
  std::string f1 = p1.stem().string();
  std::string ext1 = p1.extension();
  std::string f1WithoutExt;
  if (not p1.parent_path().string().empty())
  {
    f1WithoutExt = p1.parent_path().string() + "/" + f1; 
  }
  else
  {
    f1WithoutExt = f1;
  }

  std::string f2WithExt(argv[2]);
  std::filesystem::path p2(f2WithExt);
  std::string f2 = p2.stem().string();
  std::string ext2 = p2.extension();
  std::string f2WithoutExt;
  if (not p2.parent_path().string().empty())
  {
    f2WithoutExt = p2.parent_path().string() + "/" + f2;
  }
  else
  {
    f2WithoutExt = f2;
  }
  ssim->SetInputData(readImage(f1WithoutExt, ext1));  
  ssim->SetImageData(readImage(f2WithoutExt, ext2));  

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
