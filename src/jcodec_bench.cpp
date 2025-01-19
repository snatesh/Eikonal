#include <jcodec.hh>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <set>
#include <map>
#include <vector>


namespace fs = std::filesystem;

void getImgQuality  ( const std::string& f1, 
                      const std::string& f2,
                      double& ssim, double& psnr )
{
  std::string command = "./ssim.py " + f1 + " " + f2;
  FILE* fp = popen(command.c_str(), "r");
  char buf[256]; int line = 0;
  while (fgets(buf, sizeof(buf), fp) != NULL)
  {
    if (buf[strlen(buf) - 1] == '\n')
    {
      buf[strlen(buf) - 1] == '\0';
    } 
    if (line == 0) { ssim = atof(buf); }
    if (line == 1) { psnr = atof(buf); }
    line += 1;
  }    
  pclose(fp);
}


int main(int argc, char* argv[])
{
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] << " img_dir ext bench_out" << std::endl;
    std::cerr << "Example: ./jcodec_bench ./rgb8bit ppm rgb8bit_bench" << std::endl;
    std::cerr << "         ./jcodec_bench ./DIV2K_train_HR png DIV2K_bench" << std::endl;
    return EXIT_FAILURE; 
  }
  
  
  std::string path = argv[1];
  std::string ext = argv[2];
  std::string bpref = argv[3];    
  std::set<std::string> images;
  std::map<std::string, std::vector<int>> sizes;  

  for (const auto & entry : fs::directory_iterator(path)) 
  {
    if (!fs::is_directory(entry)) 
    {
      std::string fWithExt(entry.path().filename());
      std::filesystem::path p(fWithExt);
      std::string fNoExt = p.stem().string();
      images.insert(fNoExt);
      
      std::filesystem::path ppng(path+"/"+fNoExt+"."+ext); 
      std::filesystem::path pjpg(path+"/"+fNoExt+".jpg");
      sizes[fNoExt].push_back(std::filesystem::file_size(ppng)); 
      sizes[fNoExt].push_back(std::filesystem::file_size(pjpg)); 
    }
  }
  
  bool useMultiChannel = false;
  unsigned int Nsamps = 20;
  unsigned int Nruns = 5;
  unsigned int Norders = 21;
  unsigned int nSamps[Nsamps];
  unsigned int nRunss[Nruns];
  unsigned int orders[Norders];
 
  for (unsigned int i = 0; i < Nsamps;  ++i)  { nSamps[i] = i + 10; }
  for (unsigned int i = 0; i < Nruns;   ++i)  { nRunss[i] = i + 3;  }
  for (unsigned int i = 0; i < Norders; ++i)  { orders[i] = 8 + i; } 


  double totalBytes[Norders], ssimaves[Norders], psnrs[Norders];
  double psnrs_jpg, ssimaves_jpg;
  

  unsigned int nSamp, nRuns, order;


  auto it = images.begin();
  while (it != images.end())
  {
    bool goimg = true;
    for (unsigned int iSamp = 0; iSamp < Nsamps; ++iSamp)
    {
      for (unsigned int iRun = 0; iRun < Nruns; ++iRun)
      {
        if (goimg)
        {
          nSamp = nSamps[iSamp];
          nRuns = nRunss[iRun];
          std::string bfile = bpref + "_" + std::to_string(nSamp) 
                                    + "_" + std::to_string(nRuns) + ".txt";
          std::ofstream benchrun(bfile, std::ios::app);
          if (benchrun.is_open())
          {
            std::string f1 = path + "/" + *it + "." + ext;
            std::string f3 = path + "/" + *it + ".jpg";
            getImgQuality(f1, f3, ssimaves_jpg, psnrs_jpg);
            Triangulator* T = jcompress_triangulate ( f1.c_str(), nSamp, nRuns,
                                                      useMultiChannel, false );
            for (unsigned int iOrder = 0; iOrder < Norders; ++iOrder)
            {
              unsigned int order = orders[iOrder];
              std::string f2 = *it + "_" + std::to_string(nSamp) + "_" 
                                         + std::to_string(nRuns) + "_" 
                                         + std::to_string(order) + "_deco" + "." + ext;
              std::string vtpfile = *it + "_" + std::to_string(nSamp) + "_" 
                                              + std::to_string(nRuns) + "_" 
                                              + std::to_string(order) + ".vtp";

              // compress with current order
              totalBytes[iOrder] = jcompress  ( T, f1.c_str(),
                                                nSamp, nRuns, order,
                                                useMultiChannel, false );
              // if we compressed to a larger size than jpg, set breakflag
              if (totalBytes[iOrder] > sizes[*it][1]) 
              {
                if (iOrder == 0)
                {
                  goimg = false;
                } 
                break;
              }
              std::cout << "jacobi size: " << totalBytes[iOrder] / 1.e6 << "MB" << std::endl;
              std::cout << "jpg size: " << sizes[*it][1] / 1.e6 << " MB" << std::endl;
              // decompress and write decofile
              jdecompress  ( useMultiChannel, ext.c_str(), vtpfile.c_str() );
              getImgQuality(f1, f2, ssimaves[iOrder], psnrs[iOrder]);
              benchrun << *it << " " << order << " " 
                       << sizes[*it][0] <<  " " << sizes[*it][1] << " "
                       << ssimaves_jpg << " " << psnrs_jpg << " "
                       << totalBytes[iOrder] << " " 
                       << ssimaves[iOrder] << " " << psnrs[iOrder] << std::endl;
              // remove deco and vtp files
              fs::remove(f2); fs::remove(vtpfile);
            }
            delete T;
            benchrun.close();
          }
          else
          {
            std::cerr << "Unable to open file" << std::endl;
          }
        }
      }
    }
    it++;
  } 
  return EXIT_SUCCESS;
}
