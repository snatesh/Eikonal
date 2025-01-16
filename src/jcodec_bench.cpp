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
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " img_dir" << std::endl;
    std::cerr << "Example: ./jcodec_bench ./rgb8bit" << std::endl;
    return EXIT_FAILURE; 
  }
  
  
  //std::string path = "./rgb8bit";
  std::string path = argv[1];
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
      
      std::filesystem::path ppng(path+"/"+fNoExt+".png"); 
      std::filesystem::path pjpg(path+"/"+fNoExt+".jpg");
      sizes[fNoExt].push_back(std::filesystem::file_size(ppng)); 
      sizes[fNoExt].push_back(std::filesystem::file_size(pjpg)); 
    }
  }
  
  bool useMultiChannel = false;
  unsigned int nSamps[25] = { 50, 40, 30, 29, 
                              28, 27, 26, 25, 
                              24, 23, 22, 21, 
                              20, 19, 18, 17, 
                              16, 15, 14, 13, 
                              11, 10, 11, 12, 10 }; 
  unsigned int nRunss[5] = {0, 1, 2, 3, 4};
  unsigned int orders[21] = {5, 6, 7, 8, 9,
                             10, 11, 12, 13, 14,
                             15, 16, 17, 18, 19,
                             20, 21, 22, 23, 24, 25};
  double totalBytes[21], ssimaves[21], psnrs[21];
  double psnrs_jpg, ssimaves_jpg;
  
  std::string bpref = "benchjcodec";    
  unsigned int nSamp, nRuns, order;
  for (unsigned int iSamp = 0; iSamp < 25; ++iSamp)
  {
    nSamp = nSamps[iSamp];
    for (unsigned int iRun = 0; iRun < 5; ++iRun)
    {
      nRuns = nRunss[iRun];
      std::string bfile = bpref + "_" + std::to_string(nSamp) 
                                + "_" + std::to_string(nRuns) + ".txt";
      std::ofstream benchrun(bfile);
  
      if (benchrun.is_open())
      {
        auto it = images.begin();
        while (it != images.end())
        {
          std::string f1 = path + "/" + *it + ".png";
          std::string f3 = path + "/" + *it + ".jpg";
          getImgQuality(f1, f3, ssimaves_jpg, psnrs_jpg);
          Triangulator* T = jcompress_triangulate ( f1.c_str(), nSamp, nRuns,
                                                    useMultiChannel, false );
          for (unsigned int iOrder = 0; iOrder < 21; ++iOrder)
          {
            unsigned int order = orders[iOrder];
            std::string f2 = *it + "_" + std::to_string(order) + "_deco" + ".png";
            std::string vtpfile = *it + "_" + std::to_string(order) + ".vtp";
            // compress with current order
            totalBytes[iOrder] = jcompress ( T, f1.c_str(),
                                        nSamp, nRuns, order,
                                        useMultiChannel, false );
            if (totalBytes[iOrder] / sizes[*it][1] > 1) 
            { 
              break; 
            }
            std::cout << "jacobi size: " << totalBytes[iOrder] / 1.e6 << "MB" << std::endl;
            std::cout << "jpg size: " << sizes[*it][1] / 1.e6 << " MB" << std::endl;
            // decompress and write decofile
            jdecompress  ( useMultiChannel, "png", vtpfile.c_str() );
            getImgQuality(f1, f2, ssimaves[iOrder], psnrs[iOrder]);
  
            benchrun << *it << " " << order << " " 
                     << sizes[*it][0] <<  " " << sizes[*it][1] << " "
                     << ssimaves_jpg << " " << psnrs_jpg << " "
                     << totalBytes[iOrder] << " " 
                     << ssimaves[iOrder] << " " << psnrs[iOrder] << std::endl;
  
          }
          delete T;
          it++;
        }
        benchrun.close();
      }
      else
      {
        std::cerr << "Unable to open file" << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  return EXIT_SUCCESS;
}
