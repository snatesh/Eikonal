#include <jcodec.hh>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <set>
#include <map>
#include <vector>

namespace fs = std::filesystem;

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " mode" << std::endl;
    return EXIT_FAILURE; 
  }
  
  int mode = atoi(argv[1]);
  
  std::string path = "./rgb8bit";
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
      
      std::filesystem::path pppm(path+"/"+fNoExt+".ppm"); 
      std::filesystem::path pjpg(path+"/"+fNoExt+".jpg");
      sizes[fNoExt].push_back(std::filesystem::file_size(pppm)); 
      sizes[fNoExt].push_back(std::filesystem::file_size(pjpg)); 
    }
  }

  // generate reference ssims
  if (mode == 0)
  {
    std::ofstream benchref("benchjpg.txt");
    if (benchref.is_open())
    {
      double ssimave; auto it = images.begin();
      while (it != images.end() )
      {
        std::string f1 = "./rgb8bit/" + *it + ".ppm";
        std::string f2 = "./rgb8bit/" + *it + ".jpg";
        ssimave = ssim ( f1.c_str(), f2.c_str() );
        benchref << *it << " " << sizes[*it][0] << " " << sizes[*it][1] << " " << ssimave << std::endl;
        it++;
      }
      benchref.close();
    }
    else
    {
      std::cerr << "Unable to open file" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    bool useMultiChannel = false;
    unsigned int nSamp = 10, nRuns = 4;
    std::ofstream benchrun("benchjcodec.txt");
    if (benchrun.is_open())
    {
      auto it = images.begin();
      while (it != images.end())
      {
        unsigned int order = 15;         
        double totalBytes[1], ssimaves[1];
        unsigned int orders[1];
        std::string f1 = "./rgb8bit/" + *it + ".ppm";
        Triangulator* T = jcompress_triangulate ( f1.c_str(), 
                                                  nSamp, nRuns, order,
                                                  useMultiChannel, false );
        for (unsigned int i = 0; i < 1; ++i)
        {
          // input file name
          std::string f3 = "./rgb8bit/" + *it + ".jpg";
          std::string f2 = *it + "_" + std::to_string(order) + "_deco" + ".ppm";
          // output file names
          std::string vtpfile = *it + "_" + std::to_string(order) + ".vtp";
          // compress with current order
          totalBytes[i] = jcompress ( T, f1.c_str(),
                                      nSamp, nRuns, order,
                                      useMultiChannel, false );
          std::cout << "jpg size: " << sizes[*it][1] / 1.e6 << " MB" << std::endl;
          orders[i] = order;
          // decompress and write decofile
          jdecompress  ( useMultiChannel, "ppm", vtpfile.c_str() );
          ssimaves[i] = ssim ( f1.c_str(), f2.c_str() );
          order += 1;
        }
        benchrun << *it << " " << sizes[*it][0] <<  " " << sizes[*it][1] << " ";
        for (unsigned int i = 0; i < 1; ++i)
        {
          benchrun << totalBytes[i] << " " << ssimaves[i] << " ";
        }
        benchrun << std::endl;
        it++;
        delete T;
      }
      benchrun.close();
    }
    else
    {
      std::cerr << "Unable to open file" << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
