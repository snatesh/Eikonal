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
  
  std::string path = "./raw_images";
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

  // generate reference ssims
  if (mode == 0)
  {
    std::ofstream benchref("benchjpg.txt");
    if (benchref.is_open())
    {
      double ssim; auto it = images.begin();
      while (it != images.end() )
      {
        std::string pref = "./raw_images/" + *it;
        ssim = ssim_jpg ( pref.c_str() );
        benchref << *it << " " << sizes[*it][0] << " " << sizes[*it][1] << " " << ssim << std::endl;
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
    unsigned int nSamp = 10, nRuns = 3;
    std::ofstream benchrun("benchjcodec.txt");
    if (benchrun.is_open())
    {
      auto it = images.begin();
      while (it != images.end())
      {
        unsigned int order = 5;         
        double totalBytes[16], ssims[16];
        unsigned int orders[16];
        for (unsigned int i = 0; i < 16; ++i)
        {
          // input file name
          std::string pref = "./raw_images/" + *it;
          std::string pngfile = "./raw_images/" + *it + ".png";
          std::string decopref = *it + "_" + std::to_string(order) + "_deco";
          // output file names
          std::string vtpfile = *it + "_" + std::to_string(order) + ".vtp";
          // decompress with current order
          totalBytes[i] = jcompress ( pngfile.c_str(),
                                      nSamp, nRuns, order,
                                      useMultiChannel, false );
          orders[i] = order;
          // decompress and write decofile
          jdecompress  ( useMultiChannel, 
                         vtpfile.c_str() );
          ssims[i] = ssim_png( pref.c_str(), decopref.c_str() ); 
          order += 1;
        }
        benchrun << *it << " " << sizes[*it][0] <<  " " << sizes[*it][1] << " ";
        for (unsigned int i = 0; i < 16; ++i)
        {
          benchrun << totalBytes[i] << " " << ssims[i] << " ";
        }
        benchrun << std::endl;
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
  return EXIT_SUCCESS;
}
