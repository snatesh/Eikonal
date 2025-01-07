#include <jcodec.hh>

int usage(char* argv[])
{
  std::cout << "Usage: " << argv[0]
            << " mode multichannel args\n"
            << " if (mode=0), args = filename.jpg nSamp nRuns order \n"
            << " if (mode=1) and (multichannel=1), args = f1.vtp f2.vtp f3.vtp out_fmt)\n"
            << " if (mode=1) and (multichannel=0), args = f.vtp out_fmt\n"
            << " where out_fmt can be png, jpg, ppm, tiff\n";
  return EXIT_FAILURE;

}

int main(int argc, char* argv[])
{
  if (argc < 4) { return usage(argv); }
  
  unsigned int mode = static_cast<unsigned int>(atoi(argv[1]));
  
  if (mode != 0 && mode != 1) { return usage(argv); }
  
  bool useMultiChannel = static_cast<bool>(atoi(argv[2]));
  
  if (useMultiChannel != 0 && useMultiChannel != 1) { return usage(argv); }

  if (mode == 0)
  {
    if (argc != 7) { return usage(argv); }
    unsigned int order = static_cast<unsigned int>(atoi(argv[6]));
    // image file name
    const char* fname = argv[3];
    // triangulation settings
    unsigned int nSamp = static_cast<unsigned int>(atoi(argv[4]));
    unsigned int nRuns = static_cast<unsigned int>(atoi(argv[5]));
    jcompress  ( fname, nSamp, nRuns, order, useMultiChannel, true );
  }
  else if (mode == 1)
  {
    if (argc != 7 && argc != 5) { return usage(argv); }
    const char *channel1 = 0, *channel2 = 0, *channel3 = 0, *out_fmt;
    if (useMultiChannel)
    {
      if (argc != 7) { return usage(argv); }

      channel1 = argv[3]; channel2 = argv[4]; channel3 = argv[5];
      out_fmt = argv[6];
    }
    else
    {
      if (argc != 5) { return usage(argv); } 
      channel1 = argv[3]; out_fmt = argv[4];
    }
    jdecompress ( useMultiChannel, out_fmt, channel1, channel2, channel3 );
  }
  else
  {
    return usage(argv);
  }

  return EXIT_SUCCESS;
}
