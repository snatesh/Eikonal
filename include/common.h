#ifndef _COMMON_H
#define _COMMON_H


// flattened index into 2D array
inline unsigned int at(unsigned int i, unsigned int j, const unsigned int Nx)
{
  return i + Nx * j;
}

// flattened index into 3D array
inline unsigned int at(unsigned int i, unsigned int j,unsigned int k,\
                             const unsigned int Nx, const unsigned int Ny)
{
  return i + Nx * (j + Ny * k);
}


#endif
