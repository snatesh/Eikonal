#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include<exception>
#include<stdexcept>
#include<iostream>
#include<string.h>

/* Exception handler class for C library */

class Exception : public std::runtime_error
{
  public:
    Exception(const char* msg, const char* _func, const char* _file, int _line) 
    : std::runtime_error(msg), func(_func), file(_file), line(_line) {}
    
    void getErr() const 
    {
      std::cerr << what() << ", " << func 
                << ", " << file << ", " << line << std::endl;
    }

  protected:
    const char* func;
    const char* file;
    int line;
};


inline void exitErr(const char* msg)
{
  std::cerr << msg << " Exiting ...\n";
  exit(1); 
}

#endif
