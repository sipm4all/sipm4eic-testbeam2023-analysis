#pragma once

namespace sipm4eic {

class lightdata {

 public:

  unsigned char device = 0;
  unsigned char index = 0;
  unsigned char coarse = 0;
  unsigned char fine = 0;
  unsigned char tdc = 0;

  int chip() const { return index / 32; };
  int eoch() const { return index % 64; };
  int cindex() const { return tdc + 4 * index; };
  
  lightdata(unsigned char _device, unsigned char _index, unsigned char _coarse, unsigned char _fine, unsigned char _tdc) :
    device(_device),
    index(_index),
    coarse(_coarse),
    fine(_fine),
    tdc(_tdc) { };
    
};

} /** namespace sipm4eic **/