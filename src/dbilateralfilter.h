
#ifndef DBILATERALFILTER_H
#define DBILATERALFILTER_H

#include "dimage.h"


class DProgress; // forward declaration
///This class provides bilateral filter functionality for DImage objects
/*! 
  \verbatim

  The code for the default filter type (DBilatFilt_paris), is based
  on code provided online by Sylvain Paris and Frédo Durand.  I have
  modified their code to work within this library framework.

  **Please note** that they request you cite their paper if you use
  their code in your research, so you will want to look up the full 
  citations for the following:
     
    Sylvain Paris and Frédo Durand, "A Fast Approximation of the Bilateral
    Filter using a Signal Processing Approach"  European Conference on
    Computer Vision (ECCV'06).

    Sylvain Paris and Frédo Durand, "A Fast Approximation of the Bilateral
    Filter using a Signal Processing Approach"  MIT technical report 2006
    (MIT-CSAIL-TR-2006-073).


  Also, the following copyrights and licenses apply to various parts of
  the code:
  
    Copyright (c) 2006, Sylvain Paris and Frédo Durand

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

    ---------

    Copyright (c) 2004, Sylvain Paris and Francois Sillion
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
    
        * Redistributions in binary form must reproduce the above
        copyright notice, this list of conditions and the following
        disclaimer in the documentation and/or other materials provided
        with the distribution.

        * Neither the name of ARTIS, GRAVIR-IMAG nor the names of its
        contributors may be used to endorse or promote products derived
        from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    ---------------------------------------------

    
  \endverbatim
*/

class DBilateralFilter {
public:

  enum DBilatFiltType{
    DBilatFilt_paris, ///<Use Sylvain Paris' fast approximation method
    DBilatFilt_paris_fft, ///<Use Sylvain Paris' method with FFT convolution (Not Yet Implemented)
    DBilatFilt_slow,///<Slow brute-force method (Not Yet Implemented)
  };

  DBilateralFilter();///< defaults: sigma_s=16, sigma_r=0.1, DBilatFilt_paris
  DBilateralFilter(double sigma_s, double sigma_r, DBilatFiltType filtType);///< spatial sigma (in pixels), range sigma (intensity scale is [0,1.0])
  ~DBilateralFilter();

  void setSigmas(double sigma_s, double sigma_r);
  double getSigmaS();
  double getSigmaR();
  void setType(DBilatFiltType filtType);
  DBilatFiltType getType();
  void setNumThreads(int numThreads);
  int getNumThreads();
  
  DImage filterImage(const DImage &imgSrc, bool fAlreadyPadded = false,
		     DProgress *pProg = NULL);
  void filterImage_(DImage &imgDst,
		    const DImage &imgSrc, bool fAlreadyPadded = false,
		    DProgress *pProg = NULL);

  static void bilateralFilterImage(DImage &imgDst, const DImage &imgSrc,
				   bool fAlreadyPadded = false,
				   double sigma_s = 16., double sigma_r = 0.1,
				   DBilatFiltType filtType = DBilatFilt_paris,
				   DProgress *pProg = NULL,int numThreads = 1);

private:
  DBilatFiltType _bilatFiltType;
  double _sigma_s;
  double _sigma_r;
  int _numThreads;

//   static void bilatFilt_


//   static void medFiltHuang_u8(DImage &imgDst, const DImage &imgSrc,
// 			      int radiusX, int radiusY,
// 			      int wKern, int hKern, D_uint8 *rgKern,
// 			      int numKernPxls, int *rgRightEdge,
// 			      DProgress *pProg = NULL,
// 			      int progStart = 0, int progMax = 1,
// 			      int threadNumber = 0, int numThreads = 1);
//   static void medFiltHuang_u8_square(DImage &imgDst, const DImage &imgSrc,
// 				     int radiusX, int radiusY,
// 				     int wKern, int hKern, D_uint8 *rgKern,
// 				     int numKernPxls,
// 				     DProgress *pProg = NULL,
// 				     int progStart = 0, int progMax = 1,
// 				     int threadNumber = 0, int numThreads = 1);
//   static void* DMedianFilter_Huang8threadWrap(void* params);

  /// copy constructor is private so nobody can use it
  DBilateralFilter(const DBilateralFilter &src);
};



#endif
