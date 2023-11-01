#include "checkmatrix.hpp"

using std::size_t;

namespace code::LDPC {

const char *CheckMatrix_regular<252,504,6>::path = "H.txt";
const char *CheckMatrix_regular<256,512,6>::path = "36_512.txt";
const char *CheckMatrix_regular<512,1024,6>::path = "36_1024.txt";
const char *CheckMatrix_regular<1024,2048,6>::path = "36_2048.txt";
const char *CheckMatrix_regular<5000,10000,6>::path = "H2.txt";
const char *CheckMatrix_irregular<4999,10000>::path = "H3.txt";

template<> auto getCheckMatrix<252,504>(){return CheckMatrix_regular<252,504,6>();}
template<> auto getCheckMatrix<256,512>(){return CheckMatrix_regular<256,512,6>();}
template<> auto getCheckMatrix<512,1024>(){return CheckMatrix_regular<512,1024,6>();}
template<> auto getCheckMatrix<1024,2048>(){return CheckMatrix_regular<1024,2048,6>();}
template<> auto getCheckMatrix<5000,10000>(){return CheckMatrix_regular<5000,10000,6>();}
template<> auto getCheckMatrix<4999,10000>(){return CheckMatrix_irregular<4999,10000>();}

}
