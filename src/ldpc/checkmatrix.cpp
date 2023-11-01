#include "checkmatrix.hpp"

#define SET_CheckMatrix_regular(P,S,C,W) const char *CheckMatrix_regular<S,C,W>::path = P; template<> auto getCheckMatrix<S,C>(){return CheckMatrix_regular<S,C,W>();}
#define SET_CheckMatrix_irregular(P,S,C) const char *CheckMatrix_irregular<S,C>::path = P; template<> auto getCheckMatrix<S,C>(){return CheckMatrix_irregular<S,C>();}

namespace code::LDPC {

SET_CheckMatrix_regular("H.txt",252,504,6)
SET_CheckMatrix_regular("36_512.txt",256,512,6)
SET_CheckMatrix_regular("36_1024.txt",512,1024,6)
SET_CheckMatrix_regular("36_2048.txt",1024,2048,6)
SET_CheckMatrix_regular("H2.txt",5000,10000,6)
SET_CheckMatrix_irregular("H3.txt",4999,10000)

}

#undef SET_CheckMatrix_regular
#undef SET_CheckMatrix_irregular