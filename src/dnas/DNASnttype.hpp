#ifndef __DNASMYTYPE__
#define __DNASMYTYPE__

#include <cstdint>

namespace code::DNAS {

class nucleotide_t{
	std::uint8_t val;
public:
	explicit nucleotide_t(std::uint8_t val):val(val){}
	explicit nucleotide_t(int val):val(static_cast<std::uint8_t>(val)){}
	explicit nucleotide_t():val(0ui8){}
	inline auto operator+(nucleotide_t rhs){return nucleotide_t((val+rhs.val)&3ui8);}
	inline auto operator+(int rhs){return nucleotide_t((val+static_cast<std::uint8_t>(rhs))&3ui8);}
	inline void operator=(nucleotide_t rhs){val=rhs.val&3ui8;}
	inline void operator=(int rhs){val=static_cast<std::uint8_t>(rhs&3);}
	inline void operator+=(nucleotide_t rhs){val+=rhs.val&3ui8;}
	inline void operator+=(int rhs){val+=static_cast<std::uint8_t>(rhs&3);}
	inline auto msb(){return static_cast<bool>((val>>1)&1);}
	inline auto lsb(){return static_cast<bool>((val)&1);}
};

}

#endif