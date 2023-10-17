#ifndef __DNASMYTYPE__
#define __DNASMYTYPE__

#include <cstdint>

namespace code::DNAS {

class nucleotide_t{
	std::uint8_t val;
public:
	nucleotide_t(std::uint8_t val):val(val&3ui8){}
	nucleotide_t():val(0ui8){}
	auto operator+(nucleotide_t rhs)const{return nucleotide_t(val+rhs.val);}
	auto operator+(int rhs)const{return nucleotide_t(val+static_cast<std::uint8_t>(rhs));}
	auto operator-(nucleotide_t rhs)const{return nucleotide_t(val-rhs.val);}
	auto operator-(int rhs)const{return nucleotide_t(val-static_cast<std::uint8_t>(rhs));}
	void operator=(nucleotide_t rhs){val=rhs.val;}
	void operator=(int rhs){val=static_cast<std::uint8_t>(rhs)&3ui8;}
	void operator+=(nucleotide_t rhs){val=(val+rhs.val)&3ui8;}
	void operator+=(int rhs){val=(val+static_cast<std::uint8_t>(rhs))&3ui8;}
	void operator-=(nucleotide_t rhs){val=(val-rhs.val)&3ui8;}
	void operator-=(int rhs){val=(val-static_cast<std::uint8_t>(rhs))&3ui8;}
	auto operator^(int rhs)const{return nucleotide_t(val^static_cast<std::uint8_t>(rhs));}
	operator std::uint8_t()const{return val;}

	auto msb()const{return static_cast<bool>(val>>1);}
	auto lsb()const{return static_cast<bool>(val);}
};

}

#endif