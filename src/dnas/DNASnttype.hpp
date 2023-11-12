#ifndef INCLUDE_GUARD_dnas_DNASnttype
#define INCLUDE_GUARD_dnas_DNASnttype

#include <array>
#include <concepts>
#include <cstdint>

namespace code::DNAS {

template<std::uint8_t ATGC=0x1B>
class nucleotide_t{
	static_assert((((ATGC>>6)&3)!=((ATGC>>4)&3))&&(((ATGC>>6)&3)!=((ATGC>>2)&3))&&(((ATGC>>6)&3)!=(ATGC&3))&&(((ATGC>>4)&3)!=((ATGC>>2)&3))&&(((ATGC>>4)&3)!=(ATGC&3))&&(((ATGC>>2)&3)!=(ATGC&3)));//ATGCに重複がない
	std::uint8_t val;
public:
	constexpr nucleotide_t():val(0ui8){}
	constexpr nucleotide_t(int rhs):val(static_cast<std::uint8_t>(rhs)&3ui8){}
	constexpr nucleotide_t(const nucleotide_t &rhs):val(rhs.val){}
	constexpr auto operator+(const nucleotide_t &rhs) const{return nucleotide_t(val+rhs.val);}
	constexpr auto operator+(int rhs) const{return nucleotide_t(val+static_cast<std::uint8_t>(rhs));}
	constexpr auto operator-(const nucleotide_t &rhs) const{return nucleotide_t(val-rhs.val);}
	constexpr auto operator-(int rhs) const{return nucleotide_t(val-static_cast<std::uint8_t>(rhs));}
	void operator=(const nucleotide_t &rhs){val=rhs.val;}
	void operator=(int rhs){val=rhs&3ui8;}
	void operator+=(const nucleotide_t &rhs){val=(val+rhs.val)&3ui8;}
	void operator+=(int rhs){val=(val+static_cast<std::uint8_t>(rhs))&3ui8;}
	void operator-=(const nucleotide_t &rhs){val=(val-rhs.val)&3ui8;}
	void operator-=(int rhs){val=(val-static_cast<std::uint8_t>(rhs))&3ui8;}
	constexpr auto operator^(int rhs) const{return nucleotide_t(val^static_cast<std::uint8_t>(rhs));}
	constexpr operator std::uint8_t() const{return val;}

	constexpr bool msb() const{return static_cast<bool>(val>>1);}
	constexpr bool lsb() const{return static_cast<bool>(val&1);}
	constexpr bool is_AT() const;
	constexpr bool is_GC() const;
};

constexpr bool nucleotide_t<0x1B>::is_AT() const{return !static_cast<bool>(val>>1);}
constexpr bool nucleotide_t<0x1B>::is_GC() const{return static_cast<bool>(val>>1);}

constexpr bool nucleotide_t<0x27>::is_AT() const{return !static_cast<bool>(val&1);}
constexpr bool nucleotide_t<0x27>::is_GC() const{return static_cast<bool>(val&1);}

template<std::floating_point T, std::uint8_t ATGC=0x1B>
class nucleotide_p{
	static constexpr std::uint8_t nA = (ATGC>>6)&3, nT = (ATGC>>4)&3, nG = (ATGC>>2)&3, nC = ATGC&3;
	static_assert(nA!=nT && nA!=nG && nA!=nC && nT!=nG && nT!=nC && nG!=nC);//ATGCに重複がない
	std::array<T,4> likelihood;
public:
	auto &lhA(){return likelihood[nA];}
	auto &lhT(){return likelihood[nT];}
	auto &lhG(){return likelihood[nG];}
	auto &lhC(){return likelihood[nC];}
	auto &lh(std::size_t i){return likelihood[i];}
	const auto &lhA() const{return likelihood[nA];}
	const auto &lhT() const{return likelihood[nT];}
	const auto &lhG() const{return likelihood[nG];}
	const auto &lhC() const{return likelihood[nC];}
	const auto &lh(std::size_t i) const{return likelihood[i];}
};

}

#endif