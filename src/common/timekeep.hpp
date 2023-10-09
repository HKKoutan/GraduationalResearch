#ifndef __util_TIMEKEEP__
#define __util_TIMEKEEP__

#include <iostream>
#include <chrono>
#include <vector>

namespace util {

class Timekeep {
	std::vector<std::chrono::milliseconds> times; //記録した時間
	std::chrono::system_clock::time_point prev; //計測開始時刻

public:
	inline void start() noexcept{//計測開始
		prev = std::chrono::system_clock::now();
	}
	inline void stop(){//計測終了
		times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - prev));
		std::cout<<"time"<<std::endl;
		for(auto &t: times) std::cout<<t.count()<<"ms"<<std::endl;
	}
	inline void split(){//今の時間を記録して新たに計測
		times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - prev));
		prev = std::chrono::system_clock::now();
	}
};

}

#endif