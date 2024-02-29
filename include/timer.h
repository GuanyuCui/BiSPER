#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <chrono>
#include <map>

class Timer
{
    public:
        Timer();
        ~Timer() = default;
        
        void start();
        void stop();
        void clear();

        double get_duration(const std::string & unit = "ms") const;
        std::chrono::nanoseconds::rep get_count() const;

    private:
        decltype(std::chrono::high_resolution_clock::now()) start_timestamp;
        decltype(std::chrono::high_resolution_clock::now()) stop_timestamp;
        
        std::chrono::nanoseconds total_duration;
};

#endif