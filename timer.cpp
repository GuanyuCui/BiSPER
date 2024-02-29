#include "timer.h"
#include <iostream>

Timer::Timer() : total_duration(0){}

void Timer::start()
{
    this -> start_timestamp = std::chrono::high_resolution_clock::now();
}

void Timer::stop()
{
    this -> stop_timestamp = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop_timestamp - start_timestamp);
    this -> total_duration += duration;
}

void Timer::clear()
{
    this -> total_duration = std::chrono::nanoseconds(0);
}

std::chrono::nanoseconds::rep Timer::get_count() const
{
    return this -> total_duration.count();
}

double Timer::get_duration(const std::string & unit) const
{
    double denominator = std::map<std::string, double>({{"s", 1e9}, {"ms", 1e6}, {"us", 1e3}, {"ns", 1.0}})[unit]; 
	return static_cast<double>(this -> total_duration.count() / denominator);
}