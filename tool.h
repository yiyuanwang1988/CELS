#pragma once

#include "Graph.h"
#include <ctime>
#include <random>
#include <cstdlib>

class randTool {
public:
	std::default_random_engine randomE;
	randTool();
	unsigned int timeseed = 0;
	explicit randTool(unsigned int rseed);
	~randTool() = default;

	/*产生[low,high]的随机数*/
	int RandInt(int low, int high);
	/*产生[low,high)的随机浮点数*/
	double RandDouble(double low, double high);
};

class timeTool {
public:
	clock_t SearchStart, SearchFinish;
	double run_time;
	int CutoffTime;
	explicit timeTool(int cutoff);
	~timeTool() = default;

	void StartClock();
	void FinishClock();
	bool OverCutoff();
	void RunTime();
	double GetTime();
	double GetCurTime();
};


randTool::randTool() {
	timeseed = time(nullptr);
	randomE.seed(timeseed);
}

randTool::randTool(unsigned int rseed) {
	timeseed = rseed;
	randomE.seed(rseed);
}

timeTool::timeTool(int cutoff) {
	CutoffTime = cutoff * CLOCKS_PER_SEC;
}

void timeTool::StartClock() {
	SearchStart = clock();
	SearchFinish = clock();
}

void timeTool::FinishClock() {
	SearchFinish = clock();
}

bool timeTool::OverCutoff() {
	return (clock() - SearchStart) > CutoffTime;
}

void timeTool::RunTime() {
	run_time = static_cast<double>(SearchFinish - SearchStart) / CLOCKS_PER_SEC;
}

double timeTool::GetTime() {
	return static_cast<double>(SearchFinish - SearchStart) / CLOCKS_PER_SEC;
}

double timeTool::GetCurTime() {
	return static_cast<double>(clock() - SearchStart) / CLOCKS_PER_SEC;
}

int randTool::RandInt(int low, int high) {
	std::uniform_int_distribution<int> U(low, high);
	return U(randomE);
}

double randTool::RandDouble(double low, double high) {
	std::uniform_real_distribution<double > U(low, high);
	return U(randomE);
}




