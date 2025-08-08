#pragma once

#include <list>
#include <iostream>
#include <vector>

struct Bucket{
    struct Bucket *next;
    struct Bucket *pre;
    std::vector<int> BucketList;
    int BucketId;
    Bucket();
    explicit Bucket(int id);
    Bucket(Bucket *prePoint,int id);
    Bucket(int id,Bucket *nextPoint);
    void DeleteNode(int node,int MGIndex,int direction);
    ~Bucket();
};


Bucket::Bucket() {
	next = nullptr;
	pre = nullptr;
	BucketId = 0;
}

Bucket::Bucket(int id) {
	next = nullptr;
	pre = nullptr;
	BucketId = id;
}

Bucket::Bucket(Bucket *prePoint, int id) {
	next = prePoint->next;
	pre = prePoint;
	prePoint->next->pre = this;
	prePoint->next = this;
	BucketId = id;
}

Bucket::Bucket(int id, Bucket *nextPoint) {
	next = nextPoint;
	pre = nextPoint->pre;
	nextPoint->pre->next = this;
	nextPoint->pre = this;
	BucketId = id;
}

Bucket::~Bucket() {
	if (this->next != NULL)this->next->pre = this->pre;
	if (this->pre != NULL)this->pre->next = this->next;
}

void Bucket::DeleteNode(int node, int MGIndex, int direction) {
	Bucket *point;
	std::vector<int>::iterator iter;
	if (direction == 1) {
		point = this->next;
		while (point->BucketId < MGIndex)point = point->next;
		if (point->BucketId == MGIndex) {
			for (iter = point->BucketList.begin(); iter != point->BucketList.end(); ++iter) {
				if (*iter == node)break;
			}
			if (iter == point->BucketList.end()) {
				std::cout << "Error In Bucket! 1 " << std::endl;
				exit(-1);
			}
			else {
				point->BucketList.erase(iter);
			}
		}
		else {
			std::cout << "Error In Bucket! 2" << std::endl;
			exit(-1);
		}
	}
	else {
		point = this->pre;
		while (point->BucketId > MGIndex)point = point->pre;
		if (point->BucketId == MGIndex) {
			for (iter = point->BucketList.begin(); iter != point->BucketList.end(); ++iter) {
				if (*iter == node)break;
			}
			if (iter == point->BucketList.end()) {
				std::cout << "Error In Bucket! 3" << std::endl;
				exit(-1);
			}
			else {
				point->BucketList.erase(iter);
			}
		}
		else {
			std::cout << "Error In Bucket! 4" << std::endl;
			exit(-1);
		}
	}
	if (point->BucketList.empty()) {
		delete point;
	}
}