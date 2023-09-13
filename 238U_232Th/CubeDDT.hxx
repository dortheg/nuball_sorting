#ifndef _CubeDDT_
#define _CubeDDT_ 

#include <math.h>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include "TH2F.h"

using namespace std;

class CubeDDT
{
 public:
	CubeDDT(string name="ddtime.bin",int nchans=4096,int tchans=40,int nsperchan=10);	//@- Default constructor
	~CubeDDT(); 			//@- Normal destructor
	void Write(string name=""); //write the CubeDDT out
	void Read(string Filename=""); //Read matrix from file
	float Get(unsigned short i, unsigned short j, unsigned short t);
	void Set(unsigned short i, unsigned short j, unsigned short t, float val);
	void Fill(unsigned short i, unsigned short j, unsigned short t); //increment channel by 1 count
	void Filla(unsigned short i, unsigned short j, unsigned short t);
	double Integral();
	//void WriteProjections();

 protected:
	unsigned short 	fNChans;
	unsigned short 	fTChans;
	unsigned int 	fCubeDDTSize;
	unsigned short fnsPerChan;
	string 		fName;
	float*  fCubeDDT;
	double		fTotalCounts;
	unsigned int fIndex;
};
#endif
