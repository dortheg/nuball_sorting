#include <fstream>
#include <sstream>
#include <iomanip>
#include "CubeDDT.hxx"


CubeDDT::CubeDDT(string name, int nchans, int tchans, int nsperchan) //default constructor
{
fNChans=nchans;
fTChans=tchans;
fnsPerChan=nsperchan;
fName=name;
double checksize=nchans*nchans*tchans;
// Is the cube greater than the maximum theoretical memory size allowed with int index?
if (checksize > 4294967295) {cout << "Error. Cube index is bigger than 4billion"; exit(1);}
else {cout << "Creating cube of size " << double(checksize*4)/1000000000.0 << " Gb" << endl;}
fCubeDDTSize=nchans*nchans*tchans;
fCubeDDT=new float[fCubeDDTSize];
//cout << "Check size (elements) " << size(fCubeDDT) << endl;
//cout << "Check size (bytes) " << sizeof(fCubeDDT) << endl;
}

//________________________________________________________________________

CubeDDT::~CubeDDT() //destructor
{
 delete [] fCubeDDT; //Erase the buffer
}
//________________________________________________________________________
void CubeDDT::Fill(unsigned short i, unsigned short j, unsigned short t)
{
if ((i < fNChans) && (j < fNChans) && (t < fTChans))	
	{			
	fIndex=i+(j*fNChans)+(t*fNChans*fNChans);
	fCubeDDT[fIndex]++;
	fIndex=j+(i*fNChans)+(t*fNChans*fNChans); //Symmetric incrementation
	fCubeDDT[fIndex]++;
	}
}
//________________________________________________________________________
void CubeDDT::Filla(unsigned short i, unsigned short j, unsigned short t)
{
if ((i < fNChans) && (j < fNChans) && (t < fTChans))	
	{			
	fIndex=i+(j*fNChans)+(t*fNChans*fNChans);
	fCubeDDT[fIndex]++;
	}
}
//________________________________________________________________________
void CubeDDT::Read(string fname)
{
fName=fname;
cout << "\n";
cout << "Reading binary ddtime cube " << endl;
ifstream infile(fName.c_str(), ios::in | ios::binary);
if (!infile) {cout << "Cannot read file " << fName << endl; exit(1);}
infile.read((char *) fCubeDDT, fCubeDDTSize*sizeof(float));
cout << infile.gcount() << " bytes read\n";
infile.close();
}
//________________________________________________________________________
void CubeDDT::Write(string fname)
{
this->Integral();
fName=fname;
cout << "Writing binary ddtime cube " << fName << endl;
ofstream ofile(fName.c_str(), ios::out | ios::binary);
ofile.write((char *) fCubeDDT, fCubeDDTSize*sizeof(float)); //two bytes per entry. bsize entries
ofile.close();
}
//________________________________________________________________________
double CubeDDT::Integral()
{
double sum=0;
for (int i=0; i < fNChans; i++)
	{
	for (int j=0; j < fNChans; j++)
		{
		for (int k=0; k < fTChans; k++)
			{
			fIndex=j+(i*fNChans)+(k*fNChans*fNChans);
			sum+=fCubeDDT[fIndex];
			}
		}
	}
cout << "Cube has " << sum << " counts" << endl;
return sum;
}
//________________________________________________________________________
void CubeDDT::Set(unsigned short i, unsigned short j, unsigned short t, float val)
{
if ((i < fNChans) && (j < fNChans) && (t < fTChans))	
	{			
	fIndex=i+(j*fNChans)+(t*fNChans*fNChans);
	fCubeDDT[fIndex]=val;
	//fIndex=j+(i*fNChans)+(t*fNChans*fNChans); //Symmetric incrementation
	//fCubeDDT[fIndex]=val;
	}
}

//________________________________________________________________________
float CubeDDT::Get(unsigned short i, unsigned short j, unsigned short t)
{
if ((i < fNChans) && (j < fNChans) && (t < fTChans))	
	{			
	fIndex=i+(j*fNChans)+(t*fNChans*fNChans);
	return fCubeDDT[fIndex];
	}
else return 0;
}
//________________________________________________________________________
/*
void CubeDDT::WriteProjections() //Make an energy projection spectrum for each time value
{
Spec* proj=new Spec(fNChans);
for (int i=0; i < fTChans; i++)
	{
	string fname=Form("proj%d.spe",i);
	for (int j=0; j < fNChans; j++)
		{
		double sum=0;
			{
			for (int k=0; k <fNChans; k++){sum+=this->Get(j,k,i);}
			proj->SetChan(j,sum);
			}
		}
	proj->Write(fname);
	}
}
*/
