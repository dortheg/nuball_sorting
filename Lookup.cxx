//List of detectors to remove from the stream
//isbad[5]=true; //BGO looks like crap. Must remove this one
//isbad[134]=true;//Duplicate of 133. Bad FASTER configuration
//isbad[135]=true;//Duplicate of 133. Bad FASTER configuration
//isbad[136]=true;//Duplicate of 133. Bad FASTER configuration
//isbad[146]=true;//Duplicate of 145. Bad FASTER configuration
//isbad[147]=true;//Duplicate of 145. Bad FASTER configuration
//isbad[148]=true;//Duplicate of 145. Bad FASTER configuration
//isbad[250]=true;//EDEN. We don't need this.
//isbad[251]=true;//Beam pulse. We don't need this in the events either
isbad[112]=true;
isbad[117]=true;
isbad[177]=true;

for (int i=0; i < 255; i++) //cycle over all faster channels
{
int fasterchan=i;
int detectorposition=idmapr[i]; //The nth faster channel number in the list
	if (detectorposition > 0)
	{
	int detid=detectorposition; // The ID e.g. 40810 (phase 1 bgo)
	//cout << "Associating faster channel " << i << " with " << detid << endl;
	isthere[fasterchan]=true;
	int ringn=(detid/10000);
	int posi=((detid/100) % 100);
	int num=((detid/10) % 10); //0 = ge, 1 = bgo, 1 or 2 for LaBr3
	int segm=(detid % 10);
	ring[fasterchan]=ringn;
	alveole[fasterchan]=posi;
	if (ringn == 1) isring1[fasterchan]=true;
	if (ringn == 2) isring2[fasterchan]=true;
	if (ringn == 3) isring3[fasterchan]=true;
	if (ringn == 4) isring4[fasterchan]=true;

	// Ge and BGO first
	if ((ringn == 2) || (ringn == 3))
	{
	isclover[fasterchan]=true;
	if (num==0) {isge[fasterchan]=true; iscloverge[fasterchan]=true;}
	else {
	iscloverbgo[fasterchan]=true; isbgo[fasterchan]=true; isge[fasterchan]=false;
	//if (fasterchan == 5) {iscloverbgo[fasterchan]=false; isbgo[fasterchan]=false;} //turn off bgo 5
	}
	int snumber=(ringn-2)*12+posi;
	module[fasterchan]=snumber;	
	}
	
	if (ringn == 4)
	{
	if (num==0) {isge[fasterchan]=true; isphase1ge[fasterchan]=true;}
	isclover[fasterchan]=true;
	if (num==1) {isphase1bgo[fasterchan]=true; isbgo[fasterchan]=true;}
	int snumber=(ringn-2)*12+posi;
	module[fasterchan]=snumber;
	}
	
	if (ringn == 1) 
	{
	islabr3[fasterchan]=true;
	if (segm == 2) {ismadrid[fasterchan]=true;}
	}
	
	}	
	
}
iseden[250]=true;
ispulse[251]=true;
bigwalk[3]=true;
bigwalk[5]=true;
bigwalk[7]=true;
bigwalk[9]=true;
bigwalk[13]=true;
bigwalk[15]=true;
bigwalk[19]=true;

