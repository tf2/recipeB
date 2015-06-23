#include<set>
#include<map>
#include<math.h>
#include<vector>
#include<sstream>
#include<fstream>
#include<iostream>

using namespace std;

#include <R.h> 
#include <Rdefines.h> 
#include <Rinternals.h>
#include <Rinterface.h>

extern "C" {

/* Declarations */
const int FIELD_LENGTH=8;
double d; double dd;
double recip[] = {0, 0};
double rmer=0.5; double bmer=0.9; double fmer=0.5;
map < string, vector<string> > mdata;
map < int, vector<double> > pdata;
vector<double> all_types;

void readFile(char *fname);
void executeRecipA(char* outname, char* logname);

int breakdist1 (int a, int b);
int breakdist2 (int a, int b);
int recipoverlap(int a1, int a2, int b1, int b2);

typedef double (*recipA)[FIELD_LENGTH];
vector< vector < double > > recibreak(double Aarray[][FIELD_LENGTH], int len);

/* Functions */
int breakdist1(int a, int b) {
	return((int)(max(a,b) - abs(a-b)*bmer));
}

int breakdist2(int a, int b) {
	return((int)(min(a,b) + abs(a-b)*bmer));
}

int recipoverlap(int a1, int a2, int b1, int b2) {	
	if(a2<b1 | a1>b2) { recip[0]=0; recip[1]=0; return 0; }
	if(a1==b1 & a2 == b2) { recip[0]=1; recip[1]=1; return 1; }
	if(a1<=b1 & a2>b2) { d=b2-b1; dd=a2-a1; d=(d/dd); recip[0]=d; recip[1]=1; return 1; }
	if(a1>=b1 & a2<b2) { recip[0]=1; d=a2-a1; dd=b2-b1; d=(d/dd); recip[1]=d; return 1; }
	
	if(a2<=b2 & a1>b1) { recip[0]=1; d=a2-a1; dd=b2-b1; d=(d/dd); recip[1]=d; return 1; }
	
	if(a1==b2 & a2 > b2) { d=1; dd=a2-a1; d=(d/dd); recip[0]=d; d=1; dd=b2-b1; d=(d/dd); recip[1]=d; return 1; }
	if(a2==b1 & a1 < b1) { d=1; dd=a2-a1; d=(d/dd); recip[0]=d; d=1; dd=b2-b1; d=(d/dd); recip[1]=d; return 1; }
	if(a1<b1 & a2>b1 & a2<=b2) { d=min(a2-b1, a2-a1); dd=max(a2-b1, a2-a1); d=(d/dd); recip[0]=d; d=min(a2-b1, b2-b1); dd=max(a2-b1, b2-b1); d=(d/dd); recip[1]=d; return 1; }
	if(a2>b2 & a1<b2 & a1>=b1) { d=min(a2-a1, b2-a1);  dd=max(a2-a1, b2-a1); d=(d/dd); recip[0]=d; d=min(b2-b1, b2-a1); dd=max(b2-b1, b2-a1); d=(d/dd); recip[1]=d; return 1; }
}

vector< vector < double > > recibreak(double Aarray[][FIELD_LENGTH], int len) {	
	vector<double> chr; vector<double> start; vector<double> stop; vector<double> obs; vector<double> freq; vector<double> typ; vector<double> stud;
	for(int i=0;i<len;i++) {
		chr.push_back(Aarray[i][0]); start.push_back(Aarray[i][1]); stop.push_back(Aarray[i][2]); obs.push_back(Aarray[i][3]); freq.push_back(Aarray[i][4]); typ.push_back(Aarray[i][5]); stud.push_back(Aarray[i][6]);
	}
	int pin1=0; int first = 1; int count = 0;
	vector< vector < double > > segments;
	while(start.size()>0) {
	vector < double > features;
		if(first) { 
			recipoverlap((int)start[0],(int)stop[0],(int)start[1],(int)stop[1]); 
			first = 0;
		}	
		features.push_back(start[pin1]); features.push_back(stop[pin1]); features.push_back(obs[pin1]); features.push_back(freq[pin1]); features.push_back(typ[pin1]); features.push_back(stud[pin1]); features.push_back(recip[0]); features.push_back(recip[1]); 
		start.erase(start.begin()); stop.erase(stop.begin()); obs.erase(obs.begin()); freq.erase(freq.begin()); typ.erase(typ.begin()); stud.erase(stud.begin());
		for(int j=0;j<features.size();j+=j+FIELD_LENGTH) {
		for(int pin2=0;pin2<start.size();pin2++) {
			recipoverlap((int)start[pin2],(int)stop[pin2], (int)features[j],(int)features[j+1]);
				if(recip[0]>rmer & recip[1]>rmer) {
					features.push_back(start[pin2]); features.push_back(stop[pin2]); features.push_back(obs[pin2]); features.push_back(freq[pin2]); features.push_back(typ[pin2]); features.push_back(stud[pin2]); features.push_back(recip[1]); features.push_back(recip[0]); 
					start.erase(start.begin()+(pin2)); stop.erase(stop.begin()+(pin2)); obs.erase(obs.begin()+(pin2)); freq.erase(freq.begin()+(pin2)); typ.erase(typ.begin()+(pin2)); stud.erase(stud.begin()+(pin2)); first = 1;
					pin2--;
				}
			}
		}
		segments.push_back(features);
	}
	return(segments);
}

void readFile(char *fname) {
  	char line[100]; FILE *infile; infile=fopen(fname, "r"); 
  		while(fgets(line, sizeof(line), infile)!=NULL) {
  			line[strlen(line)-1]='\0'; string ss(line); string s=strtok(line, "\t"); mdata[s].push_back(ss);
  		}
 	fclose(infile);
}

void tologfile (vector< vector < double > > segments, char* logname) {
	map < string, vector<string> > rdata;
	for(int i=0; i<segments.size();i++) {
		string ss; std::stringstream out; vector < double > features = segments[i]; out << "cluster: " << i << endl; 
		for(int j=0;j<features.size();j++) {
			if((j+1) % FIELD_LENGTH == 0)  out << features[j] << endl;
			else out << features[j] << "\t";
		}
		ss = out.str(); rdata[ss].push_back(ss);	
	}
	ofstream filehand; filehand.open (logname, ios::out | ios::app);
	for(map<string, vector<string> >::const_iterator it = rdata.begin(); it != rdata.end(); it++) {
		string s=it->first;	filehand << s << endl;
	}
	filehand.close();		
}

void tooutfile (char* outname, double chromosome, vector < double > cluster) {	
		string ss; 
		std::stringstream out; 
		out << (int)chromosome << "\t" << (int)cluster[0] << "\t" << (int)cluster[1] << "\t";
		for(int i=2;i<cluster.size()-1;i++) {
			out <<  cluster[i] << "\t";
		}
		out <<  cluster[cluster.size()-1] << endl;
		ss = out.str();
		ofstream filehand; filehand.open (outname, ios::out | ios::app);
		filehand << ss;
		filehand.close();		
}

vector<double> kweights(vector<double> obs, vector<double> freq, vector<double> typ, vector<double> stud) {
	vector<double> type_infos; stringstream ss;
	int cobs; double cfreq; char ctypes[500]; char cstudies[500]; double samps=0;
	for(map<int, vector<double> >::const_iterator it = pdata.begin(); it != pdata.end(); it++) {
		double studset = it->first;	
		vector<double> params=it->second; int check =0;	
		for(int i=0;i<stud.size();i++) {
			if(studset == stud[i]) { 
				check = 1;
			}
		}
		if(check ==1) {
			// Here we make the comibned study tag and load params connected to each study (we dont use sensi yet).
			ss << studset;
			vector<double> params = pdata.find((int)studset)->second;
			double samplesize = params[0]; double sensitivity = params[1];
			samps+=samplesize;
		}
	}
	double studytag = atoi(ss.str().c_str()); double ttobs=0; double whattypes = typ[0];
	for(int j=0;j<all_types.size();j++) {
		double tfreq=0; double tobs=0;
		for(int i=0;i<stud.size();i++) {
		if(whattypes != typ[i]) { whattypes=0; }
			if(typ[i] == all_types[j]) {
				if(freq[i] > tfreq) { tfreq=freq[i]; }
				tobs += obs[i];
			}
		}
		// here is the standard error
		double se = sqrt((samps-tobs)/(tobs*samps));
		
		// this is to guard against the rubbish case where we have a set that has only one sample but two (or more) features were merged from that set! (i.e. were overlaping by the defined ammount).
		// Otherwise we calculate the weighted average - tobs and samps were defined above.
		if(tfreq==0 & tobs>1) { tfreq = tobs/samps; }
		if(tobs <= 1) { tfreq = 0; se = 1; } 
		else if (tobs>samps | tobs==samps) { tfreq = 0; se = 1; }
		ttobs+=tobs;
		type_infos.push_back(tobs); type_infos.push_back(tfreq); type_infos.push_back(se);
	}
	double ttfreq=ttfreq = ttobs/samps; double sse = sqrt((samps-ttobs)/(ttobs*samps));
	if (ttobs <= 1) { ttobs=1; ttfreq = 0; sse = 1; }
	if (ttobs>samps) { ttobs = samps; ttfreq = 0; sse = 1; } //ttobs = samps;
	if (ttobs==samps) { ttfreq = 0; sse = 1; }
	type_infos.push_back(ttobs); type_infos.push_back(ttfreq); type_infos.push_back(sse);
	type_infos.push_back(whattypes); type_infos.push_back(samps); type_infos.push_back(studytag);
	
	return(type_infos);
}

vector<double> ccluster(vector<double> cluster) {
	vector<double> obs; vector<double> freq; vector<double> typ; vector<double> stud;
	double inner_startpos = cluster[0]; double outer_startpos = cluster[0];
	double inner_stoppos = cluster[1];	double outer_stoppos = cluster[1];
	for(int j=0;j<cluster.size();j=j+FIELD_LENGTH) {
		if(cluster[j] > inner_startpos) inner_startpos = cluster[j];
		if(cluster[j] < outer_startpos) outer_startpos = cluster[j];
		if(cluster[j+1] < inner_stoppos) inner_stoppos = cluster[j+1];
		if(cluster[j+1] > outer_stoppos) outer_stoppos = cluster[j+1];
		obs.push_back(cluster[j+2]); freq.push_back(cluster[j+3]);
		typ.push_back(cluster[j+4]); stud.push_back(cluster[j+5]);
	}
	int break_startpos=0; int break_stoppos=0;
	if(obs.size()>1) {
		break_startpos = breakdist1((int)inner_startpos, (int)outer_startpos); 
		break_stoppos = breakdist2((int)inner_stoppos, (int)outer_stoppos);
	} else {
		break_startpos= inner_startpos;
		break_stoppos= inner_stoppos;
	}
	
	vector<double> merge_info;
	merge_info.push_back(break_startpos);
	merge_info.push_back(break_stoppos);
	
	vector<double> type_info = kweights(obs, freq, typ, stud);
	for(int i=0;i<type_info.size();i++) {
		merge_info.push_back(type_info[i]);
	}
	
	return(merge_info);
}

void cmerge(char* outname, double chromosome, vector< vector < double > > segments) {
	for(int i=0; i<segments.size();i++) {
		vector < double > features = segments[i]; vector< double > cluster;
			for(int j=0;j<features.size();j=j+FIELD_LENGTH) {
				for(int k=0;k<FIELD_LENGTH;k++) {
					cluster.push_back(features[j+k]);
				}
			}
		vector<double> merge_info = ccluster(cluster);
		tooutfile(outname, chromosome, merge_info);
	}
}

void executeRecipA(char* outname, char* logname) {	
	for(map<string, vector<string> >::const_iterator it = mdata.begin(); it != mdata.end(); it++) {
		vector<string> features=it->second;	
		int len=features.size(); double Aarray[len][FIELD_LENGTH];	
			for(int i=0;i<len;i++) {
				int count=0; string fea(features[i]); char* f=new char[fea.size()+1];
				f[fea.size()]=0; memcpy(f,fea.c_str(),fea.size()); string ss=strtok(f, "\t");	
					while (f!=NULL) {
			 			double t=atof(f); Aarray[i][count]=t; f=strtok(NULL, "\t"); count++;
  					}
		}
		vector< vector < double > > segments = recibreak(Aarray, len); 
		cmerge(outname, Aarray[0][0], segments);
		tologfile(segments, logname);
	}
}

void overlord(double rep1, double rep2, int a, int b, int *c, int *d, int *ss, double *r1, double *r2, int *indexes, int N) {	
	double r = 0; double rr=0; int pin=0;
	for(int j=0;j<N;j++) {
		recipoverlap(a,b,c[j],d[j]);
		if(recip[0]> rep1 & recip[1] > rep2) {
			r1[pin] = recip[0]; r2[pin] = recip[1]; indexes[pin] = ss[j];
			pin++;
		}
		if(c[j]>b) { break; }
	}
}	
	
void sover(int a, int b, int *c, int *d, double *r1, double *r2, int i, int N) {
	double r = 0; double rr=0;
	for(int j=0;j<N;j++) {
		recipoverlap(a,b,c[j],d[j]);
		if(recip[0]> r) { 
			r=recip[0]; rr=recip[1];
		} 
		if(c[j]>b) break;
	}
	r1[i] = r; r2[i] = rr;
}
	
void tover(int a, int b, int *c, int *d, double *r1, double *r2, int i, int N) {
	double totalforward = 0; double maxbackward = 0; double r = 0; double rr=0; vector<int> starts; vector<int> stops;
	for(int j=0;j<N;j++) {
		recipoverlap(a,b,c[j],d[j]);
		if(recip[0]>0) { 
			starts.push_back(c[j]); stops.push_back(d[j]);
			if(maxbackward<recip[1]) { maxbackward =recip[1]; }
			if(totalforward<recip[0]) { totalforward =recip[0]; }
		} 
		if(c[j]>b) break;
	}
	r1[i] = totalforward; r2[i] = maxbackward;
}
	

void rover(int a, int b, int *c, int *d, double *r1, double *r2, int i, int N) {
	double r = 0; double rr=0;
	for(int j=0;j<N;j++) {
		recipoverlap(a,b,c[j],d[j]);
		if(recip[0]> r & recip[1] >rr) { 
			r=recip[0]; rr=recip[1];
		} else if (recip[0]> r & recip[1] >=rr | recip[0] >= r & recip[1] >rr) {
			r=recip[0]; rr=recip[1];
		}
		if(c[j]>b) break;
	}
	r1[i] = r; r2[i] = rr;
}

void fover(int a, int b, int *c, int *d, double *fre1, double *fre2, double *fre3, double *r1, double *r2, double *rf1, double *rf2, double *rf3, double *ty, double *rty, int i, int N) {
	double r = 0; double rr=0; double f1=-1; double f2=-1; double f3=-1; double t1=ty[0];
	for(int j=0;j<N;j++) {
		recipoverlap(a,b,c[j],d[j]);
		if(recip[0]> r & recip[1] >rr) { 
			r=recip[0]; rr=recip[1]; f1=fre1[j]; f2=fre2[j]; f3=fre3[j]; t1=ty[j];
		} else if (recip[0]> r & recip[1] >=rr | recip[0] >= r & recip[1] >rr) {
			r=recip[0]; rr=recip[1]; f1=fre1[j]; f2=fre2[j]; f3=fre3[j]; t1=ty[j];
		}
		if(c[j]>b) break;
	}
	r1[i] = r; r2[i] = rr;
	rf1[i] = f1; rf2[i] = f2; rf3[i] = f3; rty[i]=t1;
}

void cover(int a, int b, int *c, int *d, double *r1, double *r2, double *con, double rep1, double rep2, int i, int N) {
	double r = 0; double rr=0; double count =0;
	for(int j=0;j<N;j++) {
		recipoverlap(a,b,c[j],d[j]);
		if(recip[0]> rep1 & recip[1] >rep2) {
			r = recip[0]; rr = recip[1];
			count++;
		} 
	}
	r1[i] = r; r2[i] = rr;
	con[i]=count;
}


void uoverlap(double *rep1, double *rep2, int *a, int *b, int *c, int *d, int *ss, double *r1, double *r2, int *indexes, int *size1) {
	int N=(int)*size1; 
	overlord((double)*rep1, (double)*rep2, (int)a[0], (int)b[0], c, d, ss, r1, r2, indexes, N);
}	
	
void overlap(int *a, int *b, int *c, int *d, double *r1, double *r2, int *size1, int *size2) {
	int N=(int)*size1; int NN=(int)*size2;
	for( int i=0;i<N;i++) {
		sover((int)a[i], (int)b[i], c, d, r1, r2, i, NN);
	}
}

void toverlap(int *a, int *b, int *c, int *d, double *r1, double *r2, int *size1, int *size2) {
	int N=(int)*size1; int NN=(int)*size2;
	for( int i=0;i<N;i++) {
		tover((int)a[i], (int)b[i], c, d, r1, r2, i, NN);
	}
}
	
void roverlap(int *a, int *b, int *c, int *d, double *r1, double *r2, int *size1, int *size2) {
	int N=(int)*size1; int NN=(int)*size2;
	for( int i=0;i<N;i++) {
		rover((int)a[i], (int)b[i], c, d, r1, r2, i, NN);
	}
}

void foverlap(int *a, int *b, int *c, int *d, double *fre1, double *fre2, double *fre3, double *r1, double *r2, double *rf1, double *rf2, double *rf3, double *ty, double *rty, int *size1, int *size2) {
	int N=(int)*size1; int NN=(int)*size2;
	for( int i=0;i<N;i++) {
		fover((int)a[i], (int)b[i], c, d, fre1, fre2, fre3, r1, r2, rf1, rf2, rf3, ty, rty, i, NN);
	}
}

void coverlap(int *a, int *b, int *c, int *d, double *con, double *r1, double *r2, double *rep1, double *rep2, int *size1, int *size2) {
	int N=(int)*size1; int NN=(int)*size2;
	for( int i=0;i<N;i++) {
		cover((int)a[i], (int)b[i], c, d, r1, r2, con, (double)*rep1, (double)*rep2, i, NN);
	}
}


void setparm(double study, double samplesize, double sensitivity) {
	pdata[study].push_back(samplesize); pdata[study].push_back(sensitivity);
}

void setallparms(double *studies, double *samplesizes, double *sensitivities, int size) {
	for(int i=0;i<size;i++) {
		setparm(studies[i], samplesizes[i], sensitivities[i]);
	}
}

void rbreak(char **f1, char **f2, char **f3, double *re, double *be, double *fe, double *types, double *studies, double *samplesizes, double *sensitivities, int *size1, int *size2) { 
	
	int N1=(int)*size1;
	for(int i=0;i<N1;i++) {
		all_types.push_back(types[i]);	
	}

	int N2=(int)*size2;
	setallparms(studies, samplesizes, sensitivities, N2);
		
 	rmer = *re; bmer = *be; fmer = *fe; 
 	char *file = *f1; char *outfile = *f2; char *logfile = *f3;
 	readFile(file); executeRecipA(outfile, logfile);
}

}
