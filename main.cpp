#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ran2.h"

#define MAXATOM 20000		//maximum of atoms
#define SILATICE 3.5644f	//periodic of initial latice
#define KZAR 136873.83f	//e^2/(4*Pi*Eps0)   [g*A^3/ps^2]
#define MAXATOMTYPE 4		//maximum type of atoms
#define kB	8.314462f	// Boltzman constant [J/K/mol]
#define SQRTPI 1.7724538509	//sqrt(pi())
#define POTKOF	10.0	//potential koefficient
#define DelT 0.001 //integration step time [ps]
#define QNH	1.0f
#define RDFBIN  0.05f	//rdf calc stepize
#define MAXRDFBIN  2000


//speed in A/ps


float temp;
float pressure;
int Natom;
int Latice;
int PerType;
float* ax;
float* ay;
float* az;
float* vx;
float* vy;
float* vz;
float* fx;
float* fy;
float* fz;
int* AType;
int rseed=-31425;
char* grofile;
int* TypeKol;
float HBox;
float* mass;
//
float* InitMass;
float* A_p;
float* b_p;
float* ro_p;
float* q_p;

float* la_sw;
float* ga_sw;
float* rc_sw;
float* cos_sw;
int tempi;

float pnh;
float SysTemp;
float SysTime;
float SysPress;
float la_t;
int TFtemp;
int TFpress;

float* rdf;
int rdfcount;
float* properties;

float* ndens; //number density [A^-3]
float SystemDens;

void AtomPlace(float* x, float* y,float* z, int* A,int L,int &N);
void groout(const char* filename,float* x,float* y,float* z,float* vx, float* vy, float* vz ,int* T, int N,float Hbox);
void integrate(float* x, float* y,float* z,float* fx, float* fy, float* fz,float* vx, float* vy, float* vz, float* m, float dt, int N);
void bmhpotential(float* x,float* y, float* z, float* fx, float* fy, float* fz, float* A, float* b, float* ro, float* q,int* type ,int N,float L,float* properties);
float rast (float x1,float x2,float L);
void GenVel(float* vx, float* vy, float* vz, float* m,float temp, int N);
void SetAtom(float* m, int* Atype, float* mset,int N);
void check(float* x,float* y,float* z,float L, int N);
void swpotential(float* x,float* y, float* z, float* fx, float* fy, float* fz,float* la, float* ga, float* rc, float* cosphi, int* type, int N, float L);
void bmhtest(float* A, float* b, float* ro, float* q,int vv);
float CurrentTemp(float* vx, float* vy ,float* vz,float* m, int N);
void integrateNH(float* x, float* y,float* z,float* fx, float* fy, float* fz,
	float* vx, float* vy, float* vz, float* m, float dt, int N,float CurT, float RefT, float pnh);
void TempCouple(float* vx, float* vy, float* vz,float CurTemp,float RefTemp, int N);
void PressCouple(float* x, float* y,float* z,float CurP,float RefP,int N);
void CalcRDF(float* x, float* y, float* z,int* type,float* rdf, int rdfcount, int N,float L);
void rdfOut(float* rdf, int rdfcount,int vv,float L);
float CurrentPress(float W, float T, float ro, int N);
void ResetAtoms(int* type, int N,int toset);


int main(int argc, char* argv[]){
	if (argc<1){
	printf("Set Simulation parameter files\n");
	return 1;
	}
	else {
	printf(" Iniitia parameters file %s \n", argv[1]);
	}
	FILE* InitFile=fopen(argv[1],"r");
		fscanf(InitFile,"%d",&Latice);
		fscanf(InitFile,"%f   %d",&temp,&TFtemp);
		fscanf(InitFile,"%f   %d",&pressure,&TFpress);
		fscanf(InitFile,"%d",&PerType);
	fclose(InitFile);
	printf("Temperature coupling %d Reference temperature: %f \n", TFtemp, temp);
	printf("Pressure coupling %d Reference temperature: %f \n", TFpress, pressure);
	system("sleep 1");
	// rand initial
	rseed=rand()%32000+1;
	printf("Random rseed %d \n",rseed);
	Natom=MAXATOM;
	printf("%d\n",Natom);
	ax=(float*)calloc(Natom,sizeof(float));
	ay=(float*)calloc(Natom,sizeof(float));
	az=(float*)calloc(Natom,sizeof(float));
	vx=(float*)calloc(Natom,sizeof(float));
	vy=(float*)calloc(Natom,sizeof(float));
	vz=(float*)calloc(Natom,sizeof(float));
	fx=(float*)calloc(Natom,sizeof(float));
	fy=(float*)calloc(Natom,sizeof(float));
	fz=(float*)calloc(Natom,sizeof(float));
	AType=(int*)calloc(Natom,sizeof(int));	
	mass=(float*)calloc(Natom,sizeof(float));
	ndens=(float*)calloc(MAXATOMTYPE,sizeof(float));
	printf("allocate done\n");
	//
	TypeKol=(int*)calloc(MAXATOMTYPE,sizeof(int));
	InitMass=(float*)calloc(MAXATOMTYPE,sizeof(float));
	InitMass[1]=28.085f;	//[g/mol]
	InitMass[2]=15.999f;
	InitMass[3]=14.0f;
	A_p=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	b_p=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	ro_p=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	q_p=(float*)calloc(MAXATOMTYPE,sizeof(float));
	
	la_sw=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	ga_sw=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	rc_sw=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	cos_sw=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	
	rdf=(float*)calloc(MAXRDFBIN*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	printf( "rdf diension %d \n",MAXRDFBIN*MAXATOMTYPE*MAXATOMTYPE);
	properties=(float*)calloc(20,sizeof(float));
	system("sleep 1");
	
	
	//parameters Si-Si
	int atom1=1;
	int atom2=1;
	int atom3=1;
	A_p[atom1*MAXATOMTYPE+atom2]=11321625.0f;	//[A*g/ps^2/mol]
	A_p[atom2*MAXATOMTYPE+atom1]=11321625.0f;
	b_p[atom1*MAXATOMTYPE+atom2]=2.34f;	//[A]
	b_p[atom2*MAXATOMTYPE+atom1]=2.34f;
	ro_p[atom1*MAXATOMTYPE+atom2]=0.29f;	//[A]
	ro_p[atom2*MAXATOMTYPE+atom1]=0.29f;
	
	//parameters Si-O
	atom1=1;
	atom2=2;
	A_p[atom1*MAXATOMTYPE+atom2]=18066423.87;	//
	A_p[atom2*MAXATOMTYPE+atom1]=18066423.87f;
	b_p[atom1*MAXATOMTYPE+atom2]=2.29f;	//[A]
	b_p[atom2*MAXATOMTYPE+atom1]=2.29f;
	ro_p[atom1*MAXATOMTYPE+atom2]=0.29f;	//[A]
	ro_p[atom2*MAXATOMTYPE+atom1]=0.29f;
	//paraemters O-O
	atom1=2;
	atom2=2;
	A_p[atom1*MAXATOMTYPE+atom2]=6624355.0f;	//
	A_p[atom2*MAXATOMTYPE+atom1]=6624355.0f;
	b_p[atom1*MAXATOMTYPE+atom2]=2.34f;	//[A]
	b_p[atom2*MAXATOMTYPE+atom1]=2.34f;
	ro_p[atom1*MAXATOMTYPE+atom2]=0.29f;	//[A]
	ro_p[atom2*MAXATOMTYPE+atom1]=0.29f;
	
	//paraemters Si-C
	atom1=1;
	atom2=3;
	A_p[atom1*MAXATOMTYPE+atom2]=7940000.0f;	//
	A_p[atom2*MAXATOMTYPE+atom1]=7940000.0f;
	b_p[atom1*MAXATOMTYPE+atom2]=3.5f;	//[A]
	b_p[atom2*MAXATOMTYPE+atom1]=3.5f;
	ro_p[atom1*MAXATOMTYPE+atom2]=1.899f;	//[A]
	ro_p[atom2*MAXATOMTYPE+atom1]=1.899f;
	
	//paraemters C-C
	atom1=3;
	atom2=3;
	A_p[atom1*MAXATOMTYPE+atom2]=7940000.0f;	//
	A_p[atom2*MAXATOMTYPE+atom1]=7940000.0f;
	b_p[atom1*MAXATOMTYPE+atom2]=3.0f;	//[A]
	b_p[atom2*MAXATOMTYPE+atom1]=3.0f;
	ro_p[atom1*MAXATOMTYPE+atom2]=1.54f;	//[A]
	ro_p[atom2*MAXATOMTYPE+atom1]=1.54f;
	//charges
	
	//paraemters C-O
	atom1=2;
	atom2=3;
	A_p[atom1*MAXATOMTYPE+atom2]=7940000.0f;	//
	A_p[atom2*MAXATOMTYPE+atom1]=7940000.0f;
	b_p[atom1*MAXATOMTYPE+atom2]=3.0f;	//[A]
	b_p[atom2*MAXATOMTYPE+atom1]=3.0f;
	ro_p[atom1*MAXATOMTYPE+atom2]=0.29f;	//[A]
	ro_p[atom2*MAXATOMTYPE+atom1]=0.29f;
	//charges
	q_p[1]=4.0f;
	q_p[2]=-2.0f;
	q_p[3]=0.0f;
	
	atom1=1;	//fo Si
		atom2=2;
		atom3=2;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom2+atom3;
	la_sw[tempi]=1083985.0f;
	ga_sw[tempi]=2.6f;
	rc_sw[tempi]=3.0f;
	cos_sw[tempi]=-1.0f/3.0f;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom3+atom2;
	la_sw[tempi]=1083985.0f;
	ga_sw[tempi]=2.6f;
	rc_sw[tempi]=3.0f;
	cos_sw[tempi]=-1.0f/3.0f;
	
	atom1=2;	//fo O
		atom2=1;
		atom3=1;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom2+atom3;
	la_sw[tempi]=18066.0f;
	ga_sw[tempi]=2.0f;
	rc_sw[tempi]=2.60f;
	cos_sw[tempi]=-1.0f;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom3+atom2;
	la_sw[tempi]=18066.0f;
	ga_sw[tempi]=2.0f;
	rc_sw[tempi]=2.6f;
	cos_sw[tempi]=-1.0f;
	
	atom1=3;	//fo C-Si
		atom2=1;
		atom3=1;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom2+atom3;
	la_sw[tempi]=18006600.0f;
	ga_sw[tempi]=2.0f;
	rc_sw[tempi]=4.50f;
	cos_sw[tempi]=-0.38540f;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom3+atom2;
	la_sw[tempi]=18006600.0f;
	ga_sw[tempi]=2.0f;
	rc_sw[tempi]=4.5f;
	cos_sw[tempi]=-0.3854f;
	
	atom1=1;	//fo Si-C
		atom2=1;
		atom3=3;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom2+atom3;
	la_sw[tempi]=1806600.0f;
	ga_sw[tempi]=2.0f;
	rc_sw[tempi]=4.5f;
	cos_sw[tempi]=-0.3938f;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom3+atom2;
	la_sw[tempi]=18006600.0f;
	ga_sw[tempi]=2.0f;
	rc_sw[tempi]=4.5f;
	cos_sw[tempi]=-0.3938f;
	
	atom1=3;	//fo С-C
		atom2=3;
		atom3=3;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom2+atom3;
	la_sw[tempi]=180006600.0f;
	ga_sw[tempi]=2.0f;
	rc_sw[tempi]=4.5f;
	cos_sw[tempi]=1.0f;
	tempi=MAXATOMTYPE*MAXATOMTYPE*atom1+MAXATOMTYPE*atom3+atom2;
	la_sw[tempi]=180660000.0f;
	ga_sw[tempi]=2.0f;
	rc_sw[tempi]=4.5f;
	cos_sw[tempi]=1.0f;
	
	AtomPlace(ax,ay,az,AType,Latice,Natom);
	HBox=Latice*SILATICE;
	SystemDens=0.0f;
	for(int i=1;i<MAXATOMTYPE;i++){
		ndens[i]=TypeKol[i]/HBox/HBox/HBox;
		SystemDens+=ndens[i];
	}
	printf("Total system density %f \n", SystemDens);
	//HBox=HBox*2.0;
	SetAtom(mass,AType,InitMass,Natom);
	GenVel(vx, vy, vz,mass,temp,Natom);
	bmhtest(A_p,b_p,ro_p,q_p,3);
	SysTime=0.0;
	rdfcount=0;
		//start initial 4000 K temperature XXX ps
		temp=4000.0f;
		TFpress=0;
		for(int step=1; step<300;step++){
			SysTime+=DelT;
			SysTemp=CurrentTemp(vx,vy,vz,mass,Natom);
			SysPress=CurrentPress(properties[1],SysTemp,SystemDens,Natom);
			if(TFtemp==1){
				TempCouple(vx, vy, vz, SysTemp,temp, Natom);
			}
			if(TFpress==1){
				PressCouple(ax,ay,az,SysPress,pressure,Natom);
			}
			printf("Temprerature %f  pressure %f \n", SysTemp,SysPress);
			bmhpotential(ax,ay,az,fx,fy,fz,A_p,b_p,ro_p,q_p,AType,Natom,HBox,properties);
			swpotential(ax,ay,az,fx,fy,fz,la_sw,ga_sw,rc_sw, cos_sw,AType,Natom,HBox);
			printf("step # %d \n", step);
			integrate(ax,ay,az,fx,fy,fz,vx,vy,vz,mass,DelT,Natom);
			//integrateNH(ax,ay,az,fx,fy,fz,vx,vy,vz,mass,DelT,Natom,SysTemp,temp,pnh);
			check(ax,ay,az,HBox,Natom);
			//CalcRDF(ax, ay, az,AType,rdf,rdfcount,Natom,HBox);
			//rdfcount+=1;
			//void integrate(float* x, float* y,float* z,float* fx, float* fy, float* fz, float* vx, float* vy, float* vz, float* m, float dt, int N){
		groout("step1.gro",ax,ay,az,vx,vy,vz,AType,Natom,HBox);
		}
		ResetAtoms(AType,Natom,floor(TypeKol[2]*0.1));
		SetAtom(mass,AType,InitMass,Natom);
		// PVT at 298 K  XXX ps
		temp=298.0f;
		//TFpress=1;
		pressure=10000;
		for(int step=1; step<300;step++){
			SysTime+=DelT;
			SysTemp=CurrentTemp(vx,vy,vz,mass,Natom);
			SysPress=CurrentPress(properties[1],SysTemp,SystemDens,Natom);
			if(TFtemp==1){
				TempCouple(vx, vy, vz, SysTemp,temp, Natom);
			}
			if(TFpress==1){
				PressCouple(ax,ay,az,SysPress,pressure,Natom);
			}
			printf("Temprerature %f  pressure %f \n", SysTemp,SysPress);
			bmhpotential(ax,ay,az,fx,fy,fz,A_p,b_p,ro_p,q_p,AType,Natom,HBox,properties);
			swpotential(ax,ay,az,fx,fy,fz,la_sw,ga_sw,rc_sw, cos_sw,AType,Natom,HBox);
			printf("step # %d \n", step);
			integrate(ax,ay,az,fx,fy,fz,vx,vy,vz,mass,DelT,Natom);
			//integrateNH(ax,ay,az,fx,fy,fz,vx,vy,vz,mass,DelT,Natom,SysTemp,temp,pnh);
			check(ax,ay,az,HBox,Natom);
			//CalcRDF(ax, ay, az,AType,rdf,rdfcount,Natom,HBox);
			//rdfcount+=1;
			//void integrate(float* x, float* y,float* z,float* fx, float* fy, float* fz, float* vx, float* vy, float* vz, float* m, float dt, int N){
		}
		groout("step2.gro",ax,ay,az,vx,vy,vz,AType,Natom,HBox);
		// averagin at 298 K  XXX ps
		temp=298.0f;
		for(int step=1; step<300;step++){
			SysTime+=DelT;
			SysTemp=CurrentTemp(vx,vy,vz,mass,Natom);
			SysPress=CurrentPress(properties[1],SysTemp,SystemDens,Natom);
			if(TFtemp==1){
				TempCouple(vx, vy, vz, SysTemp,temp, Natom);
			}
			if(TFpress==1){
				PressCouple(ax,ay,az,SysPress,pressure,Natom);
			}
			printf("Temprerature %f  pressure %f \n", SysTemp,SysPress);
			bmhpotential(ax,ay,az,fx,fy,fz,A_p,b_p,ro_p,q_p,AType,Natom,HBox,properties);
			swpotential(ax,ay,az,fx,fy,fz,la_sw,ga_sw,rc_sw, cos_sw,AType,Natom,HBox);
			printf("step # %d \n", step);
			integrate(ax,ay,az,fx,fy,fz,vx,vy,vz,mass,DelT,Natom);
			//integrateNH(ax,ay,az,fx,fy,fz,vx,vy,vz,mass,DelT,Natom,SysTemp,temp,pnh);
			check(ax,ay,az,HBox,Natom);
			CalcRDF(ax, ay, az,AType,rdf,rdfcount,Natom,HBox);
			rdfcount+=1;
			//void integrate(float* x, float* y,float* z,float* fx, float* fy, float* fz, float* vx, float* vy, float* vz, float* m, float dt, int N){
		}
		//groout("file.gro",ax,ay,az,vx,vy,vz,AType,Natom,HBox);
	printf("  %d \n", Natom);
	rdfOut(rdf,rdfcount,2,HBox);
	groout("file.gro",ax,ay,az,vx,vy,vz,AType,Natom,HBox);
	free(ax);
	free(ay);
	free(az);
	free(vx);
	free(vy);
	free(vz);
	free(fx);
	free(fy);
	free(fz);
	free(AType);
	free(TypeKol);
printf("sucsesfull \n %f", erfc(1.0f));

return 1;
}

void AtomPlace(float* x, float* y, float* z, int* T, int L,int &N){
	int curn=0;
	float h=SILATICE;
	printf("%d\n",L);
	for(int i=0;i<L;i++){
		for(int j=0;j<L;j++){
			for(int k=0;k<L;k++){
				//Si insertion
				if((i+j+k)%2==0){
				ax[curn]=i*h;
				ay[curn]=j*h;
				az[curn]=k*h;
				AType[curn]=1;
				//printf("%d %d %d \n", i, j, k);
				//printf("%d %f %f %f\n", curn, ax[curn],ay[curn],az[curn] );
				curn++;
				//printf("%d \n", curn);
				}
				//center Si insertion
				ax[curn]=(i+0.5f)*h;
				ay[curn]=(j+0.5f)*h;
				az[curn]=(k+0.5f)*h;
				AType[curn]=1;
				curn++;
			}
		}
	}
	//oxigen insertion
	for(int i=0;i<L;i++){
		for(int j=0;j<L;j++){
			for(int k=0;k<L;k++){
				// O insertion
				if((i+j+k)%2==0){
				ax[curn]=(i+0.5f-0.25f)*h;
				ay[curn]=(j+0.5f-0.25f)*h;
				az[curn]=(k+0.5f-0.25f)*h;
				AType[curn]=2;
				curn++;
				ax[curn]=(i+0.5f+0.25f)*h;
				ay[curn]=(j+0.5f+0.25f)*h;
				az[curn]=(k+0.5f-0.25f)*h;
				AType[curn]=2;
				curn++;
				ax[curn]=(i+0.5f+0.25f)*h;
				ay[curn]=(j+0.5f-0.25f)*h;
				az[curn]=(k+0.5f+0.25f)*h;
				AType[curn]=2;
				curn++;
				ax[curn]=(i+0.5f-0.25f)*h;
				ay[curn]=(j+0.5f+0.25f)*h;
				az[curn]=(k+0.5f+0.25f)*h;
				AType[curn]=2;
				curn++;
				}
				else{
				//center O insertion
				ax[curn]=(i+0.5f+0.25f)*h;
				ay[curn]=(j+0.5f+0.25f)*h;
				az[curn]=(k+0.5f+0.25f)*h;
				AType[curn]=2;
				curn++;
				ax[curn]=(i+0.5f-0.25f)*h;
				ay[curn]=(j+0.5f-0.25f)*h;
				az[curn]=(k+0.5f+0.25f)*h;
				AType[curn]=2;
				curn++;
				ax[curn]=(i+0.5f-0.25f)*h;
				ay[curn]=(j+0.5f+0.25f)*h;
				az[curn]=(k+0.5f-0.25f)*h;
				AType[curn]=2;
				curn++;
				ax[curn]=(i+0.5f+0.25f)*h;
				ay[curn]=(j+0.5f-0.25f)*h;
				az[curn]=(k+0.5f-0.25f)*h;
				AType[curn]=2;
				curn;
				}
			}
		}
	}
	printf("N  %d  ", curn);
	N=curn;
	for(int i=1;i<curn;i++){
		TypeKol[AType[i]]++;
	}
	for(int i=1;i<4;i++){
		printf("Particles type %d %d  \n", i, TypeKol[i]);
	}
}

void groout(const char* filename,float* x,float* y,float* z,float* vx,float* vy,float* vz, int* T, int N, float Hbox){
	const char* naz[]={"Zero","Si","O","C"};
	FILE* fileout=fopen(filename,"w");
	fprintf(fileout,"pure silica membrane\n");
	fprintf(fileout," %d \n",N);
	for (int i=0;i<N;i++){
		//printf("%s    %s\n",  naz[1],naz[2]);
		//printf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "HybSi","Si",i,x[i],y[i],z[i],0.0,0.0,0.0);
		fprintf(fileout,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i+1, "HybSi",naz[T[i]],i+1,x[i]/10.0f,y[i]/10.0f,z[i]/10.0f,vx[i]/10.0f,vy[i]/10.0f,vz[i]/10.0f);
		//fprintf(fileout)
	}
	fprintf(fileout,"%f   %f   %f \n", Hbox,Hbox,Hbox);
	fclose(fileout);
}

void integrate(float* x, float* y,float* z,float* fx, float* fy, float* fz,
	float* vx, float* vy, float* vz, float* m, float dt, int N){
	for(int i=0;i<N;i++){
//		v[i] += f[i]*dt/m;
//		r[i] += v[i]*dt;
//		f[i] = make_float3(0.0f, 0.0f, 0.0f);
		vx[i]+=fx[i]*dt/m[i];	//[A/ps]
		vy[i]+=fy[i]*dt/m[i];
		vz[i]+=fz[i]*dt/m[i];
		//printf("i %d before %f", i,x[i]);
		x[i]+=vx[i]*dt;	//[A]
		y[i]+=vy[i]*dt;
		z[i]+=vz[i]*dt;
		//printf(" speed %f after %f\n",vx[i],x[i]);
		fx[i]=0.0;	// [A*g/ps^2/mol]
		fy[i]=0.0;
		fz[i]=0.0;
	}
}

void integrateNH(float* x, float* y,float* z,float* fx, float* fy, float* fz,
	float* vx, float* vy, float* vz, float* m, float dt, int N,float CurT, float RefT, float pnh){
	for(int i=0;i<N;i++){
//		v[i] += f[i]*dt/m;
//		r[i] += v[i]*dt;
//		f[i] = make_float3(0.0f, 0.0f, 0.0f);
		vx[i]=vx[i]*(1.0f)+fx[i]*dt/m[i]-pnh/QNH*vx[i]*dt;	//[A/ps]
		vy[i]=vy[i]*(1.0f)+fy[i]*dt/m[i]-pnh/QNH*vy[i]*dt;
		vz[i]=vz[i]*(1.0f)+fz[i]*dt/m[i]-pnh/QNH*vz[i]*dt;
		//printf("i %d before %f", i,x[i]);
		x[i]+=vx[i]*dt;	//[A]
		y[i]+=vy[i]*dt;
		z[i]+=vz[i]*dt;
		//printf(" speed %f after %f\n",vx[i],x[i]);
		fx[i]=0.0;	// [A*g/ps^2/mol]
		fy[i]=0.0;
		fz[i]=0.0;
	}
}

void bmhpotential(float* x,float* y, float* z, float* fx, float* fy, float* fz, float* A, float* b, float* ro, float* q,int* type ,int N,float L, float* properties){
	float tx;
	float ty;
	float tz;
	float rz;
	int lin;
	float fbmh;
	int i;
	int j;
	//float Wtemp;
	properties[1]=0.0f;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(i!=j){
				tx=rast(x[i],x[j],L);	//minimal obraz
				ty=rast(y[i],y[j],L);
				tz=rast(z[i],z[j],L);
				rz=sqrt(tx*tx+ty*ty+tz*tz);
				if(rz<L/2.0f){
					lin=type[i]*MAXATOMTYPE+type[j];
					if(rz<ro[lin]){
						rz=ro[lin];
					}
					fbmh=0.0;
					//printf("i %d j %d type i %d type j %d , lin %d\n", i,j,type[i],type[j],lin);
					if(((type[i]==1)&&(type[j]==2))||((type[i]==2)&&(type[j]==1))){
						fbmh= -A[lin]/ro[lin]*exp(-rz/ro[lin])-KZAR*q[type[i]]*q[type[j]]/rz/rz*erfc(rz/b[lin])-2.0f*KZAR*q[type[i]]*q[type[j]]/SQRTPI/b[lin]/rz*exp(-rz*rz/b[lin]/b[lin]);
					}
					if(((type[i]==1)&&(type[j]==3))||((type[i]==3)&&(type[j]==1))){
						if (rz<b[lin]){
							fbmh=-2.0*A[lin]/10000.0*(rz-ro[lin]);
						}
					}
					if(((type[i]==3)&&(type[j]==3))||((type[i]==3)&&(type[j]==3))){
						if (rz<b[lin]){
							fbmh=-2.0*A[lin]/10000.0*(rz-ro[lin]);
						}
					}
					if(((type[i]==2)&&(type[j]==3))||((type[i]==3)&&(type[j]==2))){
						fbmh= -A[lin]/ro[lin]*exp(-rz/ro[lin])-KZAR*q[type[i]]*q[type[j]]/rz/rz*erfc(rz/b[lin])-2.0f*KZAR*q[type[i]]*q[type[j]]/SQRTPI/b[lin]/rz*exp(-rz*rz/b[lin]/b[lin]);
					}
//					printf("i %d j %d typei %d typej %d \n", i,j, type[i],type[j]);
//					printf(" A  %f ro %f b %f  q1  %f  q2  %f\n", A[lin],ro[lin],b[lin],q[type[i]],q[type[j]]);
//					printf("rz   %f   force %f\n", rz, -fbmh);
//					system("sleep 1");
					//A[lin]*exp(-rz/ro[lin])+q[i]*q[j]*KZAR*erfc(rz/b[lin]);
					fx[i]+=fbmh*tx/rz;
					fy[i]+=fbmh*ty/rz;
					fz[i]+=fbmh*tz/rz;
					properties[1]+=fx[i]*tx+fy[i]*ty+fz[i]*tz;
				}
			}
		}
	}
	//printf("W = %f ", properties[1]);
}

void swpotential(float* x,float* y, float* z, float* fx, float* fy, float* fz,float* la, float* ga, float* rc, float* cosphi, int* type, int N, float L){
	float cosal;
	float tx1;	//by i
	float ty1;
	float tz1;
	float rz1;
	float tx2;	//by j
	float ty2;
	float tz2;
	float rz2;
	float tx3;	//by k
	float ty3;
	float tz3;
	float rz3;
	float fdr1;
	float fdr2;
	float fdr3;
	int lin;
	float pexp;
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if(i!=j){
				for(int k=j;k<N;k++){
					tx1=rast(x[i],x[j],L);
					ty1=rast(y[i],y[j],L);
					tz1=rast(z[i],z[j],L);
					rz1=sqrt(tx1*tx1+tx2*tx2+tx3*tx3);
					if(rz1<rc[i]){
						tx2=rast(x[i],x[k],L);
						ty2=rast(y[i],y[k],L);
						tz2=rast(z[i],z[k],L);
						rz2=sqrt(tx2*tx2+ty2*ty2+tz2*tz2);
						if(rz2<rc[i]){
							tx3=rast(x[j],x[k],L);
							ty3=rast(y[j],y[k],L);
							tz3=rast(z[j],z[k],L);
							rz3=sqrt(tx3*tx3+ty3*ty3+tz3*tz3);
							cosal=(rz1*rz1+rz2*rz2-rz3*rz3)/2.0f/rz1/rz2;
							lin=(MAXATOMTYPE*MAXATOMTYPE*i+MAXATOMTYPE*j+k);
							pexp=exp(ga[lin]/(rz1-rc[lin])+ga[lin]/(rz2-rc[lin]));
							//dr1
							fdr1=2.0f*la[lin]*(1.0f/rz2-cosal/rz1)*(cosal-cosphi[lin])*pexp-ga[lin]*la[lin]*(cosal-cosphi[lin])*(cosal-cosphi[lin])/(rz1-rz2)/(rz1-rz2)*pexp;
							fx[i]+=fdr1*tx1/rz1;
							fy[i]+=fdr1*ty1/rz1;
							fz[i]+=fdr1*tz1/rz1;
							//dr2
							fdr2=2.0f*la[lin]*(1.0f/rz1-cosal/rz2)*(cosal-cosphi[lin])*pexp-ga[lin]*la[lin]*(cosal-cosphi[lin])*(cosal-cosphi[lin])*pexp/(rz2-rz1)/(rz2-rz1);
							fx[i]+=fdr2*tx2/rz2;
							fy[i]+=fdr2*ty2/rz2;
							fz[i]+=fdr2*tz2/rz2;
							//dr3
							fdr3=3.0f*la[lin]* rz3*rz3*(cosal-cosphi[lin])*pexp/rz1/rz2;
							fx[i]+=fdr3*tx3/rz3;
							fy[i]+=fdr3*ty3/rz3;
							fz[i]+=fdr3*tz3/rz3;
						}
					}
					
				}
			}
		}
	}
}

float rast (float x1,float x2,float L){
	float ftemp =x2-x1;
	if(ftemp>L/2.0f){
		ftemp-=L;
	}
	if(ftemp<-L/2.0f){
		ftemp+=L;
	}
	return ftemp;
}

void GenVel(float* vx, float* vy, float* vz, float* m,float temp, int N){
	float AvVel;
	for(int i=0;i<N;i++){
		AvVel=sqrt(kB*temp/m[i]/10.0f);	//[A/ps]
		//printf("i %d   AvVel %f \n", i,AvVel );
		vx[i]=AvVel*gasdev(&rseed);
		vy[i]=AvVel*gasdev(&rseed);
		vz[i]=AvVel*gasdev(&rseed);
	}
}

void SetAtom(float* m, int* Atype, float* mset,int N){
	for(int i=0;i<N;i++){
		m[i]=mset[Atype[i]];
	}

}

void check(float* x,float* y,float* z,float L, int N){
	for(int i=0;i<N;i++){
		if(x[i]<0.0f){
			x[i]+=L;
		}
		if(x[i]>L){
			x[i]-=L;
		}
		if(y[i]<0.0f){
			y[i]+=L;
		}
		if(y[i]>L){
			y[i]-=L;
		}
		if(z[i]<0.0){
			z[i]+=L;
		}
		if(z[i]>L){
			z[i]-=L;
		}
	}
}
void bmhtest(float* A, float* b, float* ro, float* q,int vv){
	FILE* TestFile=fopen("bmh.out","w");
	fprintf(TestFile," numbers of atoms %d \n",vv);
	fprintf(TestFile," parameter %d ",vv);
	int i;
	int j;
	int k;
	float rs;
	float fbmh;
	int lin;
	float pot;
		for(i=0;i<vv;i++){
			for(j=0;j<vv;j++){
				fprintf(TestFile," %1d-%1d ",i+1,j+1);
			}
		}
		fprintf(TestFile," \n");
		fprintf(TestFile,"A [] ");	//A
		for(i=0;i<vv;i++){
			for(j=0;j<vv;j++){
				fprintf(TestFile," %f ", A[MAXATOMTYPE*(i+1)+(j+1)]);
			}
		}
		fprintf(TestFile," \n");
		fprintf(TestFile,"ro [] ");	//ro
		for(i=0;i<vv;i++){
			for(j=0;j<vv;j++){
				fprintf(TestFile," %f ", ro[MAXATOMTYPE*(i+1)+(j+1)]);
			}
		}
		fprintf(TestFile," \n");
		fprintf(TestFile,"b [] ");	//B
		for(i=0;i<vv;i++){
			for(j=0;j<vv;j++){
				fprintf(TestFile," %f ", b[MAXATOMTYPE*(i+1)+(j+1)]);
			}
		}
		fprintf(TestFile," \n");
		fprintf(TestFile," \n");
		fprintf(TestFile,"z [] ");	//B
		for(i=0;i<vv;i++){
			fprintf(TestFile," %f ", q[i+1]);
		}
		fprintf(TestFile," \n \n");
				
		for(k=1;k<1000;k++){
			rs=20.0f/1000*k;
			fprintf(TestFile," %f ", rs);
			for(i=0;i<vv;i++){
				for(j=0;j<vv;j++){
				
					lin=MAXATOMTYPE*(i+1)+(j+1);
					//printf(" %d %d %f %f %f %f \n",i+1,j+1,A[lin], b[lin], ro[lin], q[i+1]);
					//system("sleep 4");
					pot=A[lin]*exp(-rs/ro[lin])+KZAR*q[i+1]*q[j+1]/rs*erfc(rs/b[lin]);
					fbmh=-(-A[lin]/ro[lin]*exp(-rs/ro[lin])-KZAR*q[i+1]*q[j+1]/rs/rs*erfc(rs/b[lin])-2.0f*KZAR*q[i+1]*q[j+1]/SQRTPI/b[lin]/rs*exp(-rs*rs/b[lin]/b[lin])); //-KZAR*q[i+1]*q[j]+1/rs/rs*erfc(rs/b[lin]); //-2.0f*KZAR*q[i]*q[j]/SQRTPI/b[lin]/rs*exp(-rs*rs/b[lin]/b[lin]); //
					fprintf(TestFile," %f  %f ", pot, fbmh);
				}
			}
			fprintf(TestFile," \n");
		}
		
	fclose(TestFile);
}

float CurrentTemp(float* vx, float* vy ,float* vz,float* m, int N){
	float Vel;
	int i;
	float sum;
	Vel=0.0f;
	sum=0.0f;
	for(i=0;i<N;i++){
		Vel=0.0f;
		Vel+=vx[i]*vx[i];	//[A^2/ps^2]
		Vel+=vy[i]*vy[i];
		Vel+=vz[i]*vz[i];
		sum+=m[i]*Vel;	//[A^2*g/ps^2/mol]
	}
	return sum/N/kB/3.0f*10.0f;	//[K]
}

float CurrentPress(float W, float T, float ro, int N){
	//printf("w %f t %f ro %f n %d \n", W,T,ro, N);
	return ro*T*kB*16.611+W/N*ro/3.0f*166.11;	//
}

void TempCouple(float* vx, float* vy, float* vz,float CurTemp,float RefTemp, int N){
	int i;
	float la;
	la=sqrt(1.0f+0.9f*(CurTemp/RefTemp-1.0f));
	//printf("la = %f \n", la);
	for(i=1;i<N;i++){
			vx[i]=vx[i]/la;
			vy[i]=vy[i]/la;
			vz[i]=vz[i]/la;
		}		
}

void PressCouple(float* x,float* y, float* z, float CurP,float RefP, int N){
	float mu;
	int i;
	mu=1.0f+5.0f/300000000000.0f*(CurP-RefP);
	if(mu>1.1f){
		mu=1.1f;
	}
	if(mu<0.9f){
		mu=0.9f;
	}
	printf("mu %f current pressure %f  reference pressure % f\n", mu, CurP, RefP);
	for(i=0;i<N;i++){
		x[i]*=mu;
		y[i]*=mu;
		z[i]*=mu;
	}
	HBox=HBox*mu;
}

void CalcRDF(float* x, float* y, float* z,int* type,float* rdf, int rdfcount, int N,float L){
	int i;
	int j;
	int lin;
	float tx;
	float ty;
	float tz;
	float rz;
	int hist;
	for(i=0;i<N;i++){
		//printf("i = %d  ", i);
		for(j=i+1;j<N;j++){
			//printf(" j = %d  ", j);
			tx=rast(x[i],x[j],L);	//minimal obraz
			ty=rast(y[i],y[j],L);
			tz=rast(z[i],z[j],L);
			rz=sqrt(tx*tx+ty*ty+tz*tz);
			//printf("rz %f  ", rz );
			if(rz<L){
				hist=floor(rz/RDFBIN);
				//printf("hist %d", hist);
				lin=MAXATOMTYPE*MAXRDFBIN*type[i]+MAXRDFBIN*type[j]+hist;
				//printf("hist %d  lin  %d type[i] %d type[j] %d rdf %f \n", hist, lin,type[i],type[j], rdf[lin]);
				rdf[lin]+=2.0f;
			}
			
		}
	}
	//system("sleep 2");
}

void rdfOut(float* rdf, int rdfcount,int vv,float L){
	int i;
	int j;
	int k;
	int lin;
	FILE* fileout=fopen("rdf.out","w");
	//fprintf(fileout,"");
	fprintf(fileout," r [A], ");
	for(i=1;i<vv+1;i++){
		for(j=1;j<vv+1;j++){
			fprintf(fileout," %d ", i);
		}
	}
	fprintf(fileout,"  \n");
	for(k=0;k<MAXRDFBIN;k++){
		if(k*RDFBIN<L){
			fprintf(fileout," %f ", k*RDFBIN);
			for(i=1;i<vv+1;i++){
				for(j=1;j<vv+1;j++){
					lin=MAXATOMTYPE*MAXRDFBIN*i+MAXRDFBIN*j+k;
					fprintf(fileout," %d %f ",lin, rdf[lin]/rdfcount);
				}
			}
		fprintf(fileout," \n");
		}
	}
	fclose(fileout);
	//printf("rdfcount %d \n", rdfcount);
}

void ResetAtoms(int* type, int N,int toset){
	int RandAtom;
	int tempn;
	int set;
	int i;
	float frand;
	tempn=N;
	set=0;
	printf("reset start N %d set %d toset %d\n", N, set, toset);
	while(set<toset){
		frand=ran2(&rseed);
		
		RandAtom=floor(frand*N);
		printf("RandAtom %d %f \n", RandAtom,frand);
		if(type[RandAtom]==2){
			type[RandAtom]=3;	//chenge O to C
			type[tempn]=3;	//add new C
			ax[tempn]=ax[RandAtom];	//set initial coordinates
			ay[tempn]=ay[RandAtom];
			az[tempn]=az[RandAtom]+0.8;
			printf("ResetAtom %d, tempn %d, set %d \n", RandAtom, tempn, set);
			tempn+=1;
			set+=1;
		}
	}
	Natom=tempn;
}
