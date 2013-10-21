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
#define MAXNEAR  30
#define PI 3.14159265359f


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

int* bond;
int* angle1;
int* angle2;
int* tors1;
int* tors2;
int* tors3;

float* PBond1;	//bonds parameters
float* PBond2;
float* PBond3;
float* PBond4;

float* PAngle1;	//angle parameters
float* PAngle2;
float* PAngle3;
float* PAngle4;
float* fi0;

float* VDWA;
float* VDWB;
float* CuQ;

float percent;



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
void GetBonds(int N,float L,float br);
void GetAngles(int N);
void GetDih(int N);
void MixingParameters();
void PotAngle(int N, int L);
void SummForce(int N, float L);
void itpout(int N,const char* filename);
void Replace(int N, float proc);

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
		fscanf(InitFile,"%f",&percent);
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
	bond=(int*)calloc(MAXATOM*MAXNEAR+MAXATOM,sizeof(int));
	angle1=(int*)calloc(MAXATOM*MAXNEAR+MAXATOM,sizeof(int));
	angle2=(int*)calloc(MAXATOM*MAXNEAR+MAXATOM,sizeof(int));
	tors1=(int*)calloc(MAXATOM*MAXNEAR*MAXNEAR+MAXATOM,sizeof(int));
	tors2=(int*)calloc(MAXATOM*MAXNEAR*MAXNEAR+MAXATOM,sizeof(int));
	tors3=(int*)calloc(MAXATOM*MAXNEAR*MAXNEAR+MAXATOM,sizeof(int));
	
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
	//system("sleep 1");
	
	//
	PBond1=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	PBond2=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	PBond3=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	PBond4=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	
	PAngle1=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	PAngle2=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	PAngle3=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	PAngle4=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	fi0=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	
	VDWA=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	VDWB=(float*)calloc(MAXATOMTYPE*MAXATOMTYPE,sizeof(float));
	CuQ=(float*)calloc(MAXATOMTYPE,sizeof(float));
	
	MixingParameters();
	AtomPlace(ax,ay,az,AType,Latice,Natom);
	HBox=Latice*SILATICE;
	SystemDens=0.0f;
	for(int i=1;i<MAXATOMTYPE;i++){
		ndens[i]=TypeKol[i]/HBox/HBox/HBox;
		SystemDens+=ndens[i];
		printf(" TypeKol[ %d ] = %d  \n", i, TypeKol[i]);
	}
	printf("Total system density %f \n", SystemDens);
	//abort;
	//return 1;
	//HBox=HBox*2.0;
	SetAtom(mass,AType,InitMass,Natom);
	GenVel(vx, vy, vz,mass,temp,Natom);
	GetBonds(Natom,HBox,1.7f);
	GetAngles(Natom);
	GetDih(Natom);
		//
	//
	groout("test1.gro",ax,ay,az,vx,vy,vz,AType,Natom,HBox);
	itpout(Natom,"test1.top");
//	
//	//bmhtest(A_p,b_p,ro_p,q_p,3);
//	SysTime=0.0;
//	rdfcount=0;
//		//start initial 4000 K temperature XXX ps
//		temp=298.0f;
//		TFtemp=1;
//		TFpress=1;
//		pressure=100.0f;
//		for(int step=1; step<1;step++){
//			SysTime+=DelT;
//			SysTemp=CurrentTemp(vx,vy,vz,mass,Natom);
//			SysPress=CurrentPress(properties[1],SysTemp,SystemDens,Natom);
//			if(TFtemp==1){
//				TempCouple(vx, vy, vz, SysTemp,temp, Natom);
//			}
//			if(TFpress==1){
//				PressCouple(ax,ay,az,SysPress,pressure,Natom);
//			}
//			printf("Temprerature %f  pressure %f \n", SysTemp,SysPress);
//			SummForce(Natom,HBox);
//			printf("step # %d \n", step);
//			integrate(ax,ay,az,fx,fy,fz,vx,vy,vz,mass,DelT,Natom);
//			check(ax,ay,az,HBox,Natom);
//		//groout("step1.gro",ax,ay,az,vx,vy,vz,AType,Natom,HBox);
//		}
//	rdfOut(rdf,rdfcount,2,HBox);
//	groout("file.gro",ax,ay,az,vx,vy,vz,AType,Natom,HBox);
	
	//
	Replace(Natom, percent);
	check(ax,ay,az,HBox,Natom);
	SetAtom(mass,AType,InitMass,Natom);
	system("sleep 1");
//	for(int mi=0;mi<Natom;mi++){
//		printf("Atom # %d  type %d \n", mi, AType[mi] );
//		for(int mj=0;mj<bond[mi];mj++){
//			printf("               %d  type %d \n",bond[mi*MAXNEAR+mj+MAXATOM],AType[bond[mi*MAXNEAR+mj+MAXATOM]]);
//		}
//	}
	printf("Natom %d\n", Natom);
	system("sleep 2");
	GetAngles(Natom);
	GetDih(Natom);
	groout("test.gro",ax,ay,az,vx,vy,vz,AType,Natom,HBox);
	itpout(Natom,"test.top");
	
	for(int tei=0;tei<Natom;tei++){
		printf("Atom %d Type %d bond %d angle %d tors %d \n",tei,AType[tei],bond[tei],angle1[tei],tors1[tei]);
		for(int tej=0;tej<bond[tei];tej++){
			printf("        bn %d , an %d , at %d \n", tej, bond[tei*MAXNEAR+tej+MAXATOM],AType[bond[tei*MAXNEAR+tej+MAXATOM]]);
		}
		for(int tek=0;tek<angle1[tei];tek++){
			printf("        an %d , a1 %d , at1 %d, a2 %d, at2 %d \n", tek, angle1[tei*MAXNEAR+tek+MAXATOM],AType[angle1[tei*MAXNEAR+tek+MAXATOM]],angle2[tei*MAXNEAR+tek+MAXATOM],AType[angle2[tei*MAXNEAR+tek+MAXATOM]]);
		}
		for(int tej=0;tej<tors1[tei];tej++){
			printf("        tn %d , a1 %d , at1 %d, a2 %d, at2 %d, a3 %d, at3 %d  \n", tej, tors1[tei*MAXNEAR+tej+MAXATOM],AType[tors1[tei*MAXNEAR+tej+MAXATOM]],tors2[tei*MAXNEAR+tej+MAXATOM],AType[tors2[tei*MAXNEAR+tej+MAXATOM]],tors3[tei*MAXNEAR+tej+MAXATOM],AType[tors3[tei*MAXNEAR+tej+MAXATOM]]);
		}
	}
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
					ax[curn]=i*h;	//в узлах решетки
					ay[curn]=j*h;
					az[curn]=k*h;
					AType[curn]=1;
					//printf("%d %d %d \n", i, j, k);
					//printf("%d %f %f %f\n", curn, ax[curn],ay[curn],az[curn] );
					curn++;
				//printf("%d \n", curn);
				}
				//center Si insertion
				if((i+j+k)%2==0){
					ax[curn]=(i+0.5f)*h;	//в центре решетки
					ay[curn]=(j+0.5f)*h;
					az[curn]=(k+0.5f)*h;
					AType[curn]=1;
					curn++;
				}
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
//				ax[curn]=(i+0.5f+0.25f)*h;
//				ay[curn]=(j+0.5f+0.25f)*h;
//				az[curn]=(k+0.5f+0.25f)*h;
//				AType[curn]=2;
//				curn++;
//				ax[curn]=(i+0.5f-0.25f)*h;
//				ay[curn]=(j+0.5f-0.25f)*h;
//				az[curn]=(k+0.5f+0.25f)*h;
//				AType[curn]=2;
//				curn++;
//				ax[curn]=(i+0.5f-0.25f)*h;
//				ay[curn]=(j+0.5f+0.25f)*h;
//				az[curn]=(k+0.5f-0.25f)*h;
//				AType[curn]=2;
//				curn++;
//				ax[curn]=(i+0.5f+0.25f)*h;
//				ay[curn]=(j+0.5f-0.25f)*h;
//				az[curn]=(k+0.5f-0.25f)*h;
//				AType[curn]=2;
//				curn++;
				}
			}
		}
	}
	printf("N  %d  ", curn);
	N=curn;
	for(int i=0;i<curn;i++){
		TypeKol[AType[i]]++;
	}
	for(int i=1;i<4;i++){
		printf("Particles type %d %d  \n", i, TypeKol[i]);
	}
}

void groout(const char* filename,float* x,float* y,float* z,float* vx,float* vy,float* vz, int* T, int N, float Hbox){
	const char* naz[]={"Zero","Si","O","C2"};
	FILE* fileout=fopen(filename,"w");
	fprintf(fileout,"pure silica membrane\n");
	fprintf(fileout," %d \n",N);
	for (int i=0;i<N;i++){
		//printf("%s    %s\n",  naz[1],naz[2]);
		//printf("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i, "HybSi","Si",i,x[i],y[i],z[i],0.0,0.0,0.0);
		fprintf(fileout,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", i+1, "HybSi",naz[T[i]],i+1,x[i]/10.0f,y[i]/10.0f,z[i]/10.0f,vx[i]/10.0f,vy[i]/10.0f,vz[i]/10.0f);
		//fprintf(fileout)
	}
	fprintf(fileout,"%f   %f   %f \n", Hbox/10.0,Hbox/10.0,Hbox/10.0);
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
	//printf("reset start N %d set %d toset %d\n", N, set, toset);
	while(set<toset){
		frand=ran2(&rseed);
		
		RandAtom=floor(frand*N);
		//printf("RandAtom %d %f \n", RandAtom,frand);
		if(type[RandAtom]==2){
			type[RandAtom]=3;	//chenge O to C
			type[tempn]=3;	//add new C
			ax[tempn]=ax[RandAtom];	//set initial coordinates
			ay[tempn]=ay[RandAtom];
			az[tempn]=az[RandAtom]+0.8;
			//printf("ResetAtom %d, tempn %d, set %d \n", RandAtom, tempn, set);
			tempn+=1;
			set+=1;
		}
	}
	Natom=tempn;
}

void GetBonds(int N,float L,float br){
	int i;
	int j;
	float tx;
	float ty;
	float tz;
	float rz;
	int tempi;
	int lin;
	for(i=0;i<MAXNEAR*MAXATOM+MAXATOM;i++){
		bond[i]=0;
	}
	
	for(i=0;i<N;i++){
		tempi=0;
		for(j=0;j<N;j++){
			if(i!=j){
				tx=rast(ax[i],ax[j],L);
				ty=rast(ay[i],ay[j],L);
				tz=rast(az[i],az[j],L);
				rz=sqrt(tx*tx+ty*ty+tz*tz);
				if (rz<br){
					//printf("i %d j %d, rz %f br %f tempi %d\n", i,j,rz,br,tempi);
					//printf("ti %d tj %d \n", AType[i],AType[j]);
					lin=i*MAXNEAR+tempi+MAXATOM;
					bond[lin]=j;
					bond[i]+=1;
					tempi+=1;
				}
			}
		}
		//system("sleep 2");
		//printf("atom %d have %d bonds\n", i, bond[i]);
	}
	//system("sleep 10");
	printf("Geting Bonds DONE\n");
}
void GetAngles(int N){
	int i;
	int j;
	int k;
	int lin;
	int lin1;
	int lin2;
	int tempi;
	for(i=0;i<N;i++){
		angle1[i]=0;
		angle2[i]=0;
	}
	
	for(i=0;i<N;i++){
		tempi=0;
		for(j=0;j<bond[i];j++){
			for(k=j+1;k<bond[i];k++){
				lin=i*MAXNEAR+tempi+MAXATOM;
				lin1=i*MAXNEAR+j+MAXATOM;
				lin2=i*MAXNEAR+k+MAXATOM;
				angle1[i]++;
				angle2[i]++;
				angle1[lin]=bond[lin1];
				angle2[lin]=bond[lin2];
				//printf("      i %d tempi %d angle1 %d angle2 %d \n",i, tempi,angle1[lin],angle2[lin]);
				tempi++;
			}
		}
		printf("i %d angle1 %d angle2 %d \n",i,angle1[i],angle2[i]);
	}
	//system("sleep 1");
	printf("GEting Angles DONE \n");
}

void GetDih(int N){
	int i;
	int j;
	int k;
	int l;
	int lin;
	int tempi;
	
	int lin1;
	int lin2;
	int lin3;
	for(i=0;i<N;i++){
		tors1[i]=0;
		tors2[i]=0;
		tors3[i]=0;
	}
	
	for(i=0;i<N;i++){
		//printf("atom number # %d bonds %d \n",i,bond[i]);
		tempi=0;
		for(j=0;j<bond[i];j++){	//bond number
			lin1=i*MAXNEAR+j+MAXATOM;	//bond massive number | bond[lin1]  -atom number
			//printf("---Bonded atom 1 %d \n", bond[lin1]);
			if(bond[lin1]!=i){
				for(k=0;k<bond[bond[lin1]];k++){
					lin2=bond[lin1]*MAXNEAR+k+MAXATOM;
					//printf("--- ---Bonded atom 2 %d \n", bond[lin2]);
					if(bond[lin2]!=i){
						for(l=0;l<bond[bond[lin2]];l++){
							lin3=bond[lin2]*MAXNEAR+l+MAXATOM;
							//printf("--- --- ---Bonded atom 3 %d \n", bond[lin3]);
							if((bond[lin3]!=bond[lin1])||(bond[lin3]!=i)){
								tors1[i]++;
								tors2[i]++;
								tors3[i]++;
								lin=i*MAXNEAR+tempi+MAXATOM;
								tors1[lin]=bond[lin1];
								tors2[lin]=bond[lin2];
								tors3[lin]=bond[lin3];
								//printf("ni %d a1 %d a2 %d a3 %d at1 %d at2 %d at3 %d\n", tempi,tors1[lin],tors2[lin],tors3[lin1],AType[tors1[lin]],AType[tors2[lin]],AType[tors3[lin]]);
								tempi++;
							}
						}
					}
				}
			}
		}
		//system("sleep 20");
	}
	printf("Getting Dihedral angles DONE \n");
}

void PotToZero(int N){
	int i;
	for(i=0;i<N;i++){
		fx[i]=0.0;
		fy[i]=0.0;
		fz[i]=0.0;
	}
}

void PotBond(int N,float L){
	int i;
	int j;
	int k;
	int lin;
	float tx;
	float ty;
	float tz;
	float rz;
	float force;
	for(i=0;i<N;i++){
		for(j=0;j<bond[i];j++){
			lin=MAXNEAR*i+j+MAXATOM;
			//printf("i %d bond[lin] %d lin %d",i,bond[lin],lin);
			tx=rast(ax[i],ax[bond[lin]],L);	//minimal obraz
			ty=rast(ay[i],ay[bond[lin]],L);
			tz=rast(az[i],az[bond[lin]],L);
			rz=sqrt(tx*tx+ty*ty+tz*tz);
			//remap
			//printf("ATYpe %d \n", AType[bond[lin]]);
			lin=AType[i]*MAXATOMTYPE+AType[bond[lin]];
			
			force=2.0f*PBond2[lin]*(rz-PBond1[lin])+3.0f*PBond3[lin]*(rz-PBond1[lin])*(rz-PBond1[lin])+4.0f*PBond4[lin]*(rz-PBond1[lin])*(rz-PBond1[lin])*(rz-PBond1[lin]);
//			printf("ti %d tj %d b0 %f k2 %f k3 %f k4 %f \n", AType[i],AType[bond[MAXNEAR*i+j+MAXATOM]], PBond1[lin], PBond2[lin],PBond3[lin],PBond4[lin]);
//			printf("force %f rz %f b0 \n",force,rz);
//			system("sleep 1");
			fx[i]+=force*tx/rz;
			fy[i]+=force*ty/rz;
			fz[i]+=force*tz/rz;
			//printf("i %d force %f\n", i,force);
		}
	}
	//system("sleep 10");
}

void PotVDW(int N,float L){
	int i;
	int j;
	int lin;
	float tx;
	float ty;
	float tz;
	float rz;
	float rz2;
	float rz7;
	float rz10;
	float force;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(i!=j){
				tx=rast(ax[i],ax[j],L);	//minimal obraz
				ty=rast(ay[i],ay[j],L);
				tz=rast(az[i],az[j],L);
				rz=sqrt(tx*tx+ty*ty+tz*tz);
				if(rz<L/2.0f){
					lin=AType[i]*MAXATOMTYPE+AType[j];
					rz2=rz*rz;
					rz7=rz2*rz2*rz2*rz;
					rz10=rz7*rz2*rz;
					force+=9.0f*VDWA[lin]/rz10-6.0f*VDWB[lin]/rz7;	//Lennara-Jones
					force+=CuQ[AType[i]]*CuQ[AType[j]]/rz2;	//Coloumb
					//
					fx[i]+=force*tx/rz;
					fy[i]+=force*ty/rz;
					fz[i]+=force*tz/rz;
				}
			}
		}
	}
}

void SummForce(int N, float L){
	printf("start summ\n");
	PotToZero(N);
	printf("set force t zero DONE \n");
	PotBond(N,L);
	printf("bond force calculating DONE\n");
	//PotVDW(N,L);
	//printf("Lennard-JOnes calculating DONE\n");
	//PotAngle(N,L);
	//printf("Angle calculating DONE\n");
}

void PotAngle(int N, int L){
	int i;
	int j;
	int k;
	int lin1;
	int lin2;
	float tx1;
	float tx2;
	float tx3;
	float ty1;
	float ty2;
	float ty3;
	float tz1;
	float tz2;
	float tz3;
	float rz1;
	float rz2;
	float rz3;
	float cosfi;
	float sinfi;
	float force;
	int lin3;
	float fi;
	for(i=0;i<N;i++){
		//printf("----------------------------------------------\ni %d angle1 %d angle2 %d\n",i,angle1[i],angle2[i]);
		for(j=0;j<angle1[i];j++){
			lin1=i*MAXNEAR+j+MAXATOM;
			//printf("i %d angle1[lin1] %d\n", j, angle1[lin1]);
			tx1=rast(ax[i],ax[angle1[lin1]],L);	//vector i-j
			ty1=rast(ay[i],ay[angle1[lin1]],L);
			tz1=rast(az[i],az[angle1[lin1]],L);
			rz1=sqrt(tx1*tx1+ty1*ty1+tz1*tz1);
			//printf("tx %f ty %f tz %f rz %f\n",tx1,ty1,tz1,rz1);
			//system("sleep 10");
			for(k=0;k<angle2[i];k++){
				lin2=i*MAXNEAR+k+MAXATOM;
				//printf("k %d lin2 %d \n",k,lin2);
				tx2=rast(ax[i],ax[angle2[lin2]],L);	//vector i-k
				ty2=rast(ay[i],ay[angle2[lin2]],L);
				tz2=rast(az[i],az[angle2[lin2]],L);
				rz2=sqrt(tx2*tx2+ty2*ty2+tz2*tz2);
				//tx3=rast(ax[j],ax[k],L);
				//ty3=rast(ax[])
				//printf("rz1 %f rz2 %f \n", rz1,rz2);
				cosfi=(tx1*tx2+ty1*ty2+tz1*tz2)/rz1/rz2;
				if(cosfi>0.999){
					cosfi=0.999;
				}
				else if(cosfi<-0.999){
					cosfi=-0.999;
				}
				sinfi=sqrt(1.0f-cosfi*cosfi);
				fi=acos(cosfi); //*180.0/PI;
			//printf("type i %d type j %d type k %d \n", AType[i],AType[angle1[lin1]],AType[angle2[lin2]]);
			printf("cosfi %f sinfi %f fi %f\n",cosfi,sinfi,fi);
				// radianes
				lin3=AType[i]*MAXATOMTYPE*MAXATOMTYPE+MAXATOMTYPE*AType[angle1[lin1]]+AType[angle2[lin2]];
			printf("phi %f  h2  %f  h3  %f  h4  %f ",PAngle1[lin3],PAngle2[lin3],PAngle3[lin3],PAngle4[lin3]);
				//f i
				//part 1 by ij
				force=-(PAngle2[lin3]*(fi-fi0[lin3])+PAngle3[lin3]*(fi-fi0[lin3])*(fi-fi0[lin3])+PAngle4[lin3]*(fi-fi0[lin3])*(fi-fi0[lin3])*(fi-fi0[lin3]))*(cosfi/sinfi/rz1);
				fx[j]+=force*tx1/rz1;
				fy[j]+=force*ty1/rz1;
				fz[j]+=force*tz1/rz1;
				
				fx[i]-=force*tx1/rz1;	//for i
				fy[i]-=force*ty1/rz1;
				fz[i]-=force*tz1/rz1;
			printf("force 11 %f ",force);
				//part 2 by kj
				force=-(PAngle2[lin3]*(fi-fi0[lin3])+PAngle3[lin3]*(fi-fi0[lin3])*(fi-fi0[lin3])+PAngle4[lin3]*(fi-fi0[lin3])*(fi-fi0[lin3])*(fi-fi0[lin3]))*(-1.0f/sinfi/rz1);
				fx[j]+=force*tx2/rz2;
				fy[j]+=force*ty2/rz2;
				fz[j]+=force*tz2/rz2;
				
				fx[i]-=force*tx2/rz2;
				fy[i]-=force*ty2/rz2;
				fz[i]-=force*tz2/rz2;
			printf("force 12 %f ",force);
				//f k
				//part 1 by ij
				force=-(PAngle2[lin3]*(fi-fi0[lin3])+PAngle3[lin3]*(fi-fi0[lin3])+PAngle4[lin3]*(fi-fi0[lin3])*(fi-fi0[lin3])*(fi-fi0[lin3])*(fi-fi0[lin3]))*(cosfi/sinfi/rz2);
				fx[k]+=force*tx1/rz1;
				fy[k]+=force*ty1/rz1;
				fz[k]+=force*tz1/rz1;
				
				fx[i]-=force*tx1/rz1;	//for i
				fy[i]-=force*ty1/rz1;
				fz[i]-=force*tz1/rz1;
			printf("force 21 %f ",force);
				//part 2 by kj
				force=-(PAngle2[lin3]*(fi-fi0[lin3])+PAngle3[lin3]*(fi-fi0[lin3])+PAngle4[lin3]*(fi-fi0[lin3])*(fi-fi0[lin3])*(fi-fi0[lin3])*(fi-fi0[lin3]))*(-1.0f/sinfi/rz2);
				fx[k]+=force*tx2/rz2;
				fy[k]+=force*ty2/rz2;
				fz[k]+=force*tz2/rz2;
				
				fx[i]-=force*tx2/rz2;
				fy[i]-=force*ty2/rz2;
				fz[i]-=force*tz2/rz2;
			printf("force 22 %f \n \n",force);
			if(AType[i]==2){
				system("sleep 3");
			}
			}
		}
	}
}

void MixingParameters(){
	int i;
	int j;
	int k;
	int lin1;
	int lin2;
	int lin3;
	//
	i=1;	//bond SI-O
	j=2;
	lin1=MAXATOMTYPE*i+j;
	lin2=MAXATOMTYPE*j+i;
	PBond1[lin1]=1.666;	// b0
	PBond2[lin2]=PBond1[lin1];
	PBond2[lin1]=106226.0f;	//k2
	PBond2[lin2]=PBond2[lin1];
	PBond3[lin1]=-307920.0f;
	PBond3[lin2]=PBond3[lin2];
	PBond4[lin1]=474813.0f;
	PBond4[lin2]=PBond4[lin1];
	
	i=1;	//Angle O-Si-O
	j=2;
	k=2;
	lin1=i*MAXATOMTYPE*MAXATOMTYPE+MAXATOMTYPE*j+k;
	lin2=i*MAXATOMTYPE*MAXATOMTYPE+MAXATOMTYPE*k+j;
	PAngle1[lin1]=110.612f/180.0f*PI;
	PAngle1[lin2]=PAngle1[lin1];
	PAngle2[lin1]=8876.0f;
	PAngle2[lin2]=PAngle2[lin1];
	PAngle3[lin1]=-3951.0f;
	PAngle3[lin2]=PAngle3[lin1];
	PAngle4[lin1]=1358.0f;
	PAngle4[lin2]=PAngle4[lin1];
	
	i=2;	//Angle Si-O-Si
	j=1;
	k=1;
	lin1=i*MAXATOMTYPE*MAXATOMTYPE+MAXATOMTYPE*j+k;
	lin2=i*MAXATOMTYPE*MAXATOMTYPE+MAXATOMTYPE*k+j;
	PAngle1[lin1]=176.265f/180.0f*PI;
	PAngle1[lin2]=PAngle1[lin1];
	PAngle2[lin1]=18358.0f;
	PAngle2[lin2]=PAngle2[lin1];
	PAngle3[lin1]=37054.0f;
	PAngle3[lin2]=PAngle3[lin1];
	PAngle4[lin1]=41783.0f;
	PAngle4[lin2]=PAngle4[lin1];
}

void itpout(int N, const char* filename){
	const char* naz[]={"Zero","Si","O","C2"};
	const char* anaz[]={"Zero","Si","O","C2"};
	int i;
	int j;
	int k;
	int lin1;
	int lin2;
	int tempi;
	FILE* fileitp=fopen(filename,"w");
		fprintf(fileitp,"[ defaults ]\n");
		fprintf(fileitp,"1  3  yes  0.0  0.0\n\n");
		fprintf(fileitp,"[ atomtypes ]\n");
		fprintf(fileitp,"Si     28.0000    0.0000 A     0.3000000000E+00    0.00000000E-00\n");
		fprintf(fileitp,"O     16.0000    0.0000 A     0.3000000000E+00    0.00000000E-00\n");
		fprintf(fileitp,"C2     14.0000    0.0000 A     0.3000000000E+00    0.00000000E-00\n\n");
		fprintf(fileitp,"[ bondtypes ]\n");
		fprintf(fileitp,"  Si   O     1    0.1666  120000.0\n");
		fprintf(fileitp,"  Si   C2    1    0.1899  106000.0\n");
		fprintf(fileitp,"  C2   C2    1    0.1899   40580.0\n\n");
		fprintf(fileitp,"[ angletypes ]\n");
//		fprintf(fileitp," O  Si  O  6 110.612  0.0  0.0   88.76  -39.53    13.58  \n");
//		fprintf(fileitp," Si  O  Si  6 176.265  0.0  0.0  183.58  370.54   417.83  \n");
//		fprintf(fileitp," C2  Si  C2  6 113.19  0.0  0.0  151.72  -85.43   83.88  \n");
//		fprintf(fileitp," O  Si  C2 6 110.612  0.0  0.0   88.76  -39.53    13.58  \n");
//		fprintf(fileitp," Si  C2  C2  6 112.67  0.0  0.0  165.58  -31.17   0.0  \n\n");
		fprintf(fileitp," O  Si  O  1 110.612  472.0  \n");
		fprintf(fileitp," O  Si  C2  1 110.612  472.0  \n");
		fprintf(fileitp," Si  O  Si  1 176.265  8794.0  \n");
		fprintf(fileitp," C2  Si  C2  1 113.19  2216.0  \n");
		fprintf(fileitp," Si  C2  C2  1 112.67  297.0  \n\n");
		fprintf(fileitp,"[ dihedraltypes ]\n");		
		fprintf(fileitp," Si  C2  C2  Si  3  0.1676  -1.7598   0.419   2.3464   0.0  0.0  \n\n");
		//[ moleculetype ]
		fprintf(fileitp,"[ moleculetype ]\n");
		fprintf(fileitp,"MEM   6 \n\n");
		//[ atoms ]
		fprintf(fileitp,"[ atoms ]\n");
		tempi=1;
		for(i=0;i<N;i++){
			fprintf(fileitp," %d %s %d  %s  %s  %d  %f  %f \n ",i+1,naz[AType[i]],1,"MEM",anaz[AType[i]],tempi,0.0, mass[i]);
			tempi++;
			if(tempi>20){
				tempi=1;
			}
		}
		fprintf(fileitp,"[ bonds ]\n");
		for(i=0;i<N;i++){
		//printf(" i %d bond [i] %d \n", i, bond[i]);
			for(j=0;j<bond[i];j++){
				lin1=i*MAXNEAR+j+MAXATOM;
				if(bond[lin1]>i){
					//printf(" lin1 %d bond[lin1] %d\n", lin1,bond[lin1]);
					fprintf(fileitp,"%d     %d    %d \n",i+1,bond[lin1]+1,1);
				}
			}
		}
		fprintf(fileitp,"[ angles ]\n");
		for(i=0;i<N;i++){
			printf("i %d angle[i] %d\n",i,angle1[i]);
			for(j=0;j<angle1[i];j++){
				lin1=i*MAXNEAR+j+MAXATOM;
				//if(angle2[lin1]>angle1[lin1]){
					fprintf(fileitp," %d   %d   %d   %d  \n",angle1[lin1]+1,i+1,angle2[lin1]+1,1);
				//}
			}
		}
		fprintf(fileitp,"[ dihedrals ]\n");
		for(i=0;i<N;i++){
			//printf("i %d tors[i] %d\n", i, tors1[i]);
			if(AType[i]==1){
				for(j=0;j<tors1[i];j++){
					lin1=i*MAXNEAR+j+MAXATOM;
					if(AType[tors1[lin1]]==3){
						if(AType[tors2[lin1]]==3){
							if(AType[tors3[lin1]]==1){
								fprintf(fileitp,"  %d   %d   %d   %d   %d \n",i+1,tors1[lin1]+1,tors2[lin1]+1,tors3[lin1]+1,3);
							}
						}
					}
				}
			}
		}
		fprintf(fileitp,"[ exclusions ]\n");
		fprintf(fileitp,"[ position_restraints ]\n");
		fprintf(fileitp,"[ constraints ]\n\n");
		fprintf(fileitp,"[ system ]\n");
		fprintf(fileitp,"generated \n");
		fprintf(fileitp,"[ molecules ] \n");
		fprintf(fileitp,"  %s  %d \n", "MEM", 1);
		
//			[ dihedrals ]
//			[ exclusions ]
//			[ position_restraints ]
//			[ constraints ]
//			[ system ]
//			generateg with membcreat v5
//			[ molecules ]
//       B     328
		
		
	fclose(fileitp);
}

void Replace(int N, float proc){
	int torep;
	int done;
	int random;
	int lin1;
	int cura;
	int i;
	int j;
	int k;
	int first;
	int second;
	torep=proc/100.0f*(TypeKol[2]+1);
	//printf("Replacing %d O atoms\n", torep);
	done=0;
	cura=Natom;
	while(done<torep){
		random=floor(ran2(&rseed)*Natom);
		if(AType[random]==2){
			//заменяем атом текущий
			AType[random]=3;
			second=bond[random*MAXNEAR+1+MAXATOM];	//запоминаем номер второго атома
			bond[random*MAXNEAR+1+MAXATOM]=cura;	//заменяем второй атом
			ax[cura]=ax[random]+0.1f;	//new coords
			ay[cura]=ay[random];
			az[cura]=az[random];
			bond[cura]=2;
			bond[cura*MAXNEAR+0+MAXATOM]=random;
			bond[cura*MAXNEAR+1+MAXATOM]=second;
			AType[cura]=3;
			//заменяем связь на втором атоме
			for(i=0;i<bond[second];i++){
				lin1=second*MAXNEAR+i+MAXATOM;
				if(bond[lin1]==random){
					bond[lin1]=cura;
				}
			}
			cura++;
			done++;
		}
	}
	Natom=cura;
}
