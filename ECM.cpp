/*
compile:
g++ -g -O2 ECM.cpp -o ECM2.0 -lgmpxx -lgmp

ATTENTION THAT the programme is just suitable for finding
factor whoes length is no more than 50-digit.
smoothness parameter B1 must be no more than 43000000,
B2 must be no more than 2.4e11.


RECOMMANDED PARAMETERS

digits  B1  		B2     			    d
15		2000		1.5e5(147396)		2*3*5*7*11;
20    	11000    	1.9e6(1873422) 		2*3*5*7*11;
25		50000		1.3e7(12746592)   	2*3*5*7*11*13;
30		250000		1.3e8(128992510)  	2*3*5*7*11*13;
35  	1000000		1.0e9(1045563762)	2*3*5*7*11*13*17;
40		3000000		1.0e9(1045563762)	2*3*5*7*11*13*17;
45		3000000		5.3e9(5706890290) 	2*3*5*7*11*13*17;
45		1100000		3.5e10(35133391030)	2*3*5*7*11*13*17;
50 		43000000	2.4e11(240491351116)2*3*5*7*11*13*17*19;

--edit by GUIWEN LUO--
    --20180806--
*/

#include <iostream>
#include <math.h>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;


//store the binary form and NAF form of an integer
int BIN[1000],LENGTHBIN;//LENGTHBIN BIN(113)==7;
int *PTBIN=BIN;

int NAF[1000],LENGTHNAF;//LENGTHNAF NAF(113)==8;
int *PTNAF=NAF;

//There are 216816 primes below 3000000.
//There are 726517 primes below 11000000.
//There are 2604535 primes below 43000000.
mpz_class PRIMEBASE[2605000];
int LENGTHPRIMEBASE;
mpz_class *PTPRIMEBASE=PRIMEBASE;

mpz_class SETS[1000000];
mpz_class SETS_Qx[1000000];
int LENGTHSETS;
mpz_class *PTSETS=SETS;

mpz_class SETT[1000000];
mpz_class SETT_Qx[1000000];
int LENGTHSETT;
mpz_class *PTSETT=SETT;

mpz_class SIGMA[1000000];
mpz_class *PTSIGMA=SIGMA;


void DEC2NAF(int* rpt,int& length,const mpz_class& op){
    // 指针 rpt不能用 int* & rpt形式，因为我们用完指针之后要归位,使得 *PTNAF=NAF,这样才能复用,
mpz_class E=op,rem;
length=0;//length 首先也要归零，不然会出错！
int ee;

    while(E>0)
    {
        if((E%2)==1)
        {
            rem=(E%4);
            ee=mpz_get_ui(rem.get_mpz_t());
            *(rpt+length)=2-ee;
            E-=*(rpt+length);
        }
        else
        {
            *(rpt+length)=0;
        }
        E/=2;
        length+=1;
    }
}

void DEC2BIN(int* rpt,int& length,const mpz_class& op){

mpz_class q=op,rem;
length=0;//length 首先也要归零，不然会出错！
int ee;
    while(q>0)
    {
        rem=q%2;
       *(rpt+length)=mpz_get_ui(rem.get_mpz_t());
       q=q/2;
       length+=1;
    }
}


class PointInAffine
{
public:
    mpz_class Xa,Ya;//仿射坐标
    PointInAffine (){};//必须要有，不然系统不会为无参数的对象调用 constructor.
    PointInAffine (const mpz_class& a,const mpz_class& b){Xa=a;Ya=b;};
    inline void Print(){cout<<"("<<Xa<<","<<endl<<" "<<Ya<<") "<<endl;};
};

class PointInMontgomery
{
public:
    mpz_class X,Z,XpZ,XmZ;//Montgomery坐标
    PointInMontgomery (){};//必须要有，不然系统不会为无参数的对象调用 constructor.
    PointInMontgomery (const mpz_class& a,const mpz_class& b){X=a;Z=b;XpZ=a+b;XmZ=a-b;};//未检查大小，初始化时一定要注意0=<a,b<q
    inline void Print(){cout<<"("<<X<<","<<endl;
                        cout<<" "<<Z<<","<<endl;
                        cout<<" "<<XpZ<<","<<endl;
                        cout<<" "<<XmZ<<")"<<endl;
                        };
};

inline void Double(PointInMontgomery& DP,const PointInMontgomery& P,const mpz_class& Aa2d4,const mpz_class& n){
mpz_class u,v;
u=P.XpZ*P.XpZ;
u%=n;
v=P.XmZ*P.XmZ;
v%=n;

DP.X=u*v;
DP.X%=n;

DP.Z=(Aa2d4*(u-v)+v)%n;
DP.Z=DP.Z*(u-v)%n;

DP.XpZ=DP.X+DP.Z;
DP.XmZ=DP.X-DP.Z;
}


inline void Double_WST(PointInAffine& DQ, const PointInAffine& Q, const mpz_class& a, const mpz_class& n){
mpz_class t1,g,inv,ss,lambda;
t1=2*Q.Ya;
mpz_gcdext(g.get_mpz_t(),inv.get_mpz_t(),ss.get_mpz_t(),t1.get_mpz_t(),n.get_mpz_t());
lambda=(3*Q.Xa*Q.Xa+a)%n;
lambda=(lambda*inv)%n;
DQ.Xa=(lambda*lambda-2*Q.Xa)%n;
DQ.Ya=(lambda*(Q.Xa-DQ.Xa)-Q.Ya)%n;
}


inline void Add(PointInMontgomery& PpQ,const PointInMontgomery& P,const PointInMontgomery& Q,const PointInMontgomery& PmQ,const mpz_class& Aa2d4,const mpz_class& n){
mpz_class u,v,w;
u=(P.XpZ*Q.XmZ)%n;
v=(P.XmZ*Q.XpZ)%n;
w=(u+v);
w=(w*w)%n;
PpQ.X=(PmQ.Z*w)%n;

w=u-v;
w=(w*w)%n;
PpQ.Z=(PmQ.X*w)%n;

PpQ.XpZ=PpQ.X+PpQ.Z;
PpQ.XmZ=PpQ.X-PpQ.Z;
}

inline void Add_WST(PointInAffine& PpQ,const PointInAffine& P,const PointInAffine& Q, const mpz_class& a, const mpz_class& n){
mpz_class t1,g,inv,ss,lambda;
t1=Q.Xa-P.Xa;
mpz_gcdext(g.get_mpz_t(),inv.get_mpz_t(),ss.get_mpz_t(),t1.get_mpz_t(),n.get_mpz_t());
lambda=((Q.Ya-P.Ya)*inv)%n;
PpQ.Xa=(lambda*lambda-P.Xa-Q.Xa)%n;
PpQ.Ya=(lambda*(P.Xa-PpQ.Xa)-P.Ya)%n;
}

void Multi(PointInMontgomery& tmp,const mpz_class& k,const PointInMontgomery& P,const mpz_class& Aa2d4,const mpz_class& n){
PointInMontgomery tmp1,tmp2,tmp3;
tmp1=P;
Double(tmp2,P,Aa2d4,n);
DEC2BIN(PTBIN,LENGTHBIN,k);

for(int i=(LENGTHBIN-2);i>=1;i--){
    if (BIN[i]==1){
    Add(tmp3,tmp2,tmp1,P,Aa2d4,n);
    tmp1=tmp3;
    Double(tmp3,tmp2,Aa2d4,n);
    tmp2=tmp3;
    }
    else{
    Add(tmp3,tmp2,tmp1,P,Aa2d4,n);
    tmp2=tmp3;
    Double(tmp3,tmp1,Aa2d4,n);
    tmp1=tmp3;
    }
}

if (BIN[0]==1){
Add(tmp3,tmp2,tmp1,P,Aa2d4,n);
tmp1=tmp3;
}
else{
Double(tmp3,tmp1,Aa2d4,n);
tmp1=tmp3;
}
tmp=tmp1;
}

void Multi_WST_NAF(PointInAffine& tmp,const mpz_class& k,PointInAffine& Q,const mpz_class& a,const mpz_class& n){
PointInAffine P,NEGA_Q,tmp1;
DEC2NAF(PTNAF,LENGTHNAF,k);
P=Q;
NEGA_Q.Xa=Q.Xa;
NEGA_Q.Ya=n-Q.Ya;
for(int i=(LENGTHNAF-2);i>=0;i--){
Double_WST(tmp1,P,a,n);
P=tmp1;
if(NAF[i]==1){
Add_WST(tmp1,P,Q,a,n);
P=tmp1;
}
else{
if(NAF[i]==-1){
Add_WST(tmp1,P,NEGA_Q,a,n);
P=tmp1;
}
}
tmp=P;
}

}
void PrimeBase(mpz_class* rpt,int& length,const mpz_class& B){
length=0;
mpz_class a,tmp,tmp1;
a=mpz_class("1",10);
mpz_nextprime(tmp.get_mpz_t(),a.get_mpz_t());
 *(rpt+length)=tmp;
while (tmp < B){
length+=1;
tmp1=tmp;
mpz_nextprime(tmp.get_mpz_t(),tmp1.get_mpz_t());
 *(rpt+length)=tmp;
}
}

void MAKE_S_T(mpz_class* rptS,mpz_class* rptT,int& lengthS,int& lengthT,mpz_class& d,const mpz_class& B1, const mpz_class& B2){
lengthS=0;
lengthT=0;
mpz_class tmp,tmp1,mini,maxi;


if(B2<=1900000){
    d=2*3*5*7*11;
}
else if(B2<=130000000){
    d=2*3*5*7*11*13;
}
else if(B2<=35000000000){
    d=2*3*5*7*11*13*17;
}
else{
    d=2*3*5*7*11*13*17*19;
}

mini=(2*B1-d)/(2*d)+1;
maxi=(2*B2+d)/(2*d);

for(mpz_class i=mini;i<=maxi;i++){
	SETS[lengthS]=i*d;
	lengthS+=1;
}

for(mpz_class j=1;j<(d/2);j=j+2){
mpz_gcd(tmp1.get_mpz_t(),j.get_mpz_t(),d.get_mpz_t());
if(tmp1==1){
    SETT[lengthT]=j;
    lengthT+=1;
}
}

}

bool ECM(mpz_class& gg, mpz_class& sigma, const mpz_class& n,const mpz_class& B1,const mpz_class& B2,const int& c0, const int& c, const int& RANDOM){
cout<<"Factorize integer n: "<<n<<"."<<endl;
cout<<"Pareparing for factorization with B1 "<<B1<<" and B2 "<<B2 <<"."<<endl;

mpz_class dS,tmprandom,randomB;
randomB=mpz_class("31415926535",10);
gmp_randstate_t State;
gmp_randinit_default(State);

PrimeBase(PTPRIMEBASE,LENGTHPRIMEBASE,B1);
MAKE_S_T(PTSETS,PTSETT,LENGTHSETS,LENGTHSETT,dS,B1,B2);

clock_t TOTAL_t;
TOTAL_t= clock();

if(RANDOM==1){
    for(int i=0;i<=c-1;i++){
        mpz_urandomm(tmprandom.get_mpz_t(),State,randomB.get_mpz_t());
        SIGMA[i]=tmprandom+7;
        }
}
else{
    if(c<=LENGTHPRIMEBASE){
        for(int i=0;i<=c-1;i++){
            SIGMA[i]=PRIMEBASE[i+3];
        }
    }
    else{
        for(int i=0;i<=LENGTHPRIMEBASE-1;i++){
            SIGMA[i]=PRIMEBASE[i+3];
        }
        for(int i=LENGTHPRIMEBASE;i<=c-1;i++){
        mpz_urandomm(tmprandom.get_mpz_t(),State,randomB.get_mpz_t());
        SIGMA[i]=tmprandom+7;
        }
    }
}

mpz_class u,v,tmp,g,r,s,Aa2d4;
PointInMontgomery Q,tmpQ;

for(int i=c0-1;i<=c-1;i++){
sigma=SIGMA[i];
cout<<"test no."<<(i+1)<<" curve with sigma:"<<sigma<<"."<<endl;

clock_t STAGE1_t;
STAGE1_t = clock();

u=(sigma*sigma-5)%n;
v=4*sigma;
Q.X=(u*u*u)%n;
Q.Z=(v*v*v)%n;
Q.XpZ=Q.X+Q.Z;
Q.XmZ=Q.X-Q.Z;
tmp=(16*Q.X*v)%n;

mpz_gcdext(g.get_mpz_t(),r.get_mpz_t(),s.get_mpz_t(),tmp.get_mpz_t(),n.get_mpz_t());
tmp=v-u;
tmp=tmp*tmp*tmp;

Aa2d4=(tmp*(3*u+v)*r)%n;


//STAGE1
int k;
double lga,lgB1;

lgB1=log(mpz_get_d(B1.get_mpz_t()));

for(int i=0; i<=LENGTHPRIMEBASE-1;i++){
    lga=log(mpz_get_d(PRIMEBASE[i].get_mpz_t()));
    k=floor(lgB1/lga);
    for(int j=1;j<=k;j++){
        Multi(tmpQ,PRIMEBASE[i],Q,Aa2d4,n);
        Q=tmpQ;
}

mpz_gcd(gg.get_mpz_t(),Q.X.get_mpz_t(),n.get_mpz_t());

if(gg>1){
cout<<"Find a prime factor in STAGE 1:"<<endl;
cout<<gg<<endl;
cout<<"sigma: "<<sigma<<endl;
cout<<"Total factorization time is: "<<(clock()-TOTAL_t)/(CLOCKS_PER_SEC)<<"s."<<endl;

return true;
}

}
cout<<"STAGE1 time is "<<(clock()-STAGE1_t)/(CLOCKS_PER_SEC)<<"s."<<endl;

//STAGE 2

clock_t STAGE2_t;
STAGE2_t = clock();

PointInAffine QA,Q1,Q2,Q3,tmpQQ;
mpz_class X,Z,Z_inv,x,A,B,a;
X=Q.X;
Z=Q.Z;
mpz_gcdext(g.get_mpz_t(),Z_inv.get_mpz_t(),s.get_mpz_t(),Z.get_mpz_t(),n.get_mpz_t());
x=(X*Z_inv)%n;
A=Aa2d4*4-2;
B=((x+A)*x*x+x)%n;

tmp=3*B*B;
mpz_gcdext(g.get_mpz_t(),r.get_mpz_t(),s.get_mpz_t(),tmp.get_mpz_t(),n.get_mpz_t());
a=((3-A*A)*r)%n;

QA.Xa=((3*x+A)*r*B)%n;
QA.Ya=(3*B*r)%n;


//100000000 Max difference bewtween two neighbor primes is 282
mpz_class ARRAYFORSMALLK_KQA[1501];
for(int i=0;i<=1500;i++){ARRAYFORSMALLK_KQA[i]=0;}

Multi_WST_NAF(Q1,dS,QA,a,n);
Multi_WST_NAF(Q2,(SETS[0]/dS),Q1,a,n);
SETS_Qx[0]=Q2.Xa;
for(int i=1;i<=LENGTHSETS-1;i++){
    Add_WST(Q3,Q2,Q1,a,n);
    Q2=Q3;
    SETS_Qx[i]=Q2.Xa;
}

SETT_Qx[0]=QA.Xa;
Q3=QA;

for(int i=1;i<=LENGTHSETT-1;i++){
    tmp=SETT[i]-SETT[i-1];
    int kk=mpz_get_ui(tmp.get_mpz_t());
    if((ARRAYFORSMALLK_KQA[3*kk])!=1){
        Multi_WST_NAF(tmpQQ,kk,QA,a,n);
        ARRAYFORSMALLK_KQA[3*kk-2]=tmpQQ.Xa;
        ARRAYFORSMALLK_KQA[3*kk-1]=tmpQQ.Ya;
        ARRAYFORSMALLK_KQA[3*kk]=1;
    }
    else{
    tmpQQ.Xa=ARRAYFORSMALLK_KQA[3*kk-2];
    tmpQQ.Ya=ARRAYFORSMALLK_KQA[3*kk-1];
    }
    Add_WST(Q2,tmpQQ,Q3,a,n);
    Q3=Q2;
    SETT_Qx[i]=Q3.Xa;
}

mpz_class res;
res=1;

for(int i=0;i<=LENGTHSETS-1;i++){
    for(int j=0;j<=LENGTHSETT-1;j++){
        tmp=SETS_Qx[i]-SETT_Qx[j];
        res=(res*tmp)%n;
    }
}

mpz_gcd(gg.get_mpz_t(),res.get_mpz_t(),n.get_mpz_t());
if(gg>1){
    cout<<"Find a prime factor in STAGE 2:"<<endl;
    cout<<gg<<endl;
    cout<<"sigma: "<<sigma<<endl;
    cout<<"Total factorization time is: "<<(clock()-TOTAL_t)/(CLOCKS_PER_SEC)<<"s."<<endl;
    return true;
}

cout<<"STAGE2 time is "<<(clock()-STAGE2_t)/(CLOCKS_PER_SEC)<<"s."<<endl;

}
cout<<"We've tested "<<(c-c0+1)<<" curves but find no factor."<<endl;
cout<<"Total factorization time is: "<<(clock()-TOTAL_t)/(CLOCKS_PER_SEC)<<"s."<<endl;
return false;
}


int main(){

mpz_class f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,N;

N=mpz_class("298777079680636209728926753957151534\
56092168433989472516365634644101537789\
60413111693109599171717007066220856768\
26928556518363105076218043402519861108\
88478565527792110944761604797925911529\
02652843841510368831007499226931739933\
55808366922676333229835998998497712492\
28784711747748003757598085124777826598\
00341062835557205582040235002076760485\
78837876927418180945214920194972851554\
78023391744685179370519119919170850660\
3506807978474027",10);

f1=mpz_class("439883",10);
f2=mpz_class("1234567891",10);
f3=mpz_class("1732792378957",10);
f4=mpz_class("832328713672817",10);
f5=mpz_class("952921534336737871",10);
f6=mpz_class("1309672216786572122317",10);
f7=mpz_class("1208925819614750508040051",10);
f8=mpz_class("170141183460469231731687303715884105727",10);
f9=mpz_class("274483146844291876982863596310017101279",10);
f10=mpz_class("274483146896724957848529135013585755443",10);//f10 is 39-digit.
f11=mpz_class("802874755424222799061978026497620694203777347803",10);//f11 is 48-digit.
f12=mpz_class("1341601698241405401448957657201625462121694828679",10);//f12 is 49-digit.
f13=mpz_class("5295922065625750333475561956908250198163190237920373195647",10);
f14=mpz_class("3457645796086684459525372196446947094107078575533908472103",10);

mpz_class B1,B2;
B1=mpz_class("50000",10);
B2=mpz_class("13000000",10);
N=N/(f1*f2*f3*f4*f5*f6);

bool b;
mpz_class g,sigma;
b=ECM(g,sigma,N,B1,B2,1,20000,0);

return 0;
}





