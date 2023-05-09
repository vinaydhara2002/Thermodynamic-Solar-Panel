#include<bits/stdc++.h>
using namespace std;

class Node {
public:
    float data;
    Node* next;
};

void Data(Node* head,float T){
  Node* headref = head;
  int i = 0;
  float A[7];
  A[0] = 560 + 18.2*pow(T,1) + 0.222*pow(T,2) + 1.09e-03*pow(T,3);//P
  A[1] = 1.40 + (3.975e-3)*pow(T,1) + (3.14e-5)*pow(T,2) + (9.14e-7)*pow(T,3);//Cpl
  A[2] = 93.4-0.333*T;//Kl
  A[3] = 1233 - 3.81*pow(T,1) - 4.08e-03*pow(T,2) + 1.38e-04*pow(T,3) + 5.37e-06*pow(T,4);//rhol
  A[4] = 18.9 + 0.644*pow(T,1) + 8.84e-03*pow(T,2) + 5.96e-05*pow(T,3) + 1.11e-07*pow(T,4);//rhov
  A[5] = 218 - 2.83*pow(T,1) + 2.20e-2*pow(T,2);//mul
  A[6] = 214 - 0.907*pow(T,1) - 8.57e-03*pow(T,2) - 1.38e-04*pow(T,3);//lambdal
  while(1){
    headref->data = A[i];
    if(i==6){break;}
    headref->next = new Node();
    headref=headref->next;
    i++;
  }
}

float Cpg(float T){
  float A[3] = {12.268, 11.701, 4.636};
  float B[3] = {-0.0699, 0.0216, 0.0618};
  float C[3] = {0.000394, 0.0000868, -0.0000309};
  float D[3] = {-0.000000837, -0.000000112, 0};
  float X[3] = {0.381,0.179,0.439};
  float Sum=0;
  for(int i=0;i<3;i++){
   float Cp = 4.184*((A[i]*T)+(B[i]*(pow(T,2)))+(C[i]*(pow(T,3)))+(D[i]*(pow(T,4))));
   Sum+=Cp*X[i];
  }
  Sum/=86.2;
  return Sum;
}

float wall_temp(float Ti,float Ts){
  float To,Tn,ri,k,ro,mcg,dmcg,C1,C3,ha,d,hc,P,sigma;
  ri=1e-2;//m
  ro=1.25e-2;//m
  k=247;//W/m.K
  ha=35;//W/m^2.K
  d=4e-2;//m
  C1=(2*3.142*k)/(ha*d*log(ro/ri)+2*3.142*k);
  //Let the fractional phase change of coolant be given by following equation (function of temperature varying along length)
  mcg = pow((exp(Ts/Ti+1)-1)/(exp(1)-1),4);
  /*The differentiation of mass changing phase with respect to temperature is given by following equation 
    (function of temperature varying along length)*/
  dmcg = 4*pow((exp(Ts/Ti+1)-1),3)*exp(Ts/Ti+1)/(Ti*pow((exp(1)-1),4));
  //Let the initial guess for wall temperature To be 0°C
  To=0;
  Node* head = new Node();
  //Import the data of coolant for temperature Ts
  Data(head,Ts);
  int i=0;
  float A[7];
  while(head!=NULL){
    Node* temp = head;
    A[i] = head->data;
    head = head->next;
    delete temp;
    i++;
  }
  sigma=11.8*1e-3;
  C3=0.00122*(pow(A[2]*1e-3,0.79)*pow(A[1]*1e3,0.45)*pow(A[3],0.49)/(pow(sigma,0.5)*pow(A[5]*1e-6,0.29)*pow(A[6]*1e3,0.24)*pow(A[4],0.24)));
  while(1){
    P = (560 + 18.2*pow(To,1) + 0.222*pow(To,2) + 1.09e-03*pow(To,3))*1e3;
    hc = C3*pow(To-Ts,0.24)*pow(P-A[0]*1e3,0.75);
    hc*=0.05;//factor fs of 0.05
    hc+=68599;//hfc'=68599 W/m^2.K
    //We get new wall temperature using following equation 
    Tn = (ri*log(ro/ri)*(Ts)*hc)/(k*(hc*ri*log(ro/ri)/k+1-C1));
    if(Tn==0){break;}
    if((((Tn/To)-1)*100>-5)&&(((Tn/To)-1)*100<5)){break;}
    To=Tn;
  }
  //Using simpson's rule for 2nd degree polynomial to solve for length
  //heat absorbed by coolant in liquid phase (CpL.dT) is given by following
  float sum = (-1)*A[1]*1e3*(36.05/(4*90))*(1-mcg)*(hc*ri*log(ro/ri)/k+1-C1)/((Ts+273)*hc*2*3.142*ri*(C1-1));
  //heat absorbed by coolant in gas phase (CpG.dT) is given by following
  sum+=(-1)*Cpg(Tn)*1e3*(36.05/(4*90))*mcg*(hc*ri*log(ro/ri)/k+1-C1)/((Ts+273)*hc*2*3.142*ri*(C1-1));
  //heat absorbed by coolant during phase transition (lambdaL.(dm/dT)) is given by following
  sum+=(-1)*A[6]*1e3*(36.05/(4*90))*dmcg*(hc*ri*log(ro/ri)/k+1-C1)/((Ts+273)*hc*2*3.142*ri*(C1-1));
  return sum;
}

void length(float Ts){
  float /*Ts,*/Ti;
  /*cout<<"Please enter the coolant inlet temperature in K"<<endl;
  cin>>Ts;
  Ts-=273;*/
  Ti=(-1)*Ts;
  //Using simpson's rule to solve with n=1000
  float p = Ti/1000;
  float C0=wall_temp(Ti,Ts);
  /*cout<<"Ts (K)"<<"   "<<"length (m)"<<endl;
  cout<<Ts+273<<"      "<<C0*p/3<<endl;*/
  int i=1;
  while(Ts<=(-1)*p){
    Ts+=p;
    if(i%2==1){
      C0 += 4*wall_temp(Ti,Ts);
      //cout<<Ts+273<<"   "<<C0*p/3<<endl;
    }
    else if(i%2==0){
      C0 += 2*wall_temp(Ti,Ts);
      //cout<<Ts+273<<"   "<<C0*p/3<<endl;
    }
    i++;
  }
  Ts=0;
  C0 += wall_temp(Ti,Ts);
  cout<</*Ts+273<<"      "<<*/C0*p/3<<endl;
  /*C0*=p/3;
  cout<<"length  =  "<<C0<<" meters";
  return 0;*/
}

int main(){
  int T=-50;
  while(T<=-10){
    length(T);
    T++;
  }
  return 0;
}
