#include<fstream>
#include<iostream>
#include<sstream>
#include<math.h>
#include<cstdlib>
#include<iomanip>
//#include<omp.h>
#include"ran.h"
#include"gasdev.h"
#include <vector>
#include <algorithm>
#include <time.h>
//#include <conio.h>

inline double function(double x1, double a1, double s1){
	return(x1*x1*x1-x1+s1*a1);
}

inline double sq(double x1){
	return(x1*x1);
}

int main(int argc,char **argv){
    //number of tubes and area of simulation
	int numx=256;
	int numy=256;
    int numt=10000000;//300image
    int datM=10000;
    int dat=0;

    double dx=0.1;
    double dt=0.0001;
    double td, usum, asum;
    int jmm, jm, jp, jpp,
    kmm, km, kp, kpp,
    i,j,k;
    int ov0;

    double *u,*unew,*a,*anew;
    u=new double[numx*numy];
    unew=new double[numx*numy];
    a=new double[numx*numy];
    anew=new double[numx*numy];

    //related to timer for the log
    clock_t start,end;
    double stime;
    time_t timer;
    int opi=0;

    //here we first put a=0.5
    double q,D, alpha,beta,conc,xc,yc,size,S,M,ep;//a,q;
    q=0.01;
    D=0.01;//tiisakusuru_kakusannkou0.1->0.01_h30626
    alpha=0.1;//0.3
    beta=0.3;
    conc=-0.6;//0.4_zenntaino_noudo_h30626
    size=3.;//0.5_syoki_size
    ep=0.1;
    S=0.5;//1.0
    M=1.;
    xc=(double)numx*dx*0.5;
    yc=(double)numy*dx*0.5;
    
    //here make a output folder
    std::ostringstream command,name;
    std::ofstream os,os2;
    name.str("");
    name<<q<<"_"<<alpha<<"_"<<beta<<"_"<<conc<<"_"<<size<<"_"<<S<<"_"<<M<<"_"<<ep<<"_"<<numx<<"_"<<numy;
    
    command.str("");
    command<<"mkdir dat";
    system((command.str()).c_str());

    command.str("");
    command<<"mkdir rawu_"<<name.str();
    system((command.str()).c_str());

    command.str("");
    command<<"mkdir rawa_"<<name.str();
    system((command.str()).c_str());
    
    command.str("");
    command<<"./dat/dat_"<<name.str()<<".txt";
    os2.open((command.str()).c_str());

    //initial condition---startup
    //later I should modify this to be packed
    //at this moment this is pure random distribution
    for(j=0; j<numx; j++){
    	for(k=0; k<numy;k++){
//    		double ang=atan2((double)k*dx-yc,(double)j*dx-xc);
//   		if(((sq((double)j*dx-xc)+sq((double)k*dx-yc))-sq(size+0.1*cos(5.*ang)+0.2*cos(7.*ang)))<0.)
//    		for(int p=0;p<3;p++){
//    			for(int q=0;q<3;q++){
    				//if((sq((double)j*dx-xc)+sq((double)k*dx-yc)-(sq((double)j*dx-xc-p+1)+sq((double)k*dx-yc-q+1)-1)<0.)){
    				if(((sq((double)j*dx-xc)+sq((double)k*dx-yc)-5)<0.)){
    					u[j+k*numx]=1.;
    				}
    				else {
    					u[j+k*numx]=conc+0.005*(2.*ran()-1.);//0.05->0.005
    					a[j+k*numx]=0.;
//    			}
//   		}
    	}
    }
}
command.str("");
command<<"./rawu_"<<name.str()<<"/"<<std::setw(4)<<std::setfill('0')<<opi<<".raw";
os.open((command.str()).c_str(),std::ios::binary);
os.write((char*)u,sizeof(double)*numx*numy);
os.close();
command.str("");
command<<"./rawa_"<<name.str()<<"/"<<std::setw(4)<<std::setfill('0')<<opi<<".raw";
os.open((command.str()).c_str(),std::ios::binary);
os.write((char*)a,sizeof(double)*numx*numy);
os.close();

usum=0.;
asum=0.;
ov0=0;
for(j=0; j<numx;j++){
	for(k=0; k<numy; k++){
		td=u[j+k*numx];
		usum+=td;
		if(td>0.){
			ov0++;
		}
		asum+=a[j+k*numx];
	}
}

std::cout<<"img:"<<opi<<"\tsumu:"<<std::setprecision(8)<<usum
<<"\tsuma:"<<std::setprecision(8)<<asum
<<"\tu>0:"<<ov0<<"\n";
os2<<opi<<"\t"<<usum<<"\t"<<asum<<"\t"<<ov0<<"\n";
opi++;

    //start of time step
start=clock();
for(i=0; i<numt; i++){
	for(j=0;j<numx;j++){
		for(k=0; k<numy;k++){
                //here peiodic boundary condition is imposed.
			if(0==j){
				jmm=numx-2;
				jm=numx-1;
			}
			else if(1==j){
				jmm=numx-1;
				jm=j-1;
			}
			else{
				jmm=j-2;
				jm=j-1;
			}

			if(j==(numx-1)){
				jpp=1;
				jp=0;
			}
			else if(j==(numx-2)){
				jpp=0;
				jp=j+1;
			}
			else{
				jpp=j+2;
				jp=j+1;
			}

			if(0==k){
				kmm=numy-2;
				km=numy-1;
			}
			else if(1==k){
				kmm=numy-1;
				km=k-1;
			}
			else{
				kmm=k-2;
				km=k-1;
			}

			if(k==(numy-1)){
				kpp=1;
				kp=0;
			}
			else if(k==(numy-2)){
				kpp=0;
				kp=k+1;
			}
			else{
				kpp=k+2;
				kp=k+1;
			}


                //here acutual simple method is going on
			unew[j+k*numx]=u[j+k*numx]
			+M*dt*(function(u[jp+k*numx],a[jp+k*numx],S)
				+function(u[jm+k*numx],a[jm+k*numx],S)
				+function(u[j+kp*numx],a[j+kp*numx],S)
				+function(u[j+km*numx],a[j+km*numx],S)
				-4.*function(u[j+k*numx],a[j+k*numx],S))/(dx*dx)
			-M*q*dt*(u[jpp+k*numx]+u[jmm+k*numx]+u[j+kpp*numx]+u[j+kmm*numx]
				+2.*(u[jp+kp*numx]+u[jp+km*numx]+u[jm+kp*numx]+u[jm+km*numx])
				-8.*(u[j+kp*numx]+u[j+km*numx]+u[jp+k*numx]+u[jm+k*numx])
				+20.*u[j+k*numx])/(dx*dx*dx*dx);
			anew[j+k*numx]=a[j+k*numx]
			+dt*D*(a[jp+k*numx]+a[jm+k*numx]+a[j+kp*numx]+a[j+km*numx]-4.*a[j+k*numx])/(dx*dx)
			-dt*alpha*0.5*(1.-u[j+k*numx])*a[j+k*numx]+dt*beta*0.5*(1.+u[j+k*numx]);
		}

	}
        //n+1->nの入れ替え
	for(j=0; j<numx;j++){
		for(k=0; k<numy; k++){
			u[j+k*numx]=unew[j+k*numx];
			a[j+k*numx]=anew[j+k*numx];
		}
	}
	dat++;
	if(dat%datM==0){

		usum=0.;
		asum=0.;
		ov0=0;
		for(int j=0; j<numx;j++){
			for(int k=0; k<numy; k++){
				td=unew[j+k*numx];
				usum+=td;
				if(td>0.){
					ov0++;
				}
				asum+=a[j+k*numx];
			}
		}
            //std::cout<<"img:"<<opi<<"\tsumu:"<<std::setprecision(8)<<usum
            //    <<"\tsuma:"<<std::setprecision(8)<<asum
            //    <<"\tu>0:"<<ov0<<"\n";
		os2<<opi<<"\t"<<usum<<"\t"<<asum<<"\t"<<ov0<<"\n";
		command.str("");
		command<<"./rawu_"<<name.str()<<"/"<<std::setw(4)<<std::setfill('0')<<opi<<".raw";
		os.open((command.str()).c_str(),std::ios::binary);
		os.write((char*)u,sizeof(double)*numx*numy);
		os.close();
		command.str("");
		command<<"./rawa_"<<name.str()<<"/"<<std::setw(4)<<std::setfill('0')<<opi<<".raw";
		os.open((command.str()).c_str(),std::ios::binary);
		os.write((char*)a,sizeof(double)*numx*numy);
		os.close();
		opi++;
	}
}
end=clock();
stime=(double)(end-start)/(double)CLOCKS_PER_SEC;
time(&timer);
command.str("");
command<<"log.txt";
os.open((command.str()).c_str(),std::ios::out | std::ios::app);
os<<name.str()<<"\t"<<"time:\t"<<stime<<"\tLocaltime:\t"<<ctime(&timer)<<"\n";
os.close();
os2.close();
delete []u;
delete []unew;
delete []a;
delete []anew;
return 0;
}













