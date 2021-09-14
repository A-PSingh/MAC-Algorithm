#include <stdio.h>
#include <math.h>
#define length 129
float error_in_u (double U[length+1][length],double V[length+1][length]);
float error_in_p (double a1[length+1][length+1],double a2[length+1][length+1]);
int main () 
{
	float Re=400,dx=1.0/(length-1),dy=1.0/(length-1),dt=0.001;
   double u[length+1][length],u_new[length+1][length],v[length][length+1],v_new[length][length+1],p[length+1][length+1],p_new[length+1][length+1];
   double error,error_p;
   double RHSU[length][length]={0},RHSV[length][length]={0};
   int i,j,iteration=0;
   double l,m,n,o;
   double uvel[length][length],vvel[length][length];
   for(j=0;j<length+1;j++){
            for(i=0;i<length;i++){
                       u[j][i]=0;
                       u_new[j][i]=0;
                  }}
                  
    for(i=0;i<length;i++){
            u[length][i]=1;
            u[length-1][i]=1;
        }
    for(j=0;j<length;j++){
            for(i=0;i<length+1;i++){
                  v[j][i]=0;
         }}
         
    for(j=0;j<length+1;j++){
            for(i=0;i<length+1;i++){
                   p[j][i]=0;
                   p_new[j][i]=0;
            }}
            
    for(j=0;j<length;j++){
            for(i=0;i<length;i++){
                    RHSU[j][i]=0;
                    RHSV[j][i]=0;
          }}
 
 //Do while LOOP 
 

do
{
    for(j=1;j<=length-1;j++){
            for(i=1;i<=length-2;i++){
                l=(dt/dx)*0.25*(pow(((u[j][i+1]+u[j][i])),2)-pow(((u[j][i]+u[j][i-1])),2));
                m=(dt/dy)*(((u[j][i]+u[j+1][i])*(v[j][i]+v[j][i+1])*0.25)-((u[j][i]+u[j-1][i])*(v[j-1][i]+v[j-1][i+1])*0.25));
                n=((dt/(dx*dx*Re))*(u[j][i-1]-(2*u[j][i])+u[j][i+1]));
                o=((dt/(dy*dy*Re))*(u[j-1][i]-(2*u[j][i])+u[j+1][i]));
                RHSU[j][i]=u[j][i]-l-m+n+o; 
           }}
           
    for(j=1;j<=length-2;j++){
            for(i=1;i<=length-1;i++){
                l=(dt/dx)*(((u[j][i]+u[j+1][i])*(v[j][i]+v[j][i+1])*0.25)-((u[j][i-1]+u[j+1][i-1])*(v[j][i]+v[j][i-1])*0.25));
                m=(dt/dy)*(pow(((v[j+1][i]+v[j][i])),2)-pow(((v[j][i]+v[j-1][i])),2))*0.25;
                n=((dt/(dx*dx*Re))*(v[j][i-1]-(2*v[j][i])+v[j][i+1]));
                o=((dt/(dy*dy*Re))*(v[j-1][i]-(2*v[j][i])+v[j+1][i]));
                RHSV[j][i]=v[j][i]-l-m+n+o;
            }}
    error=0;
do
{
   error_p=0;
    for(j=1;j<=length-1;j++){
            for(i=1;i<=length-1;i++){
                    p[j][i]=0.25*((p[j][i+1]+p[j][i-1]+p[j+1][i]+p[j-1][i])-((dx*dx/dt)
					              *(((RHSU[j][i]-RHSU[j][i-1])/dx)+((RHSV[j][i]-RHSV[j-1][i])/dy))));
         }} 
    error_p=error_in_p(p_new,p); 
    for(j=0;j<length+1;j++){
            for(i=0;i<length+1;i++){
                    p_new[j][i]=p[j][i];
            }}
    }while(error_p>0.001); 
    
    for(j=1;j<length;j++){
            for(i=1;i<length-1;i++){
                    u[j][i]=((-dt/dx)*(p[j][i+1]-p[j][i]))+RHSU[j][i];
            }}
            
    for(j=1;j<length-1;j++){
            for(i=1;i<length;i++){
                    v[j][i]=((-dt/dx)*(p[j+1][i]-p[j][i]))+RHSV[j][i];
            }}
            
    for(i=0;i<=length-1;i++){
    	//Top boundary condition
            u[length][i]=2-u[length-1][i]; 
        //Bottom boundary    
            u[0][i]=-u[1][i]; 
        }
    for(j=2;j<=length-2;j++){
        //Left boundary
	        u[j][0]=0;
		//Right boundary	 
            u[j][length-1]=0; 
        }
    for(i=0;i<=length;i++){
        //Top boundary
		    v[length-1][i]=0;
        //Bottom boundary
		   v[0][i]=0; 
        }
    for(j=2;j<=length-3;j++){ 
        //Left boundary
		   v[j][0]=-v[j][1]; 
        //Right boundary
		    v[j][length]=-v[j][length-1]; 
        }
    for(i=0;i<=length;i++){
            p[length][i]=p[length-1][i];
            p[0][i]=p[1][i];
        }
    for(j=0;j<=length-2;j++){
            p[j][0]=p[j][1];
            p[j][length]=p[j][length-1];
        }
    error=error_in_u(u_new,u);
    iteration=iteration+1;
    printf("error: %f\t\titerations: %d\n",error,iteration);
    for(j=0;j<length+1;j++){
            for(i=0;i<length;i++){
                u_new[j][i]=u[j][i];}
        }
    for(j=0;j<length;j++){
            for(i=0;i<length+1;i++){
                v_new[j][i]=v[j][i];
        }}
        
}while(error>1e-4); 
   double u_mean[length];
   for(i=0;i<length;i++){
            u_mean[i]=(u[i][64]+u[i+1][64])*0.5;
       }
    for(i=0;i<129;i++){
            printf("%lf\n",u_mean[i]);
        }
    for(i=0;i<length;i++){
           for(j=0;j<length;j++){
                uvel[i][j]=(u[i][j]+u[i+1][j])*0.5;
                vvel[i][j]=(v[i][j]+v[i][j+1])*0.5;
       }}
        
   double sf[length][length]={0};
   sf[0][0]=0;
    for(j=0;j<length-1;j++){
            sf[j+1][0]=sf[j][0]-(uvel[j+1][0]*dy);
        }
    for(j=0;j<length;j++){
            for(i=0;i<length-1;i++){
                   sf[j][i+1]=sf[j][i]+((vvel[j][i+1])*dx);
        }}
 
//PRINTING RESULTS
 
    FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6;
    fp1=fopen("u_velocity.dat","w");
    for(i=0;i<length;i++){
            for(j=0;j<length;j++){
                  fprintf(fp1,"%lf\t%lf\t%lf\t%lf\n",(i*dx),(j*dy),uvel[i][j]);
                  printf("%lf\t%lf\t%lf\t%lf\n",(i*dx),(j*dy),uvel[i][j]);
       }}
    fclose(fp1);
    fp2=fopen("v_velocity.dat","w");
    for(i=0;i<length;i++){
            for(j=0;j<length;j++){
                   fprintf(fp2,"%lf\t%lf\t%lf\t%lf\n",(i*dx),(j*dy),vvel[i][j]);
                   printf("%lf\t%lf\t%lf\t%lf\n",(i*dx),(j*dy),vvel[i][j]);
        }} 
    fclose(fp2);
    fp3=fopen("stream_function.dat","w");
    for(i=0;i<length;i++){
            for(j=0;j<length;j++){
                    fprintf(fp3,"%lf\t%lf\t%lf\n",(i*dx),(j*dy),sf[i][j]);
                    printf("%lf\t%lf\t%lf\n",(i*dx),(j*dy),sf[i][j]);
       }}
    fclose(fp3);
    fp4=fopen("pressure.dat","w"); 
    for(i=0;i<length;i++){
            for(j=0;j<length;j++){
                    fprintf(fp4,"%lf\t%lf\t%lf\n",(i*dx),(j*dy),p[i][j]);   
         }}
    fclose(fp4);
   // x=0.5 and y=0.5
   fp5=fopen("u_mid.dat","w");
    for(j=0;j<length;j++){ 
            fprintf(fp5,"%lf\t%lf\n",uvel[j][64],(j*dy));
            printf("%lf\t%lf\n",uvel[j][64],(j*dy));
       }
    fclose(fp5);
    fp6=fopen("v_mid.dat","w");
    for(j=0;j<length;j++){
            fprintf(fp6,"%lf\t%lf\n",(j*dy),vvel[64][j]);
            printf("%lf\t%lf\n",(j*dy),vvel[64][j]);
       }
    fclose(fp6);
   return 0;
}


float error_in_u (double U[length+1][length],double V[length+1][length])
{
	
float error=0;
int i,j;
for(i=1;i<length;i++){
for(j=1;j<length-1;j++){
error=error+fabs(U[i][j]-V[i][j]);
}
} 
return error;
}
float error_in_p (double a1[length+1][length+1],double a2[length+1][length+1])
 {
float error=0;
int i,j;
for(i=1;i<length;i++){
for(j=1;j<length;j++){
error=error+fabs(a1[i][j]-a2[i][j]); 
}
} 
return error;
}

