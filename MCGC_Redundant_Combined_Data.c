/*****************************Redundant Mutually Connected Giant Component ******************************************
This code can be redistributed and/or modified under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed ny the authors in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If you use this code please cite the following paper:

[1] F. Radicchi and G. Bianconi, 
 "Redundant interdependencies boost the robustness of multiplex networks" 
 Phys. Rev. X 7, 011013 (2017)

(c) G. Bianconi (email: ginestra.bianconi@gmail.com )
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


 This program reads data (edgelist) for  a multiplex network with M=3 layers 
 in the format :
 layer node1 node2 weight (weight to be disregarded)
 
 It consider a single random configuration of initial damage  
 and evaluates the Redudant Mutually Connected Giant Component as a function of thefraction of  nodes damaged
 The Message Passing algorithm for each single realization of the random damage
 is calculated on the same network.
 
 INPUT
 N total number of nodes filename of the edgelist

 
 OUTPUT
 
 The output file "Redundant_MP.txt"  contains two columns: 
 p=1-f the fraction of non damaged nodes, and the size of the RMCGC predicted by the Message Passing Algorithm
 
 The output file "Redundant_sim.txt"  contains two columns: 
 p=1-f the fraction of non damaged nodes  and the size of the RMCGC obtained by direct numerical calculations
 
*************************************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>



int *vis1,*vis2,*vis3,*size_cluster1,*size_cluster2,*size_cluster3,**knn1,**knn2,**knn3,*k1,*k2,*k3,c1,c2,c3,*occ,*dam1,*dam2,*dam3,N;


/*****************************************************************************
 Recurrence is a subrutine for calculating the giant component in each layer
 *****************************************************************************/
int Recurrence(int net, int i , int cluster_size, int ncluster){
	int j, n3,aus_cluster_size;
	
	cluster_size++;	
	if(net==1){
		vis1[i]=ncluster;
		for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
			if((vis1[j]==0)&&(dam1[j]==1)){
				aus_cluster_size=Recurrence(1,j, cluster_size, ncluster);
				cluster_size=aus_cluster_size;
			}
	
		}
	}
	if(net==2){
		vis2[i]=ncluster;
		for(n3=0;n3<k2[i];n3++){
			j=knn2[i][n3];			
			if((vis2[j]==0)&&(dam2[j]==1)){
			aus_cluster_size=Recurrence(2,j, cluster_size, ncluster);
			cluster_size=aus_cluster_size;
			}
		}
	}


	return cluster_size;
}
/*******************************************************************************************
 RecurrenceM is a subrutine for calulating the Redundant Mutually Connected Giant Component 
 *******************************************************************************************/

int RecurrenceM(int net, int i , int cluster_size, int ncluster){
	int j, n3,aus_cluster_size;
	float *xd1,*xd2,*xd3;
	
	cluster_size++;
	
	if(net==1){
		occ[i]=ncluster;
		for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
			if(dam1[j]==1){
			if(((vis2[j]==c2)&&(occ[j]==0)&&(dam2[j]==1))){
				aus_cluster_size=RecurrenceM(1,j, cluster_size, ncluster);
				cluster_size=aus_cluster_size;
			}
			}
		}
	}
	if(net==2){
		occ[i]=ncluster;
		for(n3=0;n3<k2[i];n3++){
			j=knn2[i][n3];		
			if(dam2[j]==1){
			if(((vis1[j]==c1)&&(occ[j]==0)&&(dam1[j]==1))){
				aus_cluster_size=RecurrenceM(2,j, cluster_size, ncluster);
				cluster_size=aus_cluster_size;
			}
		}
		}
	}
	

	return cluster_size;
}



int main(int argc, char** argv){
	int i,j,n,it,ncluster1,ncluster2,ncluster3, GMCC, cluster_size,m1,m2,m3,m1_aus,m2_aus,m3_aus,c1_aus,c2_aus,c3_aus,*sigma,nrun,nc;
	int s1,s2,Nc1,Nc2,Nc3,**n1,**n2,**n3,aus,aus1,aus2,aus3,**adj1,**adj2,**adj3,MCGC,nsum,nsumold,layer;
	float p,x,f,*xd1,*xd2,*xd3,*Sm;

	
	   char filec[60],**s,string1[50],string2[50];
	
	FILE *gp2,*gp,*ffp;
  

  /**************************************************************************
   open file for output 
  **************************************************************************/

		srand48(time(NULL));
	
	printf("Input: 1) N total number of nodes and 2) filename of multiplex network edgelist\n");
	scanf("%d %s",&N,filec);
	N=N+1;
	
	
	ffp=fopen(filec,"r");
	gp=fopen("MCGC.txt","w");		
	gp2=fopen("MCGC_MP.txt","w");	
	
	Sm=(float*)calloc(200,sizeof(float));
	vis1=(int*)calloc(N,sizeof(int));
	vis2=(int*)calloc(N,sizeof(int));
	vis3=(int*)calloc(N,sizeof(int));
	occ=(int*)calloc(N,sizeof(int));
	k1=(int*)calloc(N,sizeof(int));
	k2=(int*)calloc(N,sizeof(int));
	k3=(int*)calloc(N,sizeof(int));
	xd1=(float*)calloc(N,sizeof(float));
	xd2=(float*)calloc(N,sizeof(float));
	xd3=(float*)calloc(N,sizeof(float));
	dam1=(int*)calloc(N,sizeof(int));
	dam2=(int*)calloc(N,sizeof(int));
	dam3=(int*)calloc(N,sizeof(int));
	sigma=(int*)calloc(N,sizeof(int));
	knn1=(int**)calloc(N,sizeof(int*));
	knn2=(int**)calloc(N,sizeof(int*));
	knn3=(int**)calloc(N,sizeof(int*));
	n1=(int**)calloc(N,sizeof(int*));
	n2=(int**)calloc(N,sizeof(int*));
	n3=(int**)calloc(N,sizeof(int*));
	adj1=(int**)calloc(N,sizeof(int*));
	adj2=(int**)calloc(N,sizeof(int*));
	adj3=(int**)calloc(N,sizeof(int*));
		for(i=0;i<N;i++){
			knn1[i]=(int*)calloc(N,sizeof(int));
			knn2[i]=(int*)calloc(N,sizeof(int));
			knn3[i]=(int*)calloc(N,sizeof(int));
			adj1[i]=(int*)calloc(N,sizeof(int));
			adj2[i]=(int*)calloc(N,sizeof(int));
			adj3[i]=(int*)calloc(N,sizeof(int));
			n1[i]=(int*)calloc(N,sizeof(int));
			n2[i]=(int*)calloc(N,sizeof(int));
			n3[i]=(int*)calloc(N,sizeof(int));

		}
	size_cluster1=(int*)calloc(N,sizeof(int));
	size_cluster2=(int*)calloc(N,sizeof(int));
	size_cluster3=(int*)calloc(N,sizeof(int));

	for(i=0;i<N;i++){
		xd1[i]=drand48();
		xd2[i]=drand48();
		xd3[i]=drand48();
	}
	for(i=0;i<N;i++){
		k1[i]=0;
		k2[i]=0;
		k3[i]=0;
	}
	nsum=0;
	
	/***************************************************************************
	 Read data from edgelist of the  multiplex network and initialize the messages
	 **************************************************************************/
	while(fscanf(ffp,"%d %d %d %d", &layer,&i, &j,&n)!=EOF){
					if(layer==1){
				
				k1[i]++;
				k1[j]++;
				knn1[i][k1[i]-1]=j;
				knn1[j][k1[j]-1]=i;
				n1[i][j]=(int)(2.*drand48());
				n1[j][i]=(int)(2.*drand48());
				adj1[i][j]=1;
				adj1[j][i]=1;
				nsum+=(n1[i][j]+n1[j][i]);
			}
			if(layer==2){
				k2[i]++;
				k2[j]++;
				knn2[i][k2[i]-1]=j;
				knn2[j][k2[j]-1]=i;
				n2[i][j]=(int)(2.*drand48());
				n2[j][i]=(int)(2.*drand48());
				adj2[i][j]=1;
				adj2[j][i]=1;
				nsum+=(n2[i][j]+n2[j][i]);
			}
		}
					   		
		nc=0;
	
	for(f=0.;f<1.;f+=0.01){
		nc++;
		for(i=0;i<N;i++){
			if(xd1[i]<f){
				dam1[i]=0;}
			else {dam1[i]=1;}
			if(xd2[i]<f){
				dam2[i]=0;}
			else {dam2[i]=1;}
			if(xd3[i]<f){
				dam3[i]=0;}
			else {dam3[i]=1;}
			
		}
		/***********************************************************************************
		 ***Message passing algorithm
		 ***********************************************************************************/

		
		nsumold=0;
		while(fabs(nsum-nsumold)>0){
			nsumold=nsum;
			nsum=0;
			
			for (i=1;i<N;i++){
				for(j=1;j<N;j++){
					if((adj1[i][j]+adj2[i][j]+adj3[i][j])>0){
						
						aus1=1;
						for(n=0;n<k1[i];n++){
							if(knn1[i][n]!=j){
								aus1=aus1*(1-n1[knn1[i][n]][i]);}
						}
						aus2=1;
						for(n=0;n<k2[i];n++){
							if(knn2[i][n]!=j){
								aus2=aus2*(1-n2[knn2[i][n]][i]);}
						}
						aus3=1;
						for(n=0;n<k3[i];n++){
							if(knn3[i][n]!=j){
								aus3=aus3*(1-n3[knn3[i][n]][i]);}
						}
						
						aus=0;
							aus+=((1-aus1)*dam1[i])+aus1*dam1[i]*dam1[j]*adj1[i][j];
							aus+=((1-aus2)*dam2[i])+aus2*dam2[i]*dam2[j]*adj2[i][j];
						     aus+=((1-aus3)*dam3[i])+aus3*dam3[i]*dam3[j]*adj3[i][j];
						
						
						if(aus<2){
							aus=0;}
						else{
							aus=1;
						}
						n1[i][j]=adj1[i][j]*(1-aus1)*dam1[i]*dam1[j]*aus;
						n2[i][j]=adj2[i][j]*(1-aus2)*dam2[i]*dam2[j]*aus;
						n3[i][j]=adj3[i][j]*(1-aus3)*dam3[i]*dam3[j]*aus;
						nsum+=n1[i][j]+n2[i][j]+n3[i][j];
						
					}
				}
			}
		}
		MCGC=0;
		for (i=1;i<N;i++){
			aus1=1;
			for(n=0;n<k1[i];n++){
				aus1=aus1*(1-n1[knn1[i][n]][i]);}
			
			aus2=1;
			for(n=0;n<k2[i];n++){
				aus2=aus2*(1-n2[knn2[i][n]][i]);}
			
			aus3=1;
			for(n=0;n<k3[i];n++){
				aus3=aus3*(1-n3[knn3[i][n]][i]);}
			
			
			aus=(1-aus1)*dam1[i]+((1-aus2)*dam2[i])+((1-aus3)*dam3[i]);
			
			if (aus>=2){
				MCGC+=(1-aus1)*dam1[i]+(1-aus2)*dam2[i]+(1-aus3)*dam3[i];
			}
			
		}
		
	
	
	
				
	

	
	fprintf(gp,"%lf %lf \n",1-f,(float)MCGC/(3*(float)(N)));
		printf("%lf %lf \n",1-f,(float)MCGC/(3*(float)N));
		
		
		
		
		
		/***************************************************************************************
		 Calculates the giant component of each layer
		 ****************************************************************************************/
		ncluster1=0;
		ncluster2=0;
		ncluster3=0;
		for(i=0;i<N;i++){
			vis1[i]=0;
			vis2[i]=0;
			vis3[i]=0;
		}
		m1=0;
		m2=0;
		m3=0;
		for(n=0;n<N;n++){
			if((vis1[n]==0)&&(dam1[n]==1)){
				cluster_size=0;
				ncluster1++;
				cluster_size=Recurrence(1,n, cluster_size, ncluster1);
				size_cluster1[ncluster1]=cluster_size;
				if(cluster_size>m1){m1=cluster_size;
					c1=ncluster1;}
				
			}
			
			if((vis2[n]==0)&&(dam2[n]==1)){
				cluster_size=0;
				ncluster2++;
				cluster_size=Recurrence(2,n, cluster_size, ncluster2);
				size_cluster2[ncluster2]=cluster_size;
				if (cluster_size>m2) {m2=cluster_size;c2=ncluster2;}
				
			}
			
			if((vis3[n]==0)&&(dam3[n]==1)){
				cluster_size=0;
				ncluster3++;
				cluster_size=Recurrence(3,n, cluster_size, ncluster3);
				size_cluster2[ncluster2]=cluster_size;
				if (cluster_size>m3) {m3=cluster_size;c3=ncluster3;}
				
			}
		}
		
		
		Nc2=ncluster2;
		Nc1=ncluster1;
		Nc3=ncluster3;
		/*************************************************************************************
		 Calculates back and forth the giant components in each layer until the RMCGC is found 
		 GMCC is the size of the RMCGC
		 *************************************************************************************/
		GMCC=0;
		m1_aus=0;
		m2_aus=0;
		m3_aus=0;
		while(m1_aus!=m1){
			
			m1=m1_aus;
			m2=m2_aus;
			m3=m3_aus;
			m1_aus=0;
			m2_aus=0;
			m3_aus=0;
			for(i=0;i<N;i++)occ[i]=0;
			for(i=0;i<N;i++){
				if(dam2[i]==1){
					if(((vis1[i]==c1)&&(vis2[i]==c2)&&(occ[i]==0)&&(dam1[i]==1))||((vis2[i]==c2)&&(vis3[i]==c3)&&(occ[i]==0)&&(dam3[i]==1))){
						cluster_size=0;
						ncluster2++;
						cluster_size=RecurrenceM(2,i, cluster_size, ncluster2);
						if (cluster_size>m2_aus) {m2_aus=cluster_size;c2_aus=ncluster2;}
					}				
				}
			}
			c2=c2_aus;
			
			for(i=0;i<N;i++){vis2[i]=occ[i];occ[i]=0;}
			for(i=0;i<N;i++){
				if(dam3[i]==1){
					if(((vis1[i]==c1)&&(vis3[i]==c3)&&(occ[i]==0)&&(dam1[i]==1))||((vis2[i]==c2)&&(vis3[i]==c3)&&(occ[i]==0)&&(dam2[i]==1))){
						cluster_size=0;
						ncluster3++;
						cluster_size=RecurrenceM(3,i, cluster_size, ncluster3);
						if (cluster_size>m3_aus) {m3_aus=cluster_size;c3_aus=ncluster3;}
					}				
				}
			}
			c3=c3_aus;
			
			for(i=0;i<N;i++){vis3[i]=occ[i];occ[i]=0;}
			for(i=0;i<N;i++){
				if(dam1[i]==1){
					if(((vis1[i]==c1)&(vis2[i]==c2)&(occ[i]==0)&&(dam2[i]==1))||((vis3[i]==c3)&(vis1[i]==c1)&(occ[i]==0)&&(dam3[i]==1))){
						cluster_size=0;
						ncluster1++;
						cluster_size=RecurrenceM(1,i, cluster_size, ncluster1);
						if (cluster_size>m1_aus) {m1_aus=cluster_size;c1_aus=ncluster1;}
					}				
				}
			}
			c1=c1_aus;
			for(i=0;i<N;i++){vis1[i]=occ[i];occ[i]=0;}
			
		}
		GMCC=m1+m2+m3;


	fprintf(gp2,"%lf %lf \n",1-f,(float)GMCC/(3*(float)N));
}
fclose(gp2);
	
	
	fclose(gp);
			
	return 0;	
}

