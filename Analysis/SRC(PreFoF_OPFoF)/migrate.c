#include<stdio.h>
#include<math.h>
#include<mpi.h>
#include<stdlib.h>
#include<sys/times.h>
#include "fof.h"
#define max(A,B) ( (A) > (B) ? (A) : (B) )
#define min(A,B) ( (A) < (B) ? (A) : (B) )
void fresize(int);
#define TICKS 100.
void migrate_(size_t *np,int myid,int nid,int local_nz,int local_z_start,
		int nz,particle *bp, float *work)
{
  	particle tmprv;
	int ntmp,zdbound,zubound,i,j,unp,zp;
	int nget,nput,src,dest,numrecv,ierror,numsend,nend,npold;
	MPI_Status status;
	int commsize,zup,zdown;
	float *array[Na];
	int *intwork;
	float *fltwork;

	/*
	 * These two lines are to reduce the call on Operating system for
	 * more memory.  So reduce or free current massive memory usage
	 * into free list.
	*/
	/*
  	work = (float *)dlrealloc(work,densize*sizeof(float));
	den = work;
	*/

	float sstep;


	intwork = (int *) work;
	fltwork = (float *) work;

	
	zdbound = local_z_start;
	zubound = local_z_start+local_nz;
	do {
		unp = *np;
		/*
		 *  These lines are for the case
		 *  **********-----, where - is for migration.
		 */
		i = unp;
		do {
			i = i - 1;
    		zp = (int) (bp+i)->z;
	    	zup = abs(zubound-zp);
	    	zup = min(zup,nz-zup);
	    	zdown = abs(zdbound-zp);
	    	zdown = min(zdown,nz-zdown);
		} while( (zp >= zubound || zp < zdbound) && zup <= zdown );
		unp = i + 1;
		/*
		 */
		i = 0;
        while(i < unp){
			zp = (int) (bp+i)->z;
 	      	zup = abs(zubound - zp);
        	zup = min(zup,nz-zup);
			zdown = abs(zdbound-zp);
			zdown = min(zdown,nz-zdown);
			if((zp >= zubound || zp < zdbound) && zup <= zdown){
				tmprv = *(bp+i);
				do {
					unp = unp - 1;
					zp = (int) (bp+unp)->z;
					zup = abs(zubound-zp);
					zup = min(zup,nz-zup);
					zdown = abs(zdbound-zp);
					zdown = min(zdown,nz-zdown);
				} while((zp>=zubound || zp<zdbound) && zup<=zdown && unp>i);
				/* exchange two elements */
				(bp+i) = *(bp+unp);
				/*
				*(indx+i) = *(indx+unp);
				*/
				*(bp+unp) = tmprv;
				/*
				*(indx+unp) = ntmp;
				*/
			}
			i++;
		}
		src = (myid-1+nid) % nid;
		dest = (myid+1+nid) % nid;
  		numsend = *np-unp;
		numrecv = numsend;
		MPI_Sendrecv_replace(&numrecv,1,MPI_INT,dest,1,
				src,1,MPI_COMM_WORLD,&status);
		/*
		 * resize memory(x,y,z,vx,vy,vz,indx);
		 */
		/*
		 * total # of particles to send;
		 */
   		MPI_Reduce(&numsend,&commsize,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
   		MPI_Bcast(&commsize,1,MPI_INT,0,MPI_COMM_WORLD);
		if(commsize != 0){
		    nend = max(numrecv+unp,*np);
      		fresize(nend);
	     	if(myid == 0) printf("total upper sending .... %d \n",commsize);
			for(i=0;i<Na;i++){
				MPI_Sendrecv(array[i]+unp,numsend,MPI_FLOAT,dest,i,work,
						numrecv,MPI_FLOAT,src,i,MPI_COMM_WORLD,&status);
				for(j=0;j<numrecv;j++) *(array[i]+unp+j)=*(work+j);
			}
		    *np = *np-numsend+numrecv;
      		fresize(*np);
		}
	} while(commsize != 0);

	do {
		unp = *np;
		/*
		 *  These lines are for the case
		 *  **********-----, where - is for migration.
		 */
		i = unp;
		do {
			i = i - 1;
    		zp = (int) (bp+i)->z;
		} while(zp >= zubound || zp < zdbound);
		unp = i + 1;
		/*
		 */
		i = 0;
        while( i < unp){
			zp = (int) (bp+i)->z;
			if((zp >= zubound || zp < zdbound)){
				tmprv = *(bp+i);
				/*
				ntmp = *(indx+i);
				*/
				do {
					unp = unp - 1;
					zp = (int) (bp+unp)->z;
				} while((zp >= zubound || zp < zdbound)&& unp > i) ;
				/* change two elements */
				*(bp+i) = *(bp+unp);
				*(bp+unp) = tmprv;
			}
			i++;
		}
		src = (myid+1+nid)%nid;
		dest = (myid-1+nid)%nid;
  		numsend = *np-unp;
		numrecv = numsend;
		MPI_Sendrecv_replace(&numrecv,1,MPI_INT,dest,1,
		     src,1,MPI_COMM_WORLD,&status);
		/*
		 * resize memory(x,y,z,vx,vy,vz,indx);
		 */
		/*
		 * total # of particles to send;
		 */
   		MPI_Reduce(&numsend,&commsize,1,MPI_INT,MPI_SUM, 0,MPI_COMM_WORLD);
   		MPI_Bcast(&commsize,1,MPI_INT,0,MPI_COMM_WORLD);
	    fltwork = (float *) work;
		if(commsize != 0){
		    nend = max(numrecv+unp,*np);
      		fresize(nend);
	     	if(myid == 0) printf("total down sending .... %d \n",commsize);
			for(i=0;i<Na;i++){
				MPI_Sendrecv(array[i]+unp,numsend,MPI_FLOAT,dest,i,work,
						numrecv,MPI_FLOAT,src,i,MPI_COMM_WORLD,&status);
				for(j=0;j<numrecv;j++) *(array[i]+unp+j)=*(work+j);
			}
			/*
			MPI_Sendrecv((indx+unp),numsend,MPI_INT,dest,0,intwork,
					numrecv,MPI_INT,src,0,MPI_COMM_WORLD,&status);
			for(j=0;j<numrecv;j++) *(indx+unp+j)=*(intwork+j);
			*/
			/*
			 * resize memory
			 */
		    *np = *np-numsend+numrecv;
      		fresize(*np);
		}
	} while(commsize != 0);
	/*
	printf("%d %d %f %f \n",myid,*np,*(z+*np-2),*(z+*np-1));
	fflush(stdout);
	MPI_Finalize();
	exit(0);
	*/
}
#include<unistd.h>
void fresize(int np)
{
	return;
	/*
    x = (float *)realloc(x,np*sizeof(float));
    y = (float *)realloc(y,np*sizeof(float));
    z = (float *)realloc(z,np*sizeof(float));
    vx = (float *)realloc(vx,np*sizeof(float));
    vy = (float *)realloc(vy,np*sizeof(float));
    vz = (float *)realloc(vz,np*sizeof(float));
    indx = (int *)realloc(indx,np*sizeof(int));
	*/
}
