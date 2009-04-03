#include<malloc.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

#define fillincrement 0.1

float taugam,**topo,block,*heightvec,**mask,**height,delta;
int *heightvecind,*iup,*jup,*idown,*jdown,lattice_size_x,lattice_size_y,change;

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(float *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
      int  i,**m;

       /*allocate pointers to rows */
        m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
      m -= nrl;

       /*allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
      m[i] -= ncl;
      }
       /* return pointer to array of pointers to rows */
        return m;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)
);
            m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000

void indexx(n,arr,indx)
float arr[];
int indx[],n;
{
        unsigned long i,indxt,ir=n,itemp,j,k,l=1;
        int jstack=0,*istack;
        float a;

        istack=ivector(1,NSTACK);
        for (j=1;j<=n;j++) indx[j]=j;
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                indxt=indx[j];
                                a=arr[indxt];
                                for (i=j-1;i>=1;i--) {
                                        if (arr[indx[i]] <= a) break;
                                        indx[i+1]=indx[i];
                                }
                                indx[i+1]=indxt;
                        }
                        if (jstack == 0) break;
                        ir=istack[jstack--];
                        l=istack[jstack--];
                } else {
                        k=(l+ir) >> 1;
                        SWAP(indx[k],indx[l+1]);
                        if (arr[indx[l+1]] > arr[indx[ir]]) {
                                SWAP(indx[l+1],indx[ir])
                        }
                        if (arr[indx[l]] > arr[indx[ir]]) {
                                SWAP(indx[l],indx[ir])
                        }
                        if (arr[indx[l+1]] > arr[indx[l]]) {
                                SWAP(indx[l+1],indx[l])
                        }
                        i=l+1;
                        j=ir;
                        indxt=indx[l];
                        a=arr[indxt];
                        for (;;) {
                                do i++; while (arr[indx[i]] < a);
                                do j--; while (arr[indx[j]] > a);
                                if (j < i) break;
                                SWAP(indx[i],indx[j])
                        }
                        indx[l]=indx[j];
                        indx[j]=indxt;
                        jstack += 2;
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                        } else {
                                istack[jstack]=j-1;
                                istack[jstack-1]=l;
                                l=i;
                        }
                }
        }
        free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP

void setupgridneighbors()
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      {idown[i]=i-1;
       iup[i]=i+1;}
     idown[1]=1;
     iup[lattice_size_x]=lattice_size_x;
     for (j=1;j<=lattice_size_y;j++)
      {jdown[j]=j-1;
       jup[j]=j+1;}
     jdown[1]=1;
     jup[lattice_size_y]=lattice_size_y;
}

void addiflessthanthreshold(i,j)
int i,j;
{    float slope,slopex,slopey; 

     height[i][j]+=block;
     slopex=height[i][j]-height[iup[i]][j];
     if (height[i][j]-height[idown[i]][j]>slopex) slopex=height[i][j]-height[idown[i]][j];
     slopey=height[i][j]-height[i][jup[j]];
     if (height[i][j]-height[i][jdown[j]]>slopey) slopey=height[i][j]-height[i][jdown[j]];
     slope=sqrt(slopex*slopex+slopey*slopey)/delta; 
     /* this version implements slope-dependent threshold 
           tau(S)=15S^0.55 */
     if (pow(slope,0.45)*(height[i][j]-topo[i][j])>taugam*15)
      height[i][j]-=block;
     else change=1;
}

void checkmin(i,j)
int i,j;
{    float min;

     min=10000;
     if (mask[i][j]>0.1) {
     if (height[iup[i]][j]<min) min=height[iup[i]][j];
     if (height[idown[i]][j]<min) min=height[idown[i]][j];
     if (height[i][jup[j]]<min) min=height[i][jup[j]];
     if (height[i][jdown[j]]<min) min=height[i][jdown[j]];
     if (height[i][j]<min) {height[i][j]=min+0.1;
                            checkmin(iup[i],j);
                            checkmin(idown[i],j);
                            checkmin(i,jup[j]);
                            checkmin(i,jdown[j]);} }
}

main()
{    FILE *fp1,*fp2,*fp3;
     int i,j,t; 

     fp1=fopen("greenlandsmallbed","r");
     fp2=fopen("greenlandsmallmask","r");
     fp3=fopen("greenlandsmallsurf","w");
     lattice_size_x=150;
     lattice_size_y=280;
     delta=10000.0;                 /* km */
     taugam=100000.0/(920.0*9.8);   /* tau = 10^5 Pa */
     height=matrix(1,lattice_size_x,1,lattice_size_y);
     mask=matrix(1,lattice_size_x,1,lattice_size_y);
     topo=matrix(1,lattice_size_x,1,lattice_size_y);
     setupgridneighbors();
     heightvec=vector(1,lattice_size_x*lattice_size_y);
     heightvecind=ivector(1,lattice_size_x*lattice_size_y);
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {fscanf(fp1,"%f",&topo[i][j]);
        fscanf(fp2,"%f",&mask[i][j]);
        height[i][j]=topo[i][j];}
     block=100;             /* start at 100 m block size */
     while (block>0.9)      /* stop when block size < 1 m */
      {change=1;
       while (change>0)
        {change=0;
         for (j=1;j<=lattice_size_y;j++)
          for (i=1;i<=lattice_size_x;i++)
           heightvec[(j-1)*lattice_size_x+i]=height[i][j];
         indexx(lattice_size_x*lattice_size_y,heightvec,heightvecind);
         t=0;
         while (t<lattice_size_x*lattice_size_y)
          {t++;
           i=(heightvecind[t])%lattice_size_x;
           if (i==0) i=lattice_size_x;
           j=(heightvecind[t])/lattice_size_x+1;
           if (mask[i][j]>0.1)
            addiflessthanthreshold(i,j);}
         for (j=1;j<=lattice_size_y;j++)
          for (i=1;i<=lattice_size_x;i++)
           checkmin(i,j);}
       block=block/3.33;}
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {if (height[i][j]>0) 
         fprintf(fp3,"%f\n",height[i][j]);
	    else fprintf(fp3,"0\n");}
     fclose(fp1);
     fclose(fp2);
     fclose(fp3);
}  
