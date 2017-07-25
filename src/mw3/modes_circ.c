/*--------------------------- Commande MegaWave -----------------------------*/
/* mwcommand
  name = {modes_circ};
  version = {"1.0"};
  author = {"Agnes Desolneux"};
  function = {"modes detection in a Fsignal for uniform circular prior"};
  usage = {
'e':[eps=0.0]->eps  "-log10(max. number of false alarms), default 0",
'm'->m        "to have  maximal eps-meaningful modes and gaps, else only max meaningful intervals",
'f': nfa<-nfa "Flist of nfa of detected modes", 
'a'->a        "alternative definition of nfa, using intervals defined by two points",
in->in           "input Fsignal",
out<-modes_circ "output Fsignal of meaningful modes"
          };
*/
/*-- MegaWave - Copyright (C) 1994 Jacques Froment. All Rights Reserved. --*/

/*J.D : correction janvier 2007  certains creux n'�taient pas propag�s dans ifcreux*/


#include <stdio.h>
#include <math.h>
#include "mw3.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* fonction entropie relative */
float entrop(double x, double y)
{
    if (x<=y) return(0.0);
    else if (x==1.0) return (-(float) log10(y));
    else return (float) (x*log10(x/y)+(1.0-x)*log10((1.0-x)/(1.0-y)));
}


float *maxsup(float *b,
              int N)
{
    float  *c;
    int i,j;
    c=(float *) calloc(N*N , sizeof (float));


    c[N-1]=b[N-1];
    for(i=1;i<N;++i)
        c[i*N +i-1] = b[i*N +i-1];

    c[1*N+N-1]= MAX(b[1*N+N-1], MAX(c[N-1],c[1*N+0]));

    for(i=2;i<N;++i)
        for (j=i-2;j>=0;--j)
            c[i*N +j] = MAX(b[i*N +j],MAX(c[(i-1)*N+j],c[i*N+j+1])) ;

    for(j=N-2;j>=0;--j)
        c[j]=MAX(b[j],MAX(c[j+1],c[(N-1)*N+j]));

    for(i=2 ; i<N ; i++)
        c[i*N +N-1] = MAX(MAX(c[(i-1)*N +N-1] , c[i*N +0] ), b[i*N +N-1]) ;

    for (i=1 ; i<N-1 ; i++)
        for (j=N-2 ; j>= i ; j--)
            c[i*N +j] = MAX(MAX(c[(i-1)*N +j] , c[i*N +j+1] ), b[i*N +j]) ;

    return c;
    
}


float *maxinf(float *b,
              int N)
{
    float  *c;
    int i,j;

    c=(float *) calloc( N*N, sizeof (float));



    for (i=N-1 ; i>= 0 ; i--)
        for (j=i ; j<N ; j++)
        {if (j==i)
                c[i*N +j] = b[i*N +j];
            else
                c[i*N +j] = MAX(MAX(c[(i+1)*N +j] , c[i*N +j-1] ), b[i*N +j]) ;}

    c[(N-1)*N] = MAX(c[0] , MAX(c[(N-1)*N +N-1], b[(N-1)*N]));
    for (j=1; j<N-1; ++j)
        c[(N-1)*N +j] = MAX(b[(N-1)*N +j], MAX(c[j], c[(N-1)*N +j-1])) ;


    for (i=N-2;i>0;i--)
        for (j=0 ; j< i ; ++j)
        { if (j==0)
                c[i*N +j] = MAX(b[i*N +j], MAX( c[i*N+N-1] , c[(i+1)*N +j]));
            else
                c[i*N +j] = MAX(b[i*N +j], MAX( c[i*N+j-1] , c[(i+1)*N +j]));
        }

    return c;
}

/* ne garder que les intervalles qui ne contiennent pas de creux significatif */
float *ifcreux(float *b,
               float *cr,
               int N,
               float s)
{
    float  *c;
    int i,j;

    c=(float *) calloc( N*N, sizeof (float));

    for (i=0;i<N;++i)
        for (j=0;j<N; ++j)
        {if( cr[i*N +j]>=s)
            {c[i*N +j]=1;
                /*printf("creux significatif, debut=%d, fin=%d\n", i,j);*/ }
            else c[i*N +j]=0;}

    for (i=N-2;i>=0;i--)
        for (j=i+1;j<N; ++j)
            c[i*N +j] = c[i*N +j] + c[i*N+j-1] + c[(i+1)*N +j];

    c[(N-1)*N] = c[0] + c[(N-1)*N +N-1] + c[(N-1)*N];
    for (j=1; j<N-1; ++j)
        c[(N-1)*N +j] = c[(N-1)*N +j] + c[j] + c[(N-1)*N +j-1] ;

    /*rajout 09/01/07, certains creux n'�taient pas propag�s*/
    for (i=N-2; i>=0; i--)
        c[i*N] = c[i*N] + c[i*N+N-1] + c[(i+1)*N] ;


    for (i=N-2;i>0;i--)
        for (j=1 ; j< i ; ++j)
            c[i*N +j] = c[i*N +j] + c[i*N+j-1] + c[(i+1)*N +j];

    for (i=0;i<N;++i)
        for (j=0;j<N; ++j)
        {if (c[i*N +j] > 0) c[i*N +j]= 0.0;
            else c[i*N +j] = b[i*N +j];}

    return c;
}



/********** MAIN MODULE *******************************/
Fsignal modes_circ(Fsignal in,char *m,char *a,double *eps,Flist nfa)
{
    Fsignal out=0;
    float nbtotal, *H, *Hbis, *Hsup, *Hinf, seuil, *creux;
    double *r ,  *p;
    int i,j,size, nb_int,k;

    size = in->size;
    /* out=mw_new_fsignal(); */
    out=mw_change_fsignal(out,size);
    mw_clear_fsignal(out,0.0);


    H= (float *) calloc ( (size*size) , sizeof (float));
    r= (double *) calloc ( (size*size) , sizeof (double));
    p= (double *) calloc ( (size*size) , sizeof (double));
    creux = (float *) calloc ( (size*size) , sizeof (float));

    if(nfa) nfa = mw_change_flist(nfa,size,0,1);


    /* integrate signal */
    for (i=1;i<size;i++) in->values[i]+=in->values[i-1];

    nbtotal = in->values[size-1];
    /*ON NE FAIT RIEN S'IL Y A MOINS DE 2 ELEMENTS DANS LE CAS DE L'OPTION (a)*/  if ((a)&&(nbtotal<=2)) return(out);

    seuil = ((float) log10( (double) size*(size-1))+ *eps)/nbtotal ;
    /* printf("seuil=%f\n" , seuil); */
    if (a)  seuil = ((float) log10( (double) (nbtotal-1)*nbtotal )+ *eps)/(nbtotal-2) ;


    /* normalisation du  signal */
    for (i=0;i<size;i++) in->values[i]=in->values[i] /nbtotal;

    for(i=0;i<size;++i)
        for(j=0;j<size;++j)
        { if (j<i)
            {p[i*size +j]=((double) size -i+j+1)/((double) size) ;
                r[i*size +j]=(double) in->values[j] + 1 - in->values[i-1];}
            else
            { p[i*size +j]= (j-i+1)/( (double) size);
                if (i==0)
                    r[i*size +j]=(double) in->values[j];
                else
                    r[i*size +j]=(double) in->values[j]-  in->values[i-1] ;}
        }



    for(i=0;i<size;++i)
        for(j=0;j<size;++j)
        {H[i*size +j]= entrop(r[i*size +j],p[i*size +j]);
            creux[i*size +j]= entrop(1-r[i*size +j], 1-p[i*size +j]);

            if(a)
            { if (r[i*size +j]*(float)nbtotal<=2.) H[i*size+j]= 0.;
                else H[i*size +j]= entrop((r[i*size +j]*(float)nbtotal-2)/(nbtotal-2),p[i*size +j]);
                if ((1-r[i*size +j])*(float)nbtotal<=2.) creux[i*size+j]= 0.;
                else creux[i*size +j]= entrop(((1-r[i*size +j])*(float)nbtotal-2)/(nbtotal-2),1-p[i*size +j]);
            }
        }

    printf("%f,%f,%f,seuil=%f\n",creux[1*size +11],creux[0*size +11],creux[1*size +12],seuil);

    if (m)
    {
        Hbis=ifcreux(H,creux,size,seuil);
        Hsup=maxsup(Hbis,size);
        Hinf=maxinf(Hbis,size);
        nb_int =0;
        for(i=0;i<size;++i)
            for(j=0;j<size;++j)
            { printf("nb_int=%d,  Hbis=%f , nfa=%f,debut=%d, fin=%d\n", nb_int, Hbis[i*size +j],(float) log10((double) (size-1)*size) -(float)nbtotal*Hbis[i*size+j],i,j);
                if((Hsup[i*size +j]<= Hbis[i*size +j]) && (Hinf[i*size +j]<= Hbis[i*size +j]) && (Hbis[i*size+j]>= seuil))
                { nb_int++;

                    nfa->values[nfa->size++]=(float) log10((double) (size-1)*size) -(float)nbtotal*Hbis[i*size+j];
                    if(i<=j)
                    {for(k=i;k<j+1;++k)
                            out->values[k]=(float) nb_int;}
                    else
                    {for(k=i;k<size;++k)
                            out->values[k]=(float) nb_int;
                        for(k=0;k<j+1;++k)
                            out->values[k]=(float) nb_int;}
                }
            }
    }
    else
    {
        Hsup=maxsup(H,size);
        Hinf=maxinf(H,size);
        nb_int=0;
        for(i=0;i<size;++i)
            for(j=0;j<size;++j)
            {printf("i=%d j=%d r=%f p=%f H=%f, seuil = %f\n", i,j,r[i*size +j] , p[i*size +j] , H[i*size +j], seuil);
                if((Hsup[i*size +j]<= H[i*size +j]) && (Hinf[i*size +j]<= H[i*size +j]) && (H[i*size+j]>= seuil))
                {nb_int++;

                    nfa->values[nfa->size++]=(float) log10((double) (size-1)*size) -(float)nbtotal*Hbis[i*size+j];
                    if(a) nfa->values[nfa->size++]=(float) log10((double) (nbtotal-1)*nbtotal) -(float)(nbtotal-2)*Hbis[i*size+j];

                    if(i<=j)
                    {for(k=i;k<j+1;++k)
                            out->values[k]=(float) nb_int;}
                    else
                    {for(k=i;k<size;++k)
                            out->values[k]=(float) nb_int;
                        for(k=0;k<j+1;++k)
                            out->values[k]=(float) nb_int;}
                }
            }
    }
    return(out);
}







