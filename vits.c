/*********************************************************************/
/* vits:          */
/*********************************************************************/

#include "VisXV4.h"          /* VisionX structure include file       */
#include "Vutil.h"           /* VisionX utility header files         */
VisXfile_t *VXin,            /* input file structure                 */
           *VXout;           /* output file structure                */
VisXelem_t *VXlist,*VXptr;   /* VisionX data structure               */
VXparam_t par[] =            /* command line structure               */
{
{  "if=",   0,   " input file, vtpeak: threshold between hgram peaks"},
{  "of=",   0,   " output file "},
{  "-v",    0,   "(verbose) print threshold information"},
{   0,      0,   0} /* list termination */
};
#define  IVAL   par[0].val
#define  OVAL   par[1].val
#define  VFLAG  par[2].val

main(argc, argv)
int argc;
char *argv[];
{

VisXimage_t im;                    /* input image structure          */
int        i,j;                    /* index counters                 */

    int hist[256];                 /* histogram bins                 */
    int thresh;                    /* threshold                      */
    int dist;                   /* minimum distance between maxima   */
    int avg1,avg2;		/* averages of Region1 and Region2   */
    int old1,old2;         /* saved averages from previous iteration */
    int total1=0,total2=0;   		/* for computing average     */
    int number1=0,number2=0;

			     
    VXparse(&argc, &argv, par);    /* parse the command line         */
    VXin  = VXopen(IVAL, 0);       /* open input file                */
    VXout = VXopen(OVAL, 1);       /* open the output file           */

/************ End of Parameter and initialization section ************/

    while((VXlist = VXptr = VXreadframe(VXin)) != VXNIL){ /* every frame */
        VXfupdate(VXout, VXin); /* update global constants */
	/* find next byte image */
        while (VXNIL != (VXptr = VXfind(VXptr, VX_PBYTE)))  { 
            VXsetimage(&im, VXptr, VXin); /* initialize input structure */

/***************** Application specific section **********************/

            /* clear the histogram */
            for (i = 0; i < 256; i++) hist[i] = 0;
 
            /* compute the histogram */
            for (i = im.ylo; i <= im.yhi; i++){
                for (j = im.xlo; j <= im.xhi; j++){
                    hist[im.u[i][j]]++;
		}
	    }
  
            /* compute the threshold */
	    for (i=0; i <256; i++) {
   		number1=number1+hist[i];
		total1=total1+i*hist[i];
            }
	    thresh=total1/number1;
	    fprintf(stderr, "thresh = %d\n",thresh);
	    total1=0;
	    number1=0;
	    for (j=0; j<=10000; j++){
		   for (i=0; i<thresh; i++) {
			total1=total1+i*hist[i];
  			number1=number1+hist[i];
		   }
		   for (i=thresh; i<256; i++){
			total2=total2+i*hist[i];
			number2=number2+hist[i];
		   }
		   if (number1 !=0){
  		   	avg1=total1/number1;
		   }else {
			avg1=thresh/2;
		   }
		   if (number2 !=0){	
		   	avg2=total2/number2;
		   }else {
			avg2=(255-thresh)/2;
		   }
	           thresh=(avg1+avg2)/2;
		   if (avg1==old1 && avg2==old2){
			fprintf(stderr, "avg1 = %d, avg2 = %d\n", avg1, avg2);
			break;
		   }else{
			old1=avg1;
			old2=avg2;
		   }
            }
	    
	    if (VFLAG)
		fprintf(stderr, "thresh = %d\n",thresh);
	    
  
           
            /* apply the threshold */
            for (i = im.ylo; i <= im.yhi; i++) {
                for (j = im.xlo; j <= im.xhi; j++) {
                    if (im.u[i][j] >= thresh) im.u[i][j] = 255;
                    else                      im.u[i][j] = 0;
                }
            }
 
/************** End of the Application specific section **************/

            VXresetimage(&im); /* free the im image structure  */
            VXptr = VXptr->next; /* move to the next image */
        } /* end of every image section */
        VXwriteframe(VXout,VXlist); /* write frame */
        VXdellist(VXlist);          /* delete the frame */
    } /* end of every frame section */
    VXclose(VXin);  /* close files */
    VXclose(VXout);
    exit(0);
}
