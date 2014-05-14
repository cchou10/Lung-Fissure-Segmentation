/****************************************************************/
/* VisX4 program v3tpl                                          */
/* Example program to read in a 16-bit image, perform           */
/* thresholding, and output an 8-bit byte image                 */
/*                                                              */
/****************************************************************/
#include "VisXV4.h"          /* VisX structure include file     */
#include "Vutil.h"           /* VisX utility header files       */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern char *VisXhist;

char * pname = "v3tpl";

VisXfile_t *VXin,*VXin2,      /* input file structure            */
           *VXout;           /* output file structure           */
VisXelem_t *VXlist,**VXpt;   /* VisX data structure             */
VisXelem_t *mlist;           /* VisX data structure             */

VisX3dim_t  sim;             /* source image structure          */
VisX3dim_t  tim;             /* temporary image structure       */
VisX3dim_t  tim2;             /* temporary image structure       */
VisX3dim_t  rim;             /* result image structure          */
VisX3dim_t  mim;             /* mask image structure            */
VisXiinfo_t imginfo;
float xres,yres,zres,ri,rs;

void VX3frameset(VisX3dim_t *is, VisX3dim_t *ir);
//cclabel for regional minima
int label=1;
void setlabel(int i,int j,int k, int n);

VXparam_t par[] = {
    {"if=",    0,   "input file"},
    {"ig=",    0,   "binary marker file"},
    {"of=",    0,   "output file"},
    {"th=",    0,   "threshold value (default: 800)"},
    {"-v",     0,   "verbose flag"},
    { 0,       0,   0},
};

/* Command line parameters are accessed in code by vars below */
#define  IVAL    par[0].val
#define  MVAL    par[1].val
#define  OVAL    par[2].val
#define  TVAL    par[3].val
#define  VFLAG   par[4].val


int
main(argc, argv)
int argc;
char *argv[];
{
  int i,j,k;
  int xmin,xmax,ymin,ymax,zmin,zmax;
  int thresh;

  VisXelem_t *vptr = NULL, *mptr = NULL;
  
  VXparse(&argc, &argv, par);    /* parse the command line      */

  VXin  = VXopen(IVAL, 0);       /* open input file             */
  VXout = VXopen(OVAL, 1);       /* open the output file        */

  VXlist = VXread(VXin);         /* read input file             */
  VXgetresinfo( &imginfo);
  VXgetrescale( &imginfo);
  xres = imginfo.xres;
  yres = imginfo.yres;
  zres = imginfo.zres;
  ri = imginfo.ri;
  rs = imginfo.rs;
  if( VFLAG) {
    fprintf(stderr, "img res = %f x %f x %f ri %f rs %f\n",
	    imginfo.xres, imginfo.yres, imginfo.zres, imginfo.ri, imginfo.rs);
  }


  if (TVAL ) {
    thresh = atoi(TVAL);
  } else {
    thresh = 800;
  }

  if(VXNIL == (vptr = VXfind(VXlist, VX_PBYTE))){
    fprintf(stderr, "%s: no acceptable input image found, exiting.\n",pname);
    exit(1);
  }

  /* Initialize input image structure */
  VXset3dim(&sim, vptr, VXin);
  if(sim.chan != 1){
    fprintf(stderr, "%s: Multi-channel images are not supported.\n",pname);
    exit(1);
  }

  //fprintf(stderr,"check\n");

  /* Apply mask if mask image is specified*/
  if ( MVAL ) {
    if ( VFLAG ) {
      fprintf(stderr, "%s: Mask file specified, applying mask...\n",pname);
    }
    VXin2 = VXopen(MVAL, 0);       /* open mask file              */

    /* Read mask file */
    mlist = VXread(VXin2);
    if (VXNIL == (mptr = VXfind(mlist,VX_PBYTE)) ) {
      fprintf(stderr, "%s: Invalid format for mask file.\n",pname);
      exit(1);
    }

    VXset3dim(&mim, mptr, VXin2);

    /* Check if image and mask have same bounding box, warn if not 
       Note that this is not a problem for this program, so we don't
       do anything.                                                */
    if ( (sim.xlo != mim.xlo) || (sim.xhi != mim.xhi) ||
         (sim.ylo != mim.ylo) || (sim.yhi != mim.yhi) ||
         (sim.zlo != mim.zlo) || (sim.zhi != mim.zhi) ) {
      fprintf(stderr, "%s: bounding boxes do not match.\n",pname);
    }
    
    /* Determine what regions overlap */ 
    xmin = (sim.xlo>mim.xlo) ? sim.xlo : mim.xlo;
    ymin = (sim.ylo>mim.ylo) ? sim.ylo : mim.ylo;
    zmin = (sim.zlo>mim.zlo) ? sim.zlo : mim.zlo;
    xmax = (sim.xhi<mim.xhi) ? sim.xhi : mim.xhi;
    ymax = (sim.yhi<mim.yhi) ? sim.yhi : mim.yhi;
    zmax = (sim.zhi<mim.zhi) ? sim.zhi : mim.zhi;
  }
  //no mask
  else {
    VXembed3dim(&mim,&sim,1,1,1,1,1,1);  
    for (k = sim.zlo; k <= sim.zhi; k++) {
    	for (j = sim.ylo; j <= sim.yhi; j++) {
      	    for (i = sim.xlo; i <= sim.xhi; i++) {
		mim.u[k][j][i] = 255;
	    }
        }
    }
  }
  //fprintf(stderr,"checkmask\n");
  /* Create result image structure */

  VXmake3dim(&rim, VX_PBYTE, sim.bbx, sim.chan);

  /*******************************************************/
  /****************WATERSHED ALGORITHM*******************/
  VXembed3dim(&tim,&sim,1,1,1,1,1,1);  
  VXembed3dim(&tim2,&sim,1,1,1,1,1,1);  
  //fprintf(stderr,"checkbegin\n");
  //reset for labelling
  for (k = sim.zlo; k <= sim.zhi; k++) {
    for (j = sim.ylo; j <= sim.yhi; j++) {
      for (i = sim.xlo; i <= sim.xhi; i++) {
          sim.u[k][j][i] = 0;
 	  tim2.u[k][j][i] = 0;
      }
    }
  }
  //fprintf(stderr,"check\n");
  /*
  const int xsize = xmax-xmin+1;
  const int ysize = ymax-ymin+1;
  const int zsize = zmax-ymin+1;
  fprintf(stderr,"check\n");
  float gradim[zsize][ysize][xsize]; //gradient array 
  fprintf(stderr,"checkgrad\n");
  */
  float tolerance = 0.01;
  float gx, gy, gradient;
  /*******Regional Minima Calculation*********************/
  for (k = sim.zlo; k <= sim.zhi; k++) {
    for (j = sim.ylo; j <= sim.yhi; j++) {
      for (i = sim.xlo; i <= sim.xhi; i++) {
      	gx = -tim.u[k][j-1][i-1]+-2*tim.u[k][j][i-1]+-tim.u[k][j+1][i-1]
             +tim.u[k][j-1][i+1]+2*tim.u[k][j][i+1]+tim.u[k][j+1][i+1];  
        gy = -tim.u[k][j-1][i-1]+-2*tim.u[k][j-1][i]+-tim.u[k][j-1][i+1]
             +tim.u[k][j+1][i-1]+2*tim.u[k][j+1][i]+tim.u[k][j+1][i+1];
        //gradim[k][j][i] = sqrt(gx*gx + gy*gy);
        gradient = sqrt(gx*gx + gy*gy);
	//fprintf(stderr,"%f\n",gradient);
	//rim.u[k][j][i] = (int) gradient;
	/*
        if (gradim[k][j][i]<tolerance && mim.u[k][j][i]==255){
		sim.u[k][j][i]=255;
        }*/
  	if (gradient<tolerance && mim.u[k][j][i]==255){
		sim.u[k][j][i]=255;
        }
        else {
        	sim.u[k][j][i]=0;
        } 
      }
    }
  }
  //fprintf(stderr,"checkrmin\n");
  //tim original data
  //sim where regional minima exist
  /************Label Regional Minima*********************/
  for (k = sim.zlo; k <= sim.zhi; k++) {
    for (j = sim.ylo; j <= sim.yhi; j++) {
      for (i = sim.xlo; i <= sim.xhi; i++) {
             if (sim.u[k][j][i]!=0 && tim2.u[k][j][i]==0){
		setlabel(i,j,k,label);
		label=label+1;
	     }
       }
     }
     label=1;
  }
  //fprintf(stderr,"checklabel\n"); 
  //tim2 labeled regional minima
  /***********Watershed Flooding*******************************/ 
  int queuex[10000];//contains x positon of pixel
  int queuey[10000];//contains y position of pixel
  int queuep[10000];//contains pixel value (min prioiritized) 

  int count = 0;
  //fprintf(stderr,"check\n");
  //initial queue
  for (k = sim.zlo; k <= sim.zhi; k++) {
    fprintf(stderr,"check k=%d\n",k);
    for (j = sim.ylo; j <= sim.yhi; j++) {
      //fprintf(stderr,"check j=%d\n",j);
      //fprintf(stderr,"count %d\n",count); 
      for (i = sim.xlo; i <= sim.xhi; i++) {
		//fprintf(stderr,"check i=%d\n",i);
		if (tim2.u[k][j][i] != 0){//labeled
		//add neighbors (and their positions) to queue
 		   if(tim2.u[k][j-1][i]==0){
			tim2.u[k][j-1][i] = tim2.u[k][j][i];
			queuex[count] = i;
			queuey[count] = j-1;
			queuep[count] = tim.u[k][j-1][i];
		   	count=count+1;
		   }
		   if(tim2.u[k][j+1][i]==0){
			tim2.u[k][j+1][i] = tim2.u[k][j][i];
			queuex[count] = i;
			queuey[count] = j+1;
			queuep[count] = tim.u[k][j+1][i];
			count=count+1;
		   }
		   if(tim2.u[k][j][i-1]==0){
			tim2.u[k][j][i-1] = tim2.u[k][j][i];
			queuex[count] = i-1;
			queuey[count] = j;
			queuep[count] = tim.u[k][j][i-1];
			count=count+1;
		   }
		   if(tim2.u[k][j][i+1]==0){
			tim2.u[k][j][i+1] = tim2.u[k][j][i];
			queuex[count] = i+1;
			queuey[count] = j;
			queuep[count] = tim.u[k][j][i+1];
			count=count+1;
		   }
		}       		
      }
    }
    fprintf(stderr,"checkinitial\n");
    //selection sort in ascending order
    int iMin,tempx,tempy,tempp;
    int a,b;
    for (a = 0; a < count; a++){//last count++ means count = sizeof(queue)+1
	iMin = a;
	for (b = a+1; b < count; b++){
		if (queuep[b] < queuep[iMin]) iMin = b;
	}
	if (iMin != a) {
		tempx = queuex[a];
		tempy = queuey[a];
		tempp = queuep[a];
		queuep[a] = queuep[iMin];
		queuex[a] = queuex[iMin];
		queuey[a] = queuey[iMin];
		queuep[iMin] = tempp;
		queuex[iMin] = tempx;
		queuey[iMin] = tempy;
	}
    }
    /*for (a = 0; a < count; a++){
	fprintf(stderr,"queuep%d\n",queuep[a]);
    }*/
    while (count>0){
	fprintf(stderr,"checkcount%d\n",count);
	//add neighbors (and their positions) to queue
 	if(tim2.u[k][queuey[0]-1][queuex[0]]=0){
		//fprintf(stderr,"check1\n");
		tim2.u[k][queuey[0]-1][queuex[0]] = tim2.u[k][queuey[0]][queuex[0]];
		queuex[count] = i;
		queuey[count] = j-1;
		queuep[count] = tim.u[k][j-1][i];
		count=count+1;
	}
	if(tim2.u[k][queuey[0]+1][queuex[0]]==0){
		//fprintf(stderr,"check2\n");
		tim2.u[k][queuey[0]+1][queuex[0]] = tim2.u[k][queuey[0]][queuex[0]];
		queuex[count] = i;
		queuey[count] = j+1;
		queuep[count] = tim.u[k][j+1][i];
		count=count+1;
	}
	if(tim2.u[k][queuey[0]][queuex[0]-1]==0){
		//fprintf(stderr,"check3\n");
		tim2.u[k][queuey[0]][queuex[0]-1] = tim2.u[k][queuey[0]][queuex[0]];
		queuex[count] = i-1;
		queuey[count] = j;
		queuep[count] = tim.u[k][j][i-1];
		count=count+1;
	}
	if(tim2.u[k][queuey[0]][queuex[0]+1]==0){
		//fprintf(stderr,"check4\n");
		tim2.u[k][queuey[0]][queuex[0]+1] = tim2.u[k][queuey[0]][queuex[0]];
		queuex[count] = i+1;
		queuey[count] = j;
		queuep[count] = tim.u[k][j][i+1];
		count=count+1;
	}
  	//fprintf(stderr,"checkadd\n");
	//remove first pixel from queue
	for (a = 0; a < count-1; a++){ 
		queuex[a] = queuex[a+1];
		queuey[a] = queuey[a+1];
		queuep[a] = queuep[a+1];
	}
	count = count - 1;
 	//fprintf(stderr,"checkremove\n");
	//resort queue
	
	for (a = 0; a < count; a++){//last count++ means count = sizeof(queue)+1
		iMin = a;
		for (b = a+1; b < count; b++){
			if (queuep[b] < queuep[iMin]) iMin = b;
		}
		if (iMin != a) {
			tempx = queuex[a];
			tempy = queuey[a];//fprintf(stderr,"count %d\n",count);
			tempp = queuep[a];
			queuep[a] = queuep[iMin];
			queuex[a] = queuex[iMin];
			queuey[a] = queuey[iMin];
			queuep[iMin] = tempp;
			queuex[iMin] = tempx;
			queuey[iMin] = tempy;
		}
   	 } 
	 //fprintf(stderr,"checksort\n");
	
    } 	
    //fprintf(stderr,"count %d\n",count);
    //end of kth slice
    
  }
  fprintf(stderr,"checkwend\n");
  for (k = sim.zlo; k <= sim.zhi; k++) {
    for (j = sim.ylo; j <= sim.yhi; j++) {
      for (i = sim.xlo; i <= sim.xhi; i++) {
	rim.u[k][j][i] = (tim2.u[k][j][i]==0) ? 255 : 0;
      }
    }
  }
  /*******************END OF ALGORITHM*********************/
  /******************************************************/

  VX3frameset(&sim,&rim);
  
  VXwrite(VXout, rim.list);     /* write data                   */
  VXclose(VXin);                /* close files                  */
  VXclose(VXout);

  exit(0);
}


void setlabel(int i,int j, int k, int n){
	tim2.u[k][j][i]=n;

	if (sim.u[k][j][i+1]!=0 && tim2.u[k][j][i+1]==0){
		setlabel(i+1,j,k,n);
	        //fprintf(stderr, "check1\n");
	}
	if (sim.u[k][j][i-1]!=0 && tim2.u[k][j][i-1]==0){
		setlabel(i-1,j,k,n);	
		//fprintf(stderr, "check2\n");
	}
	if (sim.u[k][j-1][i]!=0 && tim2.u[k][j-1][i]==0){
		setlabel(i,j-1,k,n);
		//fprintf(stderr, "check3\n");
	}
	if (sim.u[k][j+1][i]!=0 && tim2.u[k][j+1][i]==0){
		setlabel(i,j+1,k,n);
		//fprintf(stderr, "check4\n");
	}
	
}


void VX3frameset(VisX3dim_t *is, VisX3dim_t *ir) {
/* VX3frameset: inserts frame markers in ir in the same location
   as in is

   This function assumes that both 3D images have the same number
   of images. If not, it will print an error and quit.
*/
    VisXelem_t *src, *dest;     
    VisXelem_t *fptr,*sptr;    /* ptrs into src list */
    VisXelem_t *dptr;          /* ptrs into dst list */

    /* ensure at the beginning of each list */
    src=VXfirst(is->list);
    dest=VXfirst(ir->list);

    if ( VXNIL == (fptr = VXfind(src, VX_FRAME)) ) {
        /* No frames, don't do anything */
        return;
    }
    /* start searching after the first frame marker */
    sptr = src;
    dptr = dest;
    while(VXNIL != (sptr = VXfind(sptr, VX_BBX)) ) {
        if (VXNIL == (dptr = VXfind(dptr, VX_BBX)) ) {
            fprintf(stderr, "Error: Image count not equal!\n");
            exit(1);
        }
        fptr = VXbfind(sptr, VX_FRAME);
        /* Go back one element from BBX to add frame marker */
        dptr = VXaddelem(dptr->prev,VX_FRAME,"",fptr->size);

        fptr = VXfind(fptr, VX_EFRAME);
        /* Set pointer to the image (bbx, then image) */
        /* TODO: Implement variable number of element skipping */
        dptr = dptr->next->next;
        dptr = VXaddelem(dptr, VX_EFRAME, "", fptr->size);
        sptr = sptr->next;
    }
    if (VXNIL != VXfind(dptr, VX_BBX)) {
        fprintf(stderr, "Error: Image count not equal at end!\n");
        exit(1);
    }
}
