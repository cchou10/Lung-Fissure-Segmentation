/****************************************************************/
/* VisX4 program v3tpl                                          */
/* Example program to read in a 16-bit image, perform           */
/* thresholding, and output an 8-bit byte image                 */
/*                                                              */
/****************************************************************/
#include "VisXV4.h"          /* VisX structure include file     */
#include "Vutil.h"           /* VisX utility header files       */

extern char *VisXhist;

char * pname = "v3tpl";

VisXfile_t *VXin,*VXin2,      /* input file structure            */
           *VXout;           /* output file structure           */
VisXelem_t *VXlist,**VXpt;   /* VisX data structure             */
VisXelem_t *mlist;           /* VisX data structure             */

VisX3dim_t  sim;             /* source image structure          */
VisX3dim_t  rim;             /* result image structure          */
VisX3dim_t  mim;             /* mask image structure            */

VisX3dim_t exm;		     /* external image */
VisX3dim_t tm;               /* temp image     */

VisXiinfo_t imginfo;
float xres,yres,zres,ri,rs;

void VX3frameset(VisX3dim_t *is, VisX3dim_t *ir);

VXparam_t par[] = {
    {"if=",    0,   "input file"},
    {"ig=",    0,   "binary mask file"},
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



const int thresh; // threshold for connected components
int count = 0; // count for the SSRA function
int removeQueuex[thresh]; // queue for the small segment removal
int removeQueuey[thresh];

void SSRA(int i, int j, int k); // initialize functions
void SSRAr();
int i,j,k;



int
main(argc, argv)
int argc;
char *argv[];
{
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
    for (k = sim.zlo; k <= sim.zhi; k++) {
    	for (j = sim.ylo; j <= sim.yhi; j++) {
      	    for (i = sim.xlo; i <= sim.xhi; i++) {
		mim.u[k][j][i] = 255;
	    }
        }
    }
  }
  /* Create result image structure */
  VXmake3dim(&rim, VX_PBYTE, sim.bbx, sim.chan);

  /**********************************************************/
  /*******************ALGORITHM: EDIT THIS*******************/

  
	        /*********************** Small Segment Removal Algorithm **********************/
  for (k = sim.zlo; k <= sim.zhi; k++){
   	for (i = sim.ylo; i<= sim.yhi; i++) {
		for (j = sim.xlo; j <= sim.xhi; j++) {
			if (sim.u[k][i][j] != 0) { SSRA(i,j,k); }
			if (count <= thresh) { SSRAr();	}
			count = 0;
		}
  	}


  }

  /**********************END OF ALGORITHM********************/
  /***********************************************************/


  VX3frameset(&sim,&rim);
  
  VXwrite(VXout, rim.list);     /* write data                   */
  VXclose(VXin);                /* close files                  */
  VXclose(VXout);

  exit(0);
}

/* SSRA: recursively counts connected components for size */ 
void SSRA(int i, int j, int k) { 
	int breakf = 0;
	if (count < thresh) {	
		removeQueuex[count] = i;
		removeQueuey[count] = j;
	}
	else {
		memset(removeQueuex, 0, thresh);
		memset(removeQueuey, 0, thresh);
		breakf=1;
	}
	
	if (breakf==0){
	if (tm.u[k][i][j+1] == 255 && sim.u[k][i][j+1] == 255) { 
		count = count + 1;
		SSRA(i, j+1,k); 
	} 
	if (tm.u[k][i][j-1] == 255 && sim.u[k][i][j-1] == 255) { 
		count = count + 1;
		SSRA(i,j-1,k); 
	} 
	if (tm.u[k][i+1][j] == 255 && sim.u[k][i+1][j] == 255) { 
		count = count + 1;
		SSRA(i+1,j,k); 
	} 
	if (tm.u[k][i-1][j] == 255 && sim.u[k][i-1][j] == 255) { 
		count = count + 1;
		SSRA(i-1,j,k); 
	} 
	}
}

/* SSRA: removes the connected components that are under threshold */
void SSRAr() {
	for (i = 0; i < thresh; i++) {
		sim.u[k][removeQueuex[i]][removeQueuey[i]] = 0;
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
