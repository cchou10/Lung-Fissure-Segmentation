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
VisXiinfo_t imginfo;
float xres,yres,zres,ri,rs;

void VX3frameset(VisX3dim_t *is, VisX3dim_t *ir);

VXparam_t par[] = {
    {"if=",    0,   "input file"},
    {"ig=",    0,   "ground truth file"},
    {"of=",    0,   "output file"},
    {"s=",    0,   "fissure size (default: 10000)"},
    {"-v",     0,   "verbose flag"},
    { 0,       0,   0},
};

/* Command line parameters are accessed in code by vars below */
#define  IVAL    par[0].val
#define  MVAL    par[1].val
#define  OVAL    par[2].val
#define  SVAL    par[3].val
#define  VFLAG   par[4].val

int
main(argc, argv)
int argc;
char *argv[];
{
  int i,j,k;
  int xmin,xmax,ymin,ymax,zmin,zmax;

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
  // no mask
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

  /* set history to history of source image file */
  /* this is not necessary since the result image was opened
     before the mask image 
     VXhistset(VXin->list, 1);
  */

  /**************************************************/
  /********EVAL ALGORITHM****************************/
  int size;
  if (SVAL ) {
    size = atoi(SVAL);
  } else {
    size = 10000;
  }
  int basex[size];
  int basey[size];
  int testx[size];
  int testy[size];
  int count = 0;
  int counti = 0;
  int xdiff = 0;
  int ydiff = 0;
  int hyp = 0;
  int distTotal = 0;
  for (k = sim.zlo; k <= sim.zhi; k++) {
    for (j = sim.ylo; j <= sim.yhi; j++) {
      for (i = sim.xlo; i <= sim.xhi; i++) {
        /* Do processing stuff here */
        if (sim.u[k][j][i] == 255) { 
		basex[count] = i;
		basey[count] = j;
	}
	if (mim.u[k][j][i] == 255) { 
	        rim.u[k][j][i] = 255;
		testx[counti] = i;
		testy[counti] = j;
	}
        else {
		rim.u[k][j][i] = 0;
        }
      }
    }
  }
  for (k = 0; k < size; k++) {
	xdiff = (basex[k] - testx[k])*(basex[k] - testx[k]);
	ydiff = (basey[k] - testy[k])*(basey[k] - testy[k]);
	hyp = sqrt(xdiff*xdiff + ydiff*ydiff);
	distTotal = distTotal + hyp;
  }
  distTotal = distTotal / size;
  printf("RMSD is %d", distTotal); 
  /*************************************************/
  VX3frameset(&sim,&rim);
  
  VXwrite(VXout, rim.list);     /* write data                   */
  VXclose(VXin);                /* close files                  */
  VXclose(VXout);

  exit(0);
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
