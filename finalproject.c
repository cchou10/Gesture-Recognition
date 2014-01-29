/*********************************************************************/ 
/* finalproject: */
/* Does image analysis on various hands to identify */
/* the number of fingers they are holding up */ 
/*********************************************************************/
#include "VisXV4.h" /* VisionX structure include file */ 
#include "Vutil.h" /* VisionX utility header files */ 
VisXfile_t *VXin, /* input file structure */ 
*VXout; /* output file structure */ 
VisXelem_t *VXlist,*VXptr; /* VisionX data structure */ 
VXparam_t par[] = /* command line structure */ 
{ 
{ "if=", 0, " input file, vtpeak: threshold between hgram peaks"}, 
{ "of=", 0, " output file "}, 
{ "d=", 0, " min dist between hgram peaks (default 10)"}, 
{ "-v", 0, "(verbose) print threshold information"}, 
{ "x=", 0, " center of mass x value"},
{ "y=", 0, " center of mass y value"},
{ 0, 0, 0} /* list termination */ }; 

#define IVAL par[0].val 
#define OVAL par[1].val 
#define DVAL par[2].val 
#define VFLAG par[3].val
#define X par[4].val
#define Y par[5].val

/* pulling these from within main method (to use in method setlabel) */ 
VisXimage_t im; // i/o image structure 
VisXimage_t tm; // temporary image structure 
int l = 1; // index counters 
int range; // range value 
int first = 0;

void setlabel(int i, int j, int n); // initialize function

main(argc, argv) 
int argc; 
char *argv[]; 
{
	VisXimage_t im; /* input image structure */ 
	int i,j; /* index counters */
	int hist[256]; /* histogram bins */ 
	int thresh; /* threshold */ 
	int dist; 
	VXparse(&argc, &argv, par); /* parse the command line */ 
	VXin = VXopen(IVAL, 0); /* open input file */ 
	VXout = VXopen(OVAL, 1); /* open the output file */

	/************ Parameter and initialization section *******************/
	dist = 10; /* default dist */ 
	if (DVAL) dist = atoi(DVAL); /* if d= was specified, get value */ 
	if (dist < 0 || dist > 255) { 
		fprintf(stderr, "d= must be between 0 and 255\nUsing d=10\n"); 
		dist = 10; 
	}
	
	/************ threshold the image to find the hand ************/
	while((VXlist = VXptr = VXreadframe(VXin)) != VXNIL){ /* every frame */ 
		VXfupdate(VXout, VXin); /* update global constants */ 
		/* find next byte image */ 
		while (VXNIL != (VXptr = VXfind(VXptr, VX_PBYTE))) { 
			VXsetimage(&im, VXptr, VXin); /* initialize input structure */
			int total = 0; /* used to calculate the total pixel value */
			int number = 0; 
			int maxloops = 0; /* maximinum number of iterations to prevent infinitelooping */ 
			int vglow; /* takes average of pixel values */ 
			int avghi; 
			int preAvglow; /* saves value to compare to current average value */ 
			int preAvghi;
			
			/* clear the histogram */ 
			for (i = 0; i < 256; i++) hist[i] = 0;
			
			/* compute the histogram */ 
			for (i = im.ylo; i <= im.yhi; i++) { 
				for (j = im.xlo; j <= im.xhi; j++) { 
					hist[im.u[i][j]]++; 
					total = total + im.u[i][j]; 
				} 
			}
			
			/* compute the threshold */ 
			thresh = total / ((im.yhi-im.ylo)*(im.xhi-im.xlo));
			
			/* iterates through until convergence (or if it has too many iterations) */ 
			while (maxloops <= 1000 && avglow != preAvglow && avghi != preAvghi) { 
				preAvghi = avghi; 
				preAvglow = avglow; 
				for (i = 0; i < thresh; i++) { 
					total = total + hist[i]*i; 
				} 
				avglow = total / thresh; 
				total = 0; 
				for (i = thresh; i <= 255; i++) { 
					total = total + hist[i]*i; 
				} 
				avghi = total / (255 - thresh); total = 0;
				thresh = round(avghi + avglow) / 2);
				maxloops = maxloops + 1; }
				if(VFLAG) fprintf(stderr, "Average 1 = %d Average 2 = %d thresh = %d\n", avglow, avghi, thresh);
				
				/* apply the threshold for a light background image */ 
				for (i = im.ylo; i <= im.yhi; i++) { 
					for (j = im.xlo; j <= im.xhi; j++) { 
						if (im.u[i][j] >= thresh) im.u[i][j] = 0; 
						else im.u[i][j] = 255; 
					} 
				}	
	/*********************** close images *********************************************/				
				VXresetimage(&im); /* free the im image structure */
				VXptr = VXptr->next; /* move to the next image */ 
			} /* end of every image section */ 	
			
	/*********************** find center of mass using moments ************************/
			
	
	/*********************** add region growing to determine number of fingers ********/
			i = 0; /* index counters */
			j = 0;
			range = atoi(RVAL); 
			if (VXNIL != (VXptr = VXfind(VXlist, VX_PBYTE))){ /* find image */ 
				VXsetimage(&im, VXptr, VXin); /* initialize input structure */ 
				VXembedimage(&tm,&im,1,1,1,1); /* image structure with border */
				
				/* sets entire original image to be zeros */ 
				for (i = im.ylo; i <= im.yhi; i++) { 
					for (j = im.xlo; j <= im.xhi; j++) { 
						im.u[i][j] = 0; 
					} 
				}
				
				/* goes through image and sets the label */ 
				for (i = im.ylo; i<= im.yhi; i++) { 
					for (j = im.xlo; j <= im.xhi; j++) { 
						if (tm.u[i][j] != 0 && im.u[i][j] == 0) { 
							first = tm.u[i][j]; 
							if (PFLAG) { 
								setlabel(i, j, first); 
							} 
							else { 
								setlabel(i, j, l); 
							} 
							if (l == 255) { 
								l = 1; 
							} 
							else { 
								l = l + 1; 
							} 
						} 
					} 
				}
				VXwrite(VXout, VXlist); /* write data */ 
			} 
			else { 
				fprintf(stderr, "vtemp: no byte image data in input file\n"); 
				exit(-1); 
			} 
			
			VXwriteframe(VXout,VXlist); /* write frame */ 
			VXdellist(VXlist); /* delete the frame */ 
		} /* end of every frame section */ 
		VXclose(VXin); /* close files */ 
		VXclose(VXout); 
		exit(0); 
	}
	
	/* setlabel: recursively sets label based on requirements documented in lab writeup */ 
	void setlabel(int i, int j, int l) { 
		im.u[i][j] = l; 
		if (tm.u[i][j+1] != 0 && im.u[i][j+1] == 0 && tm.u[i][j+1]-first < range) { setlabel(i, j+1, l); } 
		if (tm.u[i][j-1] != 0 && im.u[i][j-1] == 0 && tm.u[i][j-1]-first < range) { setlabel(i,j-1,l); } 
		if (tm.u[i+1][j] != 0 && im.u[i+1][j] == 0 && tm.u[i+1][j]-first < range) { setlabel(i+1,j,l); } 
		if (tm.u[i-1][j] != 0 && im.u[i-1][j] == 0 && tm.u[i-1][j]-first < range) { setlabel(i-1,j,l); } 
	}
