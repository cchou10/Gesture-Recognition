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
{ "th=", 0, "threshold of number of pixels inside bounding box"},
{ "l=", 0, "length of bounding box"},
{ "w=", 0, "width of bounding box"},
{ "-v", 0, "(verbose) print threshold information"}, 
{ "x=", 0, " center of mass x value"},
{ "y=", 0, " center of mass y value"},
{ 0, 0, 0} /* list termination */ }; 

#define IVAL par[0].val 
#define OVAL par[1].val 
#define DVAL par[2].val 
#define TVAL par[3].val
#define LVAL par[4].val
#define WVAL par[5].val
#define VFLAG par[6].val
#define X par[7].val
#define Y par[8].val

/* pulling these from within main method (to use in method setlabel) */ 
VisXimage_t im; // i/o image structure 
VisXimage_t tm; // temporary image structure
VisXimage_t mm; //temporary mask image to hold bounding boxes
int l = 1; // index counters 
int range; // range value 
int first = 0;

void setlabel(int i, int j, int n); // initialize function

main(argc, argv) 
int argc; 
char *argv[]; 
{
	
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

    int length; /*threshold in bounding box*/
	if (LVAL) length = atoi(LVAL); /* if d= was specified, get value */
    if (length < im.ylo || length > im.yhi) {
        fprintf(stderr, "l= must be between %d and %d\nUsing l=10\n",im.ylo,im.yhi);
        length = 10;
    }
    int width; /*threshold in bounding box*/
	if (WVAL) width = atoi(WVAL); /* if d= was specified, get value */
    if (dist < im.xlo || dist > im.xhi) {
        fprintf(stderr, "w= must be between %d and %d\nUsing w=10\n",im.xlo,im.xhi);
        width = 10;
    }
    int fthresh; /*threshold of count in bounding box*/
	if (TVAL) fthresh = atoi(TVAL); /* if d= was specified, get value */
    if (dist < 0 || dist > length*width) {
        fprintf(stderr, "th= must be between 0 and %d\nUsing th=20\n",length*width);
        fthresh = 30;
    }
    /************ Bounding Box Location Specification *******************/
    int xpos[5]; /*x positions of bounding boxes relative to COM*/
    int ypos[5]; /*y positions of bounding boxes relative to COM*/
    
    
	/************ threshold the image to find the hand ************/
	while((VXlist = VXptr = VXreadframe(VXin)) != VXNIL){ /* every frame */ 
		VXfupdate(VXout, VXin); /* update global constants */ 
		/* find next byte image */ 
		while (VXNIL != (VXptr = VXfind(VXptr, VX_PBYTE))) { 
			VXsetimage(&im, VXptr, VXin); /* initialize input structure */
            VXembedimage(&tm,&im,1,1,1,1); /* image structure with border */
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
						if (tm.u[i][j] >= thresh) im.u[i][j] = 0;
						else tm.u[i][j] = 255;
					} 
				}	
	/*********************** close images *********************************************/				
				VXresetimage(&im); /* free the im image structure */
				VXptr = VXptr->next; /* move to the next image */ 
			} /* end of every image section */ 	
			
	/*********************** find center of mass using moments ************************/
            VXparse(&argc, &argv, par); /* parse the command line */
            VXin = VXopen(IVAL, 0); /* open input file */
            VXout = VXopen(OVAL, 1); /* open the output file */
            VXlist = VXread(VXin); /* read file */
        
            VXembedimage(&tm,&im,1,1,1,1); /* image structure with border */
            total = 0;
            int xtotal = 0;
            int ytotal = 0;
            int xcom;
            int ycom;
            for (i = im.ylo; i <= im.yhi; i++) {
                for (j = im.xlo; j <= im.xhi; j++) {
                    total = total + tm.u[i][j];
                    xtotal = xtotal + j*tm.u[i][j];
                    ytotal = ytotal + i*tm.u[i][j];
                }
            }
            xcom = xtotal/total;
            ycom = ytotal/total;
    /************ Bounding Box Mask Construction **************************************/
        VXembedimage(&mm,&im,1,1,1,1); /* image structure with border */
        for (i = im.ylo; i <= im.yhi; i++) {
            for (j = im.xlo; j <= im.xhi; j++) {
                mm.u[i][j] = 0;
            }
        }
        int lhalf = length/2, whalf = width/2;
        if (length%2 == 0) lhalf = length/2 -1;
        if (width%2 ==0) whalf = width/2 -1;
        for (int k = 0; k < 5; k++){
            for (i = -length/2; i <= length/2; i++){
                for (j = -width/2; j <= width/2; j++){
                    mm.u[ycom+y[k]+i][xcom+x[k]+j] = 255;
                }
                    
            }
        }
        l = 1;
        for (j = im.xlo; j <= im.xhi; j++) {
            for (i = im.ylo; i <= im.yhi; i++) {
                if (mm.u[i][j] != 0 && mm.u[i][j] == 0) {
                    
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
    /************ Finger Detection ****************************************************/
        int f1 = 0, f2 = 0, f3 = 0, f4 = 0, f5 = 0; /* counts for each finger in order thumb, pointer, middle, ring, pinky */
        for (i = im.ylo; i <= im.yhi; i++){
            for (j = im.xlo; j <= im.xhi; j++){
                if (tm.u[i][j] != 0 && mm.u[i][j] != 0){
        
                    switch (mm.u[i][j]){
                        case 1:
                            f1++;
                            break;
                        case 2:
                            f2++;
                            break;
                        case 3:
                            f3++;
                            break;
                        case 4:
                            f4++;
                            break;
                        case 5:
                            f5++;
                            break;
                    }
                }
            }
        }
        bool thumb = FALSE, pointer = FALSE, middle = FALSE, ring = FALSE, pinky = FALSE;
        if (f1 > fthresh) thumb = TRUE;
        if (f2 > fthresh) pointer = TRUE;
        if (f3 > fthresh) middle = TRUE;
        if (f4 > fthresh) ring = TRUE;
        if (f5 > fthresh) pinky = TRUE;
        
        /******** ASL NUMBER DEFINITION **********/
        int ASLNUM = 0;
        if (!thumb && pointer && !middle && !ring && !pinky) ASLNUM = 1;
        if (!thumb && pointer && middle && !ring && !pinky) ASLNUM = 2;
        if (thumb && pointer && middle && !ring && !pinky) ASLNUM = 3;
        if (!thumb && pointer && middle && ring && pinky) ASLNUM = 4;
        if (thumb && pointer && middle && ring && pinky) ASLNUM = 5;
        if (!thumb && pointer && middle && ring && !pinky) ASLNUM = 6;
        if (!thumb && pointer && middle && !ring && pinky) ASLNUM = 7;
        if (!thumb && pointer && !middle && ring && pinky) ASLNUM = 8;
        if (!thumb && !pointer && middle && ring && pinky) ASLNUM = 9;
        if (thumb && !pointer && !middle && !ring && !pinky) ASLNUM = 10;
        
	/*********************** add region growing to determine number of fingers ********/
			i = 0; /* index counters */
			j = 0;
			VXparse(&argc, &argv, par); /* parse the command line */ 
			VXin = VXopen(IVAL, 0); /* open input file */ 
			VXout = VXopen(OVAL, 1); /* open the output file */ 
			VXlist = VXread(VXin); /* read file */ 
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