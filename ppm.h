#include <stdio.h>

// Output in P6 format, a binary file containing:
// P6
// ncolumns nrows
// Max colour value
// colours in binary format thus unreadable
void save_imageP6(int Width, int Height, char* fname,unsigned char* pixels) {
  FILE *fp;
  const int maxVal=255; 
  
  printf("Saving image %s: %d x %d\n", fname,Width,Height);
  fp = fopen(fname,"wb");
  if (!fp) {
        printf("Unable to open file '%s'\n",fname);
        return;
  }
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", Width, Height);
  fprintf(fp, "%d\n", maxVal);

  for(int j = 0; j < Height; j++) {
		  fwrite(&pixels[j*Width*3], 3,Width,fp);
  }

  fclose(fp);
}

// Output in P3 format, a text file containing:
// P3
// ncolumns nrows
// Max colour value (for us, and usually 255)
// r1 g1 b1 r2 g2 b2 .....
void save_imageP3(int Width, int Height, char* fname,unsigned char* pixels) {
	FILE *fp;
	const int maxVal=255;
	
	printf("Saving image %s: %d x %d\n", fname,Width,Height);
	fp = fopen(fname,"w");
	if (!fp) {
		printf("Unable to open file '%s'\n",fname);
		return;
	}
	fprintf(fp, "P3\n");
	fprintf(fp, "%d %d\n", Width, Height);
	fprintf(fp, "%d\n", maxVal);
	
	int k = 0 ;
	for(int j = 0; j < Height; j++) {
		
		for( int i = 0 ; i < Width; i++)
		{
			fprintf(fp," %d %d %d", pixels[k],pixels[k+1],pixels[k+2]) ;
			k = k + 3 ;
		}
		fprintf(fp,"\n") ;
	}
	fclose(fp);
}

/**

// This main function is meant only to illustrate how to use the save_imageXX functions.
// You should get rid of this code, and just paste the save_imageXX functions into your
// raytrace.cpp code. 
int main() {
	int Width = 128;	// Move these to your setup function. The actual resolution will be
	int Height= 128;	// specified in the input file
    char fname3[20] = "sceneP3.ppm"; //This should be set based on the input file
	char fname6[20] = "sceneP6.ppm"; //This should be set based on the input file
	unsigned char *pixels;
	// This will be your image. Note that pixels[0] is the top left of the image and
	// pixels[3*Width*Height-1] is the bottom right of the image.
    pixels = new unsigned char [3*Width*Height];

	// This loop just creates a gradient for illustration purposes only. You will not use it.
	float scale = 128.0 / (float) Width ;
	int k = 0 ;
	for(int i = 0; i < Height; i++) {
		for (int j = 0; j < Width; j++) {
			int c = (i+j)*scale ;
			pixels[k] = c;
			pixels[k+1] = c;
			pixels[k+2] = c;
			k = k + 3 ;
		}
	}
	save_imageP3(Width, Height, fname3, pixels);
	save_imageP6(Width, Height, fname6, pixels);
}

**/
