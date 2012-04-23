#include <fitsio.h> 

#include "histosize.h"

int saveFITS(char * filename,
	     long nx, 
	     long ny,
	     double array[__ORSA_FITS_N__][__ORSA_FITS_N__]) {
  
  int status = 0;
  long naxis = 2;
  long naxes[2] = {nx,ny};
  
  // create FITS file
  fitsfile * fptr; 
  fits_create_file(&fptr, filename, &status); // use "!filename" to force overwrite
  fits_report_error(stderr, status);
  
  // Create the primary array image
  // fits_create_img(fptr, SHORT_IMG, naxis, naxes, &status);
  fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
  fits_report_error(stderr, status);
  
  /* Write the array of integers to the image */
  long fpixel    = 1;
  long nelements = naxes[0] * naxes[1];
  // fits_write_img(fptr, TSHORT, fpixel, nelements, &array, &status);
  fits_write_img(fptr, TDOUBLE, fpixel, nelements, array, &status);
  fits_report_error(stderr, status);
  
  fits_close_file(fptr, &status);
  fits_report_error(stderr, status);
  
  return status;
}

int main() {
    
    // long nx = __ORSA_FITS_N__;
    // long ny = __ORSA_FITS_N__;
    double array[__ORSA_FITS_N__][__ORSA_FITS_N__];

    const double emptyVal = 0.0;
    for (size_t i=0; i<__ORSA_FITS_N__; ++i) {
        for (size_t j=0; j<__ORSA_FITS_N__; ++j) {
            array[i][j] = emptyVal;
        }
    }
    
    FILE * fp = fopen("histo_2D.out","r");
    char line[4096];
    int i,j;
    double val;
    double minVal = 1.0e99;
    while (fgets(line,4096,fp)) {
        sscanf(line,
               "%i %i %lf",
               &i,
               &j,
               &val);
        array[i][j] = val;
        if (val<minVal) minVal=val;
    }
    
    /* 
       for (size_t i=0; i<__ORSA_FITS_N__; ++i) {
       for (size_t j=0; j<__ORSA_FITS_N__; ++j) {
       if (array[i][j] == emptyVal) array[i][j] = 0.01*minVal;
       }
       }
    */
    
    saveFITS("!histo_2D.fits", __ORSA_FITS_N__, __ORSA_FITS_N__, array);
    
    return 0;
}
