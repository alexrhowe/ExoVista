#include "specincld.h"
#include "Image.h"

Image::Image(double rs, double Te, double rdb, double ts, vector<vector<double> > di, vector<double> rd, vector<double> drd, vector<vector<double> > Qs)
{
  rstar = rs;
  Teff = Te;
  rdust_blowout = rdb;
  tsublimate = ts;
  disk = di;
  rdust = rd;
  drdust = drd;
  Qsca = Qs;

  nsizes = rdust.size();
  nlambda = Qsca[0].size();
  if(nlambda==0) nlambda = 1;

  // Here is the normalization factor that gives the proper
  // V band surface brightness for 1 zodi solar system twin.
  // I play pretty fast and loose with units in this routine, 
  // so this includes all of the unit conversions that I skip.
  
  // (Based on a calibration of exozodiacal fluxes output by the code.)
  n0 = 4.739e-8; // normalization using Hedman & Stark SPF, multiplied by average SPF(90) of the exoVista SPF distribution

  rstarAU = rstar * 0.00465047;
}

vector<vector<vector<double> > > Image::disk_imager(vector<vector<vector<double> > > x, vector<vector<vector<double> > > y, vector<vector<vector<double> > > z, vector<vector<vector<double> > > r, vector<vector<vector<double> > > dv, vector<vector<vector<double> > > cosscattang)
{
  // Reference: d = ('n', 'longnode', 'i', 'nzodis', 'r', 'dror', 'rinner', 'eta', 'hor', 'g0', 'g1', 'g2', 'w0', 'w1', 'w2')
  
  double pi = 4.*atan(1.);
  
  vector<int> dims(3);
  dims[0] = x.size();
  dims[1] = x[0].size();
  dims[2] = x[0][0].size();

  vector<vector<vector<double> > > flux(dims[0], vector<vector<double> >(dims[1], vector<double>(nlambda,0)));
  vector<vector<vector<double> > > t(dims[0], vector<vector<double> >(dims[1], vector<double>(dims[2],0)));
  vector<vector<vector<double> > > n(dims[0], vector<vector<double> >(dims[1], vector<double>(dims[2],0)));

  for(int i=0; i<disk.size(); i++){
    vector<double> d = disk[i];
    if(d[3] > 0){
      double g0 = d[9];
      double g1 = d[10];
      double g2 = d[11];
      double w0 = d[12];
      double w1 = d[13];
      double w2 = d[14];
	    
      double r0 = d[4];
      double dr = d[5] * r0;
      double rinner = d[6];
      double oneodr = 1./dr;
      double oneor0pdr = 1./(r0+dr);
      double oneor0mdr = 1./(r0-dr);
      
      for(int j=0; j<dims[0]; j++){
	for(int k=0; k<dims[1]; k++){
	  for(int l=0; l<dims[2]; l++){
	    
	    t[j][k][l] = Teff * sqrt(sqrt(1.0-0.5) * rstarAU / (2*r[j][k][l])) * 1.5; // 1.5x factor is a fudge factor for grain super-heating
            // Calculate the scattering phase function for this disk component (multi-component HG)
            double temp0 = 1.+g0*g0 - 2.*g0*cosscattang[j][k][l];
            double temp1 = 1.+g1*g1 - 2.*g1*cosscattang[j][k][l];
            double temp2 = 1.+g2*g2 - 2.*g2*cosscattang[j][k][l];
            double pfunc = w0*(1.-g0*g0)/(4*pi)/pow(temp0,1.5) + w1*(1.-g1*g1)/(4*pi)/pow(temp1,1.5) + w2*(1.-g2*g2)/(4*pi)/pow(temp2,1.5);
            
            // Now calculate disk density at each point
            // n = np.zeros(dims)        # all points set to zero
            n[j][k][l] = (t[j][k][l] < tsublimate) * (r[j][k][l] > rstarAU) * n0 * d[3]; // grid outside of sublimation radius set to nzodis
            
            // Multiply by differential volume
            n[j][k][l] *= dv[j][k][l];
            
            // Modify by the variation w/ z
            n[j][k][l] *= exp(-0.5 * pow(z[j][k][l]/(d[8] * r[j][k][l]),2));
            n[j][k][l] /= (d[8] * sqrt(2*pi)) * r[j][k][l]; // normalize by the integral over z to maintain same density vs r
            
            // The Gaussian ring component
            if(abs(r[j][k][l]-r0) < dr) n[j][k][l] *= exp(-0.5 * pow((r[j][k][l]-r0)*oneodr,2));
            
	    // The blowout component
            // Make sure it's continuous at +dr from r0
            if(r[j][k][l] > r0+dr){
	      double temp = r[j][k][l] * oneor0pdr;
	      n[j][k][l] *= exp(-0.5) / pow(temp,1.5);
	    }
            
            // The inner truncation radius
            if(r[j][k][l] <= rinner){
	      double rorinner = r[j][k][l]/rinner;
	      n[j][k][l] *= pow(rorinner,3);
	    }
            
            // Multiply by phase function and 1/r^2     
            n[j][k][l] *= pfunc/(4.*pi*r[j][k][l]*r[j][k][l]);
	  }
	}
      }
      
      for(int isize=0; isize<nsizes; isize++){
	double grainfactors = pi * pow(rdust[isize],-3.5+2.0) * drdust[isize];
	double tempeta = d[7] * rdust[isize]/rdust_blowout; // this creates a radial color gradient in disk
	for(int j=0; j<dims[0]; j++){
	  for(int k=0; k<dims[1]; k++){
	    
	    double tempflux = 0.;
	    for(int l=0; l<dims[2]; l++){
	      
	      double temp = 1. - sqrt(r[j][k][l] * oneor0mdr);
	      double tempn = n[j][k][l] * grainfactors; // multiply by Donanyi size distribution and cross section
	      
	      // Finally take care of the inward decay component
              // It can change depending on grain size
	      // Make sure it's continuous at -dr from r0
	      if(r[j][k][l] < r0-dr) tempn *= exp(-0.5) / (1. + 4.*tempeta * temp);
	      
	      if(!isfinite(tempn)){
		printf("tempn NaN!");
		vector<vector<vector<double> > > fluxnull(dims[0], vector<vector<double> >(dims[1], vector<double>(nlambda,-1)));
		return fluxnull;
	      }
	      
	      // Now we can integrate along z-axis
	      tempflux += tempn;
	    }
	    for(int m=0; m<nlambda; m++){
	      flux[j][k][m] += tempflux * Qsca[isize][m];
	    }
	  }
	}
      } // End isize loop
    } // End if
  } // End i loop
  return flux;
}
