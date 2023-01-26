#ifndef IMAGE_H
#define IMAGE_H

#include "specincld.h"

class Image
{
 public:
  Image(double rs, double Te, double rdb, double ts, vector<vector<double> > di, vector<double> rd, vector<double> drd, vector<vector<double> > Qs);
  
  double rstar;
  double Teff;
  double rdust_blowout;
  double tsublimate;
  vector<vector<double> > disk;
  vector<double> rdust;
  vector<double> drdust;
  vector<vector<double> > Qsca;

  int nsizes;
  int nlambda;
  double n0;
  double rstarAU;
  
  vector<vector<vector<double> > > x;
  vector<vector<vector<double> > > y;
  vector<vector<vector<double> > > z;
  vector<vector<vector<double> > > r;
  vector<vector<vector<double> > > dv;
  vector<vector<vector<double> > > cosscattang;

  vector<vector<vector<double> > > disk_imager(vector<vector<vector<double> > > x, vector<vector<vector<double> > > y, vector<vector<vector<double> > > z, vector<vector<vector<double> > > r, vector<vector<vector<double> > > dv, vector<vector<vector<double> > > cosscattang);
};

#endif // IMAGE_H
