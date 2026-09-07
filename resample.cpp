#include<vector>
#include<cmath>

std::vector<double>resample(const std::vector<double>&x,int oldFs,int newFs)
{
 if(oldFs==newFs)
 {
  return x;
 }
 int oldSize=(int)x.size();
 if(oldSize==0)
 {
  return {};
 }

 int newSize=(int)std::round((double)oldSize*newFs/oldFs);
 std::vector<double>y(newSize);
 double ratio=(double)(oldSize-1)/(newSize-1);
 for(int i=0;i<newSize;i++)
 {
  double position=i*ratio;
  int index=(int)std::floor(position);
  double fraction=position-index;
  if(index>=oldSize-1)
  {
   y[i]=x[oldSize-1];
  }
  else
  {
   y[i]=x[index]*(1.0-fraction)+x[index+1]*fraction;
  }
 }
 return y;
}