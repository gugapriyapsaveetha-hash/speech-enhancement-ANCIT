#include<vector>
#include<random>
#include<cmath>
std::vector<double>addAWGN(const std::vector<double>&x,double SNR)
{
 std::vector<double>y=x;
 double signalPower=0.0;
 for(int i=0;i<static_cast<int>(x.size());i++)
 {
  signalPower+=x[i]*x[i];
 }
 if(x.size()>0)
 {
  signalPower/=static_cast<double>(x.size());
 }
 double snrLinear=std::pow(10.0,SNR/10.0);
 double noisePower=0.0;
 if(snrLinear>0.0)
 {
  noisePower=signalPower/snrLinear;
 }
 double noiseStd=std::sqrt(noisePower);
 std::random_device rd;
 std::mt19937 generator(rd());
 std::normal_distribution<double>distribution(0.0,noiseStd);
 for(int i=0;i<static_cast<int>(x.size());i++)
 {
  y[i]=x[i]+distribution(generator);
 }
 return y;
}