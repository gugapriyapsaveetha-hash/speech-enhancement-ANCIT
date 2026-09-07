#include<vector>
#include<complex>
#include<cmath>
#include<algorithm>

const double PI=std::acos(-1.0);
std::vector<std::complex<double>>fft(const std::vector<std::complex<double>>&x)
{
 int N=(int)x.size();
 std::vector<std::complex<double>>X(N);
 for(int k=0;k<N;k++)
 {
  for(int n=0;n<N;n++)
  {
   double angle=-2.0*PI*k*n/N;
   X[k]+=x[n]*std::complex<double>(std::cos(angle),std::sin(angle));
  }
 }
 return X;
}
std::vector<std::complex<double>>ifft(const std::vector<std::complex<double>>&X);
std::vector<double>hanning(int N)
{
 std::vector<double>w(N);
 for(int i=0;i<N;i++)
 {
  w[i]=0.5-0.5*std::cos(2.0*PI*i/(N-1));
 }
 return w;
}
std::vector<double>hamming(int N)
{
 std::vector<double>w(N);
 for(int i=0;i<N;i++)
 {
  w[i]=0.54-0.46*std::cos(2.0*PI*i/(N-1));
 }
 return w;
}
std::vector<double>gaincontrol(const std::vector<double>&Gain,int ConstraintInLength)
{
 int NFFT=(int)Gain.size();
 double meanGain=0.0;
 for(int i=0;i<NFFT;i++)
 {
  meanGain+=Gain[i]*Gain[i];
 }
 meanGain/=NFFT;
 int L2=ConstraintInLength;
 std::vector<double>win=hamming(L2);
 std::vector<std::complex<double>>gain(NFFT);
 for(int i=0;i<NFFT;i++)
 {
  gain[i]=std::complex<double>(Gain[i],0.0);
 }
 std::vector<std::complex<double>>impulse=ifft(gain);
 std::vector<double>ImpulseR(NFFT);
 for(int i=0;i<NFFT;i++)
 {
  ImpulseR[i]=impulse[i].real();
 }
 std::vector<double>ImpulseR2(NFFT,0.0);
 for(int i=0;i<L2/2;i++)
 {
  ImpulseR2[i]=ImpulseR[i]*win[i+L2/2];
 }
 for(int i=NFFT-L2/2;i<NFFT;i++)
 {
  int j=i-(NFFT-L2/2);
  ImpulseR2[i]=ImpulseR[i]*win[j];
 }
 std::vector<std::complex<double>>impulse2(NFFT);
 for(int i=0;i<NFFT;i++)
 {
  impulse2[i]=std::complex<double>(ImpulseR2[i],0.0);
 }
 std::vector<std::complex<double>>temp=fft(impulse2);
 std::vector<double>NewGain(NFFT);
 for(int i=0;i<NFFT;i++)
 {
  NewGain[i]=std::abs(temp[i]);
 }
 double meanNewGain=0.0;
 for(int i=0;i<NFFT;i++)
 {
  meanNewGain+=NewGain[i]*NewGain[i];
 }
 meanNewGain/=NFFT;
 if(meanNewGain>0)
 {
  double factor=std::sqrt(meanGain/meanNewGain);
  for(int i=0;i<NFFT;i++)
  {
   NewGain[i]*=factor;
  }
 }
 return NewGain;
}

void WienerNoiseReduction(const std::vector<double>&ns,int fs,double IS,std::vector<double>&esTSNR,std::vector<double>&esHRNR)
{
 int l=(int)ns.size();
 int wl=(int)std::floor(0.020*fs);
 int NFFT=2*wl;
 std::vector<double>hanwin=hanning(wl);
 if(IS<=0)
 {
  IS=10*wl;
 }
 std::vector<double>nsum(NFFT,0.0);
 int count=0;
 int silence=(int)IS;
 for(int m=0;m<=silence-wl;m++)
 {
  std::vector<std::complex<double>>nwin(NFFT);
  for(int i=0;i<wl;i++)
  {
   nwin[i]=std::complex<double>(ns[m+i]*hanwin[i],0.0);
  }
  std::vector<std::complex<double>>temp=fft(nwin);
  for(int i=0;i<NFFT;i++)
  {
   nsum[i]+=std::norm(temp[i]);
  }
  count++;
 }
 if(count==0)
 {
  esTSNR=ns;
  esHRNR=ns;
  return;
 }
 std::vector<double>d(NFFT);
 for(int i=0;i<NFFT;i++)
 {
  d[i]=nsum[i]/count;
  if(d[i]<=0)
  {
   d[i]=1e-12;
  }
 }
 double SP=0.25;
 double normFactor=1.0/SP;
 int overlap=(int)std::floor((1.0-SP)*wl);
 int offset=wl-overlap;
 int max_m=(int)std::floor((double)(l-NFFT)/offset);

 if(max_m<0)
 {
  esTSNR=ns;
  esHRNR=ns;
  return;
 }
 std::vector<double>oldmag(NFFT,0.0);
 std::vector<double>news(l,0.0);
 std::vector<std::vector<double>>phasea(max_m+1,std::vector<double>(NFFT,0.0));
 std::vector<std::vector<double>>xmaga(max_m+1,std::vector<double>(NFFT,0.0));
 std::vector<std::vector<double>>tsnra(max_m+1,std::vector<double>(NFFT,0.0));
 std::vector<std::vector<double>>newmags(max_m+1,std::vector<double>(NFFT,0.0));
 double alpha=0.98;
 for(int m=0;m<=max_m;m++)
 {
  int begin=m*offset;
  std::vector<std::complex<double>>winy(NFFT);
  for(int i=0;i<wl;i++)
  {
   winy[i]=std::complex<double>(hanwin[i]*ns[begin+i],0.0);
  }
  std::vector<std::complex<double>>ffty=fft(winy);
  std::vector<double>phasey(NFFT);
  std::vector<double>magy(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   phasey[i]=std::arg(ffty[i]);
   magy[i]=std::abs(ffty[i]);
   phasea[m][i]=phasey[i];
   xmaga[m][i]=magy[i];
  }
  std::vector<double>postsnr(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   postsnr[i]=((magy[i]*magy[i])/d[i])-1.0;
   postsnr[i]=std::max(postsnr[i],0.1);
  }
  std::vector<double>eta(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   eta[i]=alpha*((oldmag[i]*oldmag[i])/d[i])+(1.0-alpha)*postsnr[i];
  }
  std::vector<double>newmag(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   newmag[i]=(eta[i]/(eta[i]+1.0))*magy[i];
  }
  std::vector<double>tsnr(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   tsnr[i]=(newmag[i]*newmag[i])/d[i];
  }
  std::vector<double>Gtsnr(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   Gtsnr[i]=tsnr[i]/(tsnr[i]+1.0);
   tsnra[m][i]=Gtsnr[i];
  }
  Gtsnr=gaincontrol(Gtsnr,NFFT/2);
  std::vector<std::complex<double>>output(NFFT);
 for(int i=0;i<NFFT;i++)
  {
   newmag[i]=Gtsnr[i]*magy[i];
   newmags[m][i]=newmag[i];
   output[i]=newmag[i]*std::complex<double>(std::cos(phasey[i]),std::sin(phasey[i]));
   oldmag[i]=newmag[i];
  }
  std::vector<std::complex<double>>temp=ifft(output);
  for(int i=0;i<NFFT;i++)
  {
   int pos=begin+i;
   if(pos<l)
   {
    news[pos]+=temp[i].real()/normFactor;
   }
  }
 }
 esTSNR=news;
 std::vector<double>newharm(l);
 for(int i=0;i<l;i++)
 {
  newharm[i]=std::max(esTSNR[i],0.0);
 }
 news.assign(l,0.0);
 for(int m=0;m<=max_m;m++)
 {
  int begin=m*offset;
  std::vector<double>nharm(wl);
  for(int i=0;i<wl;i++)
  {
   nharm[i]=hanwin[i]*newharm[begin+i];
  }
  std::vector<std::complex<double>>harmonic(NFFT);
  for(int i=0;i<wl;i++)
  {
   harmonic[i]=std::complex<double>(nharm[i],0.0);
  }
  std::vector<std::complex<double>>ffth=fft(harmonic);
  std::vector<double>ffthMag(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   ffthMag[i]=std::abs(ffth[i]);
  }
  std::vector<double>snrham(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   snrham[i]=(tsnra[m][i]*(newmags[m][i]*newmags[m][i])+(1.0-tsnra[m][i])*(ffthMag[i]*ffthMag[i]))/d[i];
  }
  std::vector<double>newgain(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   newgain[i]=snrham[i]/(snrham[i]+1.0);
  }
  newgain=gaincontrol(newgain,NFFT/2);
  std::vector<std::complex<double>>output(NFFT);
  for(int i=0;i<NFFT;i++)
  {
   double newmag=newgain[i]*xmaga[m][i];
   double phase=phasea[m][i];
   output[i]=newmag*std::complex<double>(std::cos(phase),std::sin(phase));
  }
  std::vector<std::complex<double>>temp=ifft(output);
  for(int i=0;i<NFFT;i++)
  {
   int pos=begin+i;
   if(pos<l)
   {
    news[pos]+=temp[i].real()/normFactor;
   }
  }
 }
 esHRNR=news;
}