#include<vector>
#include<complex>
#include<cmath>
#include<algorithm>
#include<numeric>
const double FX_PI=std::acos(-1.0);
std::vector<std::complex<double>>fx_fft(const std::vector<std::complex<double>>&x)
{
 int N=static_cast<int>(x.size());
 std::vector<std::complex<double>>X(N);
 for(int k=0;k<N;k++)
 {
  for(int n=0;n<N;n++)
  {
   double angle=-2.0*FX_PI*k*n/N;
   X[k]+=x[n]*std::complex<double>(std::cos(angle),std::sin(angle));
  }
 }
 return X;
}
std::vector<std::complex<double>>fx_ifft(const std::vector<std::complex<double>>&X)
{
 int N=static_cast<int>(X.size());
 std::vector<std::complex<double>>x(N);
 for(int n=0;n<N;n++)
 {
  for(int k=0;k<N;k++)
  {
   double angle=2.0*FX_PI*k*n/N;
   x[n]+=X[k]*std::complex<double>(std::cos(angle),std::sin(angle));
  }
  x[n]/=static_cast<double>(N);
 }
 return x;
}
std::vector<double>fx_hamming(int N)
{
 std::vector<double>w(N);
 for(int i=0;i<N;i++)
 {
  w[i]=0.54-0.46*std::cos(2.0*FX_PI*i/N);
 }
 return w;
}
double fx_median(std::vector<double>x)
{
 if(x.empty())
  return 0.0;
 std::sort(x.begin(),x.end());
 int n=static_cast<int>(x.size());
 if(n%2==0)
  return 0.5*(x[n/2-1]+x[n/2]);
 return x[n/2];
}

std::vector<double>fx_smooth(const std::vector<double>&x,int n)
{
 int N=static_cast<int>(x.size());
 if(N==0)
  return x;
 if(n<1)
  return x;

 std::vector<double>y(N,0.0);
 int half=n/2;
 for(int i=0;i<N;i++)
 {
  int a=std::max(0,i-half);
  int b=std::min(N-1,i+half);
  double sum=0.0;
  for(int j=a;j<=b;j++)
   sum+=x[j];
  y[i]=sum/(b-a+1);
 }

 return y;
}

double fx_gaussian2(const std::vector<double>&x,const std::vector<double>&mean,const std::vector<double>&var)
{
 double det=var[0]*var[1]-var[2]*var[2];
 if(det<=0.0)
  det=1e-12;

 double a=x[0]-mean[0];
 double b=x[1]-mean[1];
 double q=(var[1]*a*a-2.0*var[2]*a*b+var[0]*b*b)/det;
 return -0.5*(2.0*std::log(2.0*FX_PI)+std::log(det)+q);
}

double fx_logsum(const std::vector<double>&v)
{
 if(v.empty())
  return -1e100;
 double mx=*std::max_element(v.begin(),v.end());
 double sum=0.0;
 for(double x:v)
  sum+=std::exp(x-mx);
 return mx+std::log(sum);
}

double fx_gaussmixp(const std::vector<double>&x,const std::vector<std::vector<double>>&means,const std::vector<std::vector<double>>&vars,const std::vector<double>&weights)
{
 std::vector<double>p(means.size());
 for(size_t i=0;i<means.size();i++)
 {
  std::vector<double>v=vars[i];
  double lp=fx_gaussian2(x,means[i],v);
  p[i]=std::log(std::max(weights[i],1e-15))+lp;
 }
 return fx_logsum(p);
}

std::vector<double>fx_findpeaks(const std::vector<double>&x,int minDistance)
{
 std::vector<double>pos;
 int N=static_cast<int>(x.size());
 if(N<3)
  return pos;
 int last=-minDistance;
 for(int i=1;i<N-1;i++)
 {
  if(x[i]>x[i-1] && x[i]>=x[i+1])
  {
   if(i-last>=minDistance)
   {
    pos.push_back(static_cast<double>(i));
    last=i;
   }
   else if(!pos.empty())
   {
    int old=static_cast<int>(pos.back());
    if(x[i]>x[old])
    {
     pos.back()=static_cast<double>(i);
     last=i;
    }
   }
  }
 }
 return pos;
}

void fx_filterbank(const std::vector<double>&freq,double lowFreq,double highFreq,std::vector<double>&cf,std::vector<std::vector<double>>&trans)
{
 int N=static_cast<int>(freq.size());
 if(N<2)
  return;
 int bins=static_cast<int>(std::round(std::log(highFreq/lowFreq)/std::log(2.0)*24.0));
 bins=std::max(20,bins);
 cf.resize(bins);
 for(int i=0;i<bins;i++)
 {
  double r=static_cast<double>(i)/(bins-1);
  cf[i]=lowFreq*std::pow(highFreq/lowFreq,r);
 }

 trans.assign(bins,std::vector<double>(N,0.0));
 for(int j=0;j<bins;j++)
 {
  double left;
  double centre;
  double right;
  if(j==0)
   left=cf[0];
  else
   left=0.5*(cf[j-1]+cf[j]);
  centre=cf[j];
  if(j==bins-1)
   right=cf[bins-1];
  else
   right=0.5*(cf[j]+cf[j+1]);
  for(int k=0;k<N;k++)
  {
   double f=freq[k];
   if(f>=left && f<=centre && centre>left)
    trans[j][k]=(f-left)/(centre-left);
   else if(f>centre && f<=right && right>centre)
    trans[j][k]=(right-f)/(right-centre);
  }
 }
}

std::vector<double>fx_filter_spectrum(const std::vector<double>&spectrum,const std::vector<double>&filter)
{
 int N=static_cast<int>(spectrum.size());
 int M=static_cast<int>(filter.size());

 std::vector<double>result(N,0.0);
 int centre=M/2;
 for(int i=0;i<N;i++)
 {
  for(int j=0;j<M;j++)
  {
   int k=i+j-centre;
   if(k>=0 && k<N)
    result[i]+=spectrum[k]*filter[j];
  }
 }
 return result;
}

std::vector<double>fx_logfilter(const std::vector<double>&cf)
{
 int N=static_cast<int>(cf.size());
 std::vector<double>f(N);
 if(N==0)
  return f;
 int ini=0;
 while(ini<N && cf[ini]<=3.0*cf[0])
  ini++;
 if(ini>=N)
  ini=0;
 double scale=cf[ini];
 for(int i=0;i<N;i++)
  f[i]=cf[i]/scale;
 std::vector<double>result;
 for(double v:f)
 {
  if(v>0.5 && v<10.5)
   result.push_back(v);
 }
 return result;
}
void fxpefac(const std::vector<double>&s,int fs,double tinc,double fres,std::vector<double>&fx,std::vector<double>&tx,std::vector<double>&pv)
{
 double fstep=5.0;
 double fmax=4000.0;
 double fbanklo=10.0;
 if(fs<=0 || s.empty())
 {
  fx.clear();
  tx.clear();
  pv.clear();
  return;
 }
 if(tinc<=0.0)
  tinc=0.01;
 int winlen2=static_cast<int>(std::round(1.81*fs/fres));
 if(winlen2<16)
  winlen2=16;

 int hop=static_cast<int>(std::round(tinc*fs));
 if(hop<1)
  hop=1;

 int nfft=1;
 while(nfft<4*winlen2)
  nfft*=2;

 std::vector<double>w=fx_hamming(winlen2);
 int nframes=0;
 if(static_cast<int>(s.size())>=winlen2)
  nframes=1+(static_cast<int>(s.size())-winlen2)/hop;

 if(nframes<=0)
 {
  fx.clear();
  tx.clear();
  pv.clear();
  return;
 }

 std::vector<double>freq;
 for(int k=0;k<=nfft/2;k++)
 {
  double f=static_cast<double>(k)*fs/nfft;
  if(f>=0.0 && f<=fmax)
   freq.push_back(f);
 }

 int nf=static_cast<int>(freq.size());
 std::vector<std::vector<double>>MIX(nframes,std::vector<double>(nf,0.0));
 tx.resize(nframes);
 for(int frame=0;frame<nframes;frame++)
 {
  int start=frame*hop;
  std::vector<std::complex<double>>in(nfft);
  for(int i=0;i<winlen2;i++)
  {
   if(start+i<static_cast<int>(s.size()))
    in[i]=std::complex<double>(s[start+i]*w[i],0.0);
  }

  std::vector<std::complex<double>>spec=fx_fft(in);
  for(int k=0;k<nf;k++)
   MIX[frame][k]=std::norm(spec[k]);

  tx[frame]=(start+winlen2/2.0)/fs;
 }

 std::vector<double>cf;
 std::vector<std::vector<double>>trans;
 fx_filterbank(freq,fbanklo,freq.back(),cf,trans);
 int nbands=static_cast<int>(cf.size());
 if(nbands==0)
 {
  fx.clear();
  tx.clear();
  pv.clear();
  return;
 }

 std::vector<std::vector<double>>O(nframes,std::vector<double>(nbands,0.0));
 for(int i=0;i<nframes;i++)
 {
  for(int j=0;j<nbands;j++)
  {
   double value=0.0;
   for(int k=0;k<nf;k++)
    value+=MIX[i][k]*trans[j][k];
   O[i][j]=value;
  }
 }

 std::vector<double>ltass(nbands,1.0);
 for(int j=0;j<nbands;j++)
 {
  double f=cf[j];
  if(f<20.0)
   ltass[j]=1.0;
  else
   ltass[j]=1.0/std::sqrt(f);
 }

 double maxlt=*std::max_element(ltass.begin(),ltass.end());
 for(double&v:ltass)
  v/=maxlt;
 std::vector<double>auxf(nbands+1);
 if(nbands==1)
 {
  auxf[0]=cf[0];
  auxf[1]=cf[0];
 }
 else
 {
  auxf[0]=cf[0];
  for(int i=1;i<nbands;i++)
   auxf[i]=0.5*(cf[i-1]+cf[i]);
  auxf[nbands]=cf.back();
 }

 std::vector<double>binwidth(nbands);
 for(int i=0;i<nbands;i++)
  binwidth[i]=std::max(auxf[i+1]-auxf[i],1e-12);
 for(int i=0;i<nframes;i++)
 {
  for(int j=0;j<nbands;j++)
   O[i][j]*=binwidth[j];
 }

 std::vector<double>eltass(nbands,0.0);
 for(int j=0;j<nbands;j++)
 {
  for(int i=0;i<nframes;i++)
   eltass[j]+=O[i][j];
  eltass[j]/=nframes;
 }

 eltass=fx_smooth(eltass,21);
 for(int j=0;j<nbands;j++)
 {
  double alpha=ltass[j]/std::max(eltass[j],1e-12);
  for(int i=0;i<nframes;i++)
   O[i][j]*=alpha;
 }

 std::vector<double>sca;
 int ini=0;
 while(ini<nbands && cf[ini]<=3.0*cf[0])
  ini++;

 if(ini>=nbands)
  ini=nbands/4;
 double scale=cf[ini];
 for(int i=0;i<nbands;i++)
 {
  double value=cf[i]/scale;
  if(value>0.5 && value<10.5)
   sca.push_back(value);
 }

 std::vector<double>filh(sca.size());
 double pefact=1.8;
 for(size_t i=0;i<sca.size();i++)
  filh[i]=1.0/(pefact-std::cos(2.0*FX_PI*sca[i]));
 double avg=0.0;
 for(double v:filh)
  avg+=v;
 if(!filh.empty())
  avg/=filh.size();
 for(double&v:filh)
  v-=avg;

 double maxfil=0.0;
 for(double v:filh)
  maxfil=std::max(maxfil,std::abs(v));
 if(maxfil>0.0)
 {
  for(double&v:filh)
   v/=maxfil;
 }

 std::vector<std::vector<double>>B(nframes,std::vector<double>(nbands,0.0));

 for(int i=0;i<nframes;i++)
  B[i]=fx_filter_spectrum(O[i],filh);
 int numopt=3;

 std::vector<std::vector<double>>ff(nframes,std::vector<double>(numopt,0.0));

 std::vector<std::vector<double>>amp(nframes,std::vector<double>(numopt,0.0));

 double flim1=60.0;
 double flim2=400.0;
 int first=0;
 int last=nbands-1;

 while(first<nbands && cf[first]<=flim1)
  first++;
 while(last>=0 && cf[last]>=flim2)
  last--;

 for(int i=0;i<nframes;i++)
 {
  std::vector<double>peakpos;
  for(int j=std::max(1,first);
      j<std::min(last,nbands-1);
      j++)
  {
   if(B[i][j]>B[i][j-1] &&
      B[i][j]>=B[i][j+1])
   {
    peakpos.push_back(static_cast<double>(j));
   }
  }

  std::sort(peakpos.begin(),peakpos.end(),[&](double a,double b){return B[i][static_cast<int>(a)]>B[i][static_cast<int>(b)];});
  int count=std::min(numopt,static_cast<int>(peakpos.size()));

  for(int j=0;j<count;j++)
  {
   int p=static_cast<int>(peakpos[j]);
   ff[i][j]=cf[p];
   amp[i][j]=B[i][p];
  }
 }

 std::vector<std::vector<double>>meansU=
 {
  {13.38533,0.4199435},
  {12.23505,0.1496836},
  {12.76646,0.2581733},
  {13.69822,0.6893078},
  {9.804372,0.02786567},
  {11.03848,0.07711229}
 };

 std::vector<double>weightsU=
 {
  0.1461799,
  0.3269458,
  0.2632178,
  0.02331986,
  0.06360947,
  0.1767271
 };

 std::vector<std::vector<double>>varsU=
 {
  {0.4575519,0.002619074,0.002619074,0.01262138},
  {0.7547719,0.008568089,0.008568089,0.001933864},
  {0.5770533,0.003561592,0.003561592,0.00527957},
  {0.3576287,0.01388739,0.01388739,0.04742106},
  {0.9049906,0.01033191,0.01033191,0.0001887114},
  {0.637969,0.009936445,0.009936445,0.0007082946}
 };

 std::vector<std::vector<double>>meansV=
 {
  {15.36667,0.8961554},
  {13.52718,0.4809653},
  {13.95531,0.8901121},
  {14.56318,0.6767258},
  {14.59449,1.190709},
  {13.11096,0.2861982}
 };

 std::vector<double>weightsV=
 {
  0.1391365,
  0.221577,
  0.2214025,
  0.1375109,
  0.1995124,
  0.08086066
 };

 std::vector<std::vector<double>>varsV=
 {
  {0.196497,-0.002605404,-0.002605404,0.05495016},
  {0.6054919,0.007776652,0.007776652,0.01899244},
  {0.5944617,0.0485788,0.0485788,0.03511229},
  {0.3871268,0.0292966,0.0292966,0.02046839},
  {0.3377683,0.02839657,0.02839657,0.04756354},
  {1.00439,0.03595795,0.03595795,0.006737475}
 };

 pv.resize(nframes);
 std::vector<std::vector<double>>vuvfea(nframes,std::vector<double>(2,0.0));
 for(int i=0;i<nframes;i++)
 {
  double power=0.0;
  for(double v:O[i])
   power+=v;
 power/=nbands;
  power=std::max(power,1e-12);
  double peakSum=amp[i][0]+amp[i][1]+amp[i][2];
  vuvfea[i][0]=std::log(power);
  vuvfea[i][1]=1e-3*peakSum/(power+1.75e5);
  double pru=fx_gaussmixp(vuvfea[i],meansU,varsU,weightsU);
  double prv=fx_gaussmixp(vuvfea[i],meansV,varsV,weightsV);
  double z=pru-prv;

  if(z>50.0)
   pv[i]=0.0;
  else if(z<-50.0)
   pv[i]=1.0;
  else
   pv[i]=1.0/(1.0+std::exp(z));
 }

 std::vector<std::vector<double>>cost(nframes,std::vector<double>(numopt,0.0));
 std::vector<std::vector<int>>prev(nframes,std::vector<int>(numopt,0));
 std::vector<double>medfx(nframes,0.0);

 double w1=1.0;
 double w2=0.8250;
 double w3=0.01868;
 double w4=0.006773;
 double w5=98.9;
 double w6=-0.4238;
 double txinc=tinc;

 int inmf=static_cast<int>(std::round(2.0/txinc));
 inmf=std::max(1,inmf);
 std::vector<double>firstFreq;
 for(int i=0;i<std::min(inmf,nframes);i++)
 {
  if(ff[i][0]>0.0 && pv[i]>0.6)
   firstFreq.push_back(ff[i][0]);
 }

 double mf=fx_median(firstFreq);
 if(mf==0.0)
 {
  firstFreq.clear();
  for(int i=0;i<std::min(inmf,nframes);i++)
  {
   if(ff[i][0]>0.0 && pv[i]>0.5)
    firstFreq.push_back(ff[i][0]);
  }
  mf=fx_median(firstFreq);
 }
 if(mf==0.0)
 {
  firstFreq.clear();
  for(int i=0;i<std::min(inmf,nframes);i++)
  {
   if(ff[i][0]>0.0 && pv[i]>0.4)
    firstFreq.push_back(ff[i][0]);
  }
  mf=fx_median(firstFreq);
 }
 if(mf==0.0)
 {
  firstFreq.clear();
  for(int i=0;i<std::min(inmf,nframes);i++)
  {
   if(ff[i][0]>0.0 && pv[i]>0.3)
    firstFreq.push_back(ff[i][0]);
  }
  mf=fx_median(firstFreq);
 }
 medfx[0]=mf;
 for(int j=0;j<numopt;j++)
 {
  if(amp[0][j]==0.0)
   cost[0][j]=w5;
  else
  {
   double mx=*std::max_element(amp[0].begin(),amp[0].end());
   cost[0][j]=-amp[0][j]/std::max(mx,1e-12);
  }
 }
 double dffact=2.0/txinc;
 for(int i=1;i<nframes;i++)
 {
  if(i>inmf)
  {
   std::vector<double>temp;
   for(int k=i-inmf;k<=i;k++)
   {
    if(ff[k][0]>0.0 && pv[k]>0.6)
     temp.push_back(ff[k][0]);
   }
   mf=fx_median(temp);
   if(mf==0.0)
  {
    temp.clear();
   for(int k=i-inmf;k<=i;k++)
    {
     if(ff[k][0]>0.0 && pv[k]>0.5)
      temp.push_back(ff[k][0]);
    }
    mf=fx_median(temp);
   }
   if(mf==0.0)
   {
    temp.clear();
    for(int k=i-inmf;k<=i;k++)
    {
     if(ff[k][0]>0.0 && pv[k]>0.4)
      temp.push_back(ff[k][0]);
    }
    mf=fx_median(temp);
   }
   if(mf==0.0)
   {
    temp.clear();
    for(int k=i-inmf;k<=i;k++)
    {
     if(ff[k][0]>0.0 && pv[k]>0.3)
      temp.push_back(ff[k][0]);
    }
    mf=fx_median(temp);
   }
  }
  medfx[i]=mf;
  for(int j=0;j<numopt;j++)
  {
   double bestCost=1e100;
   int bestPrev=0;
   for(int k=0;k<numopt;k++)
   {
    double a=ff[i][j];
    double b=ff[i-1][k];
    double df=0.0;
    if(a+b>0.0)
     df=dffact*(a-b)/(a+b);
    double costdf=w3*std::min((df-w6)*(df-w6),w4);
    double value=cost[i-1][k]+costdf;
    if(value<bestCost)
    {
     bestCost=value;
     bestPrev=k;
    }
   }
   cost[i][j]=bestCost;
   double costf=0.0;
   if(mf!=0.0)
    costf=std::abs(ff[i][j]-mf)/mf;
   double mx=*std::max_element(amp[i].begin(),amp[i].end());
   double camp;
   if(amp[i][j]==0.0)
    camp=w5;
   else
    camp=-amp[i][j]/std::max(mx,1e-12);
   prev[i][j]=bestPrev;
 cost[i][j]+=w2*costf+w1*camp;
  }
 }
 fx.assign(nframes,0.0);
 int best=0;
 if(nframes>0)
 {
  best=static_cast<int>(
       std::min_element(cost.back().begin(),cost.back().end())-cost.back().begin());
  fx[nframes-1]=ff[nframes-1][best];
  for(int i=nframes-1;i>=1;i--)
  {
   best=prev[i][best];
   fx[i-1]=ff[i-1][best];
  }
 }
}