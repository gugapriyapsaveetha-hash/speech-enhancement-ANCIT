#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>

const double PI=std::acos(-1.0);
double princ(double angle);
void WienerNoiseReduction(const std::vector<double>&ns,int fs,double IS,std::vector<double>&esTSNR,std::vector<double>&esHRNR);

std::vector<double> addAWGN(const std::vector<double>&x,double SNR);
std::vector<double> resample(const std::vector<double>&x,int fs,int fsr);

bool readWavFile(const std::string&filename,std::vector<double>&data,int&fs);
void writeWavFile(const std::string&filename,const std::vector<double>&data,int fs);
void fxpefac(const std::vector<double>&x,int fs,double tinc,double fres,std::vector<double>&f0,std::vector<double>&tv,std::vector<double>&pv);

std::vector<double> fftReal(const std::vector<double>&x)
{
    int N=static_cast<int>(x.size());
    std::vector<double>result(N);
    for(int k=0;k<N;k++)
    {
     std::complex<double>sum(0.0,0.0);
     for(int n=0;n<N;n++)
     {
     double angle=-2.0*PI*k*n/N;
     sum+=x[n]*std::complex<double>(std::cos(angle),std::sin(angle));
     }
     result[k]=std::abs(sum);
    }
    return result;
}
std::vector<std::complex<double>> fft(const std::vector<double>&x)
{
    int N=static_cast<int>(x.size());
    std::vector<std::complex<double>>X(N);
    for(int k=0;k<N;k++)
    {
     for(int n=0;n<N;n++)
     {
     double angle=-2.0*PI*k*n/N;
        X[k]+=x[n]*std::complex<double>(
        std::cos(angle),std::sin(angle));
        }
    }
   return X;
}

std::vector<std::complex<double>> ifft(
    const std::vector<std::complex<double>>&X)
{
    int N=static_cast<int>(X.size());
    std::vector<std::complex<double>>x(N);
    for(int n=0;n<N;n++)
    {
     for(int k=0;k<N;k++)
     {
      double angle=2.0*PI*k*n/N;
        x[n]+=X[k]*std::complex<double>(
        std::cos(angle),std::sin(angle));
      }
       x[n]/=static_cast<double>(N);
    }
    return x;
}

std::vector<double> hammingPeriodic(int N)
{
    std::vector<double>w(N);
    for(int i=0;i<N;i++)
    {
     w[i]=0.54-0.46*
        std::cos(2.0*PI*i/N);
    }
    return w;
}
int main()
{
    std::string filename;
    std::cout<<"\n\nFilename: ";
    std::cin>>filename;
    std::vector<double>x;
    int fs=0;
    if(!readWavFile(filename,x,fs))
    {
     std::cout<<"Could not read audio file.\n";
     return 1;
    }
    int fsr=16000;
    double M_default=0.032*fsr;
    double overlap=0.875;
    double H=(1.0-overlap)*M_default;
    int K=6;
    if(fs!=fsr)
    {
        x=resample(x,fs,fsr);
        fs=fsr;
    }
    double SNR=30.0;
    x=addAWGN(x,SNR);
    double fres=1.81*fs/M_default;
    double tinc=H/static_cast<double>(fs);

    std::vector<double>f0;
    std::vector<double>tv;
    std::vector<double>pv;
    fxpefac(x,fs,tinc,fres,f0,tv,pv);
    int nframes=static_cast<int>(tv.size());
    std::vector<int>voiced(nframes,0);
    for(int i=0;i<nframes;i++)
    {
     if(pv[i]>=0.9)
     {
        voiced[i]=1;
     }
   }
    double IS=M_default;
    int i=1;
    while(i<nframes && voiced[i]!=1)
    {
     IS+=H;
     i++;
   }
    std::vector<double>x_noisy=x;
    std::vector<double>esTSNR;
    std::vector<double>x_filtered;

    WienerNoiseReduction(x,fs,IS,esTSNR,x_filtered);
    x=x_filtered;
    double min_overlap=0.73;
    int M_min=static_cast<int>(
        std::round(H/(1.0-min_overlap)));
    int M_max=1024;

    std::vector<int>M(nframes,0);
    for(int l=0;l<nframes;l++)
    {
     if(voiced[l]==1)
     {
      if(f0[l]>0.0)
       {
       M[l]=static_cast<int>(
       std::round(K*fs/f0[l]));
       }
     else
       {
        M[l]=static_cast<int>(
        std::round(M_default));
       }
           if(M[l]<M_min)
            {
            M[l]=M_min;
           }
            else if(M[l]>M_max)
            {
            M[l]=M_max;
            }
        }
      else
       {
        M[l]=static_cast<int>(
        std::round(M_default));
      }
    }
    std::vector<std::vector<std::complex<double>>>X(nframes);
    std::vector<int>Mh(nframes);
    std::vector<std::vector<double>>w(nframes);
    double offx=M_default/2.0;
    for(int l=0;l<nframes;l++)
    {
        int start;
        if(M[l]%2==0)
       {
        Mh[l]=M[l]/2;
         start=static_cast<int>(
          std::round(offx))-Mh[l];
       }
        else
        {
         Mh[l]=(M[l]-1)/2;
          start=static_cast<int>(
            std::round(offx))-Mh[l]-1;
        }
        if(start<0)
        {
        start=0;
        }
       w[l]=hammingPeriodic(M[l]);
        std::vector<double>xf(M[l],0.0);
        for(int j=0;j<M[l];j++)
        {
            int pos=start+j;
            if(pos>=0 &&
             pos<static_cast<int>(x.size()))
           {
            xf[j]=x[pos]*w[l][j];
          }
       }
        X[l]=fft(xf);
        offx+=H;
    }
    std::vector<std::vector<std::complex<double>>>Y=X;
    int nharm_prev=0;
    std::vector<double>phaY_prev;
    for(int l=1;l<nframes;l++)
    {
        if(voiced[l]==1 && f0[l]>0.0)
       {
        int ML=M[l];
         std::vector<double>magY(ML);
         std::vector<double>phaY(ML);
         for(int k=0;k<ML;k++)
          {
            magY[k]=std::abs(X[l][k]);
            phaY[k]=std::arg(X[l][k]);
            }
            int nbands;
            if(ML%2==0)
            {
            nbands=1+ML/2;
            }
            else
            {
           nbands=(1+ML)/2;
            }
            phaY.resize(nbands);
            int nharm=static_cast<int>(std::floor(static_cast<double>(fs)/(2.0*f0[l])));
            std::vector<int>bh(nharm);
            for(int r=1;r<=nharm;r++)
            {
             bh[r-1]=r*K;
            }
            if(voiced[l-1]!=0)
            {
            if(nharm_prev>=nharm)
             {
               for(int j=0;j<nharm;j++)
               {
                int index=bh[j];
                  if(index<nbands &&index<static_cast<int>(phaY_prev.size()))
                        {
                            double wh=2.0*PI*(j+1)*f0[l]/fs;
                            phaY[index]=princ(phaY_prev[index]+wh*H);
                        }
                    }
                }
                else
                {
                  for(int r=1;r<=nharm_prev;r++)
                   {
                    int index=r*K;
                        if(index<nbands &&index<static_cast<int>(phaY_prev.size()))
                        {
                         double wh=2.0*PI*r*f0[l-1]/fs;
                        phaY[index]=princ(phaY_prev[index]+wh*H);
                        }
                    }
                }
            }
            int delk=static_cast<int>(
                std::ceil(K/2.0));
            std::vector<std::complex<double>>
                windowFFT=fft(w[l]);

            for(int ii=-delk;
                ii<=delk;
                ii++)
            {
                if(ii==0)
                {
                    continue;
                }
                int absii=std::abs(ii);
                if(absii>=ML)
                {
                    continue;
                }
                double diff=
                    (ii>0 ? 1.0:-1.0)*
                    std::arg(windowFFT[absii]);
               for(int j=0;j<nharm;j++)
                {
                    int index=bh[j];
                    int target=index+ii;
                    if(target>=0 && target<nbands && index<nbands)
                    {
                    phaY[target]=princ(phaY[index]+diff);
                    }
                }
            }
            nharm_prev=nharm;
            phaY_prev=phaY;
            std::vector<double>phaY2(ML,0.0);
            for(int k=0;k<nbands;k++)
            {
             phaY2[k]=phaY[k];
           }
            std::vector<double>temp;
            for(int k=nbands-1;k>=0;k--)
            {
             temp.push_back(phaY[k]);
           }
            if(ML%2==0)
            {
                if(temp.size()>=2)
                {
                temp.erase(temp.begin());
                 temp.pop_back();
               }
            }
            else
            {
               if(!temp.empty())
                {
                temp.pop_back();
               }
            }
            int remaining=ML-nbands;
            for(int k=0;k<remaining &&k<static_cast<int>(temp.size());k++)
            {
             phaY2[nbands+k]=temp[k];
            }
            Y[l].resize(ML);
            for(int k=0;k<ML;k++)
            {
            Y[l][k]=magY[k]*std::complex<double>(std::cos(phaY2[k]),std::sin(phaY2[k]));
            }
        }
    }
    int lastMh=0;
    if(!Mh.empty())
    {
     lastMh=Mh.back();
    }
    int leny=std::max(static_cast<int>(x.size()),static_cast<int>(std::round(M_default/2.0+(nframes-1)*H+lastMh)));
    std::vector<double>y(leny,0.0);
    std::vector<double>NF(leny,0.0);

    double offy=M_default/2.0;
    for(int l=0;l<nframes;l++)
    {
        int start;
        if(M[l]%2==0)
        {
         start=static_cast<int>(std::round(offy))-Mh[l];
        }
        else
        {
         start=static_cast<int>(std::round(offy))-Mh[l]-1;
        }
        std::vector<std::complex<double>>timeSignal=ifft(Y[l]);
        for(int k=0;k<M[l];k++)
        {
            int pos=start+k;
            if(pos>=0 && pos<leny)
            {
                double value=timeSignal[k].real()*w[l][k];
                y[pos]+=value;
                NF[pos]+=w[l][k]*w[l][k];
            }
        }
        offy+=H;
    }
  for(int k=0;k<leny;k++)
    {
        if(NF[k]!=0.0)
        {
         y[k]/=NF[k];
        }
    }
    writeWavFile("enhanced_output.wav",y,fs);
    std::cout<<"\nProcessing completed.\n";
    std::cout<<"Output: enhanced_output.wav\n";
    return 0;
}