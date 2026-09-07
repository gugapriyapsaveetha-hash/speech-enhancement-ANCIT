#include<vector>
#include<string>
#include<fstream>
#include<cstdint>
bool readWavFile(const std::string&filename,std::vector<double>&data,int&fs)
{
 std::ifstream file(filename,std::ios::binary);
 if(!file)
 {
  return false;
 }

 char riff[4];
 char wave[4];
 file.read(riff,4);
 file.seekg(4,std::ios::cur);
 file.read(wave,4);
 if(riff[0]!='R'||riff[1]!='I'||riff[2]!='F'||riff[3]!='F')
 {
  return false;
 }
 if(wave[0]!='W'||wave[1]!='A'||wave[2]!='V'||wave[3]!='E')
 {
  return false;
 }
 uint16_t channels=0;
 uint16_t bits=0;
 uint32_t sampleRate=0;
 uint32_t dataSize=0;
 std::streampos dataPosition=0;
 while(file)
 {
  char id[4];
  uint32_t size=0;
  file.read(id,4);
  file.read(reinterpret_cast<char*>(&size),4);
  if(!file)
  {
   break;
  }
  if(id[0]=='f'&&id[1]=='m'&&id[2]=='t'&&id[3]==' ')
  {
   uint16_t format=0;
   uint32_t byteRate=0;
   uint16_t blockAlign=0;

   file.read(reinterpret_cast<char*>(&format),2);
   file.read(reinterpret_cast<char*>(&channels),2);
   file.read(reinterpret_cast<char*>(&sampleRate),4);
   file.read(reinterpret_cast<char*>(&byteRate),4);
   file.read(reinterpret_cast<char*>(&blockAlign),2);
   file.read(reinterpret_cast<char*>(&bits),2);

   if(size>16)
   {
    file.seekg(size-16,std::ios::cur);
   }
   if(format!=1)
   {
    return false;
   }
  }
  else if(id[0]=='d'&&id[1]=='a'&&id[2]=='t'&&id[3]=='a')
  {
   dataSize=size;
   dataPosition=file.tellg();
   file.seekg(size,std::ios::cur);
  }
  else
  {
   file.seekg(size,std::ios::cur);
  }
 }
 if(dataSize==0||channels==0||bits==0)
 {
  return false;
 }
 fs=static_cast<int>(sampleRate);
 file.clear();
 file.seekg(dataPosition);
 int bytesPerSample=bits/8;
 if(bytesPerSample<=0)
 {
  return false;
 }
 int totalSamples=static_cast<int>(dataSize/(bytesPerSample*channels));
 data.resize(totalSamples);
 for(int i=0;i<totalSamples;i++)
 {
  double value=0.0;
  if(bits==16)
  {
   int sum=0;
   for(int c=0;c<channels;c++)
   {
    int16_t sample=0;
    file.read(reinterpret_cast<char*>(&sample),2);
    if(!file)
    {
     return false;
    }
    sum+=sample;
   }
   value=static_cast<double>(sum)/(channels*32768.0);
  }
  else if(bits==8)
  {
   int sum=0;
   for(int c=0;c<channels;c++)
   {
    unsigned char sample=0;
    file.read(reinterpret_cast<char*>(&sample),1);
    if(!file)
    {
     return false;
    }
    sum+=sample;
   }
   value=(static_cast<double>(sum)/channels-128.0)/128.0;
  }
  else
  {
   return false;
  }
  data[i]=value;
 }
 file.close();
 return true;
}
void writeWavFile(const std::string&filename,const std::vector<double>&data,int fs)
{
 std::ofstream file(filename,std::ios::binary);
 if(!file)
 {
  return;
 }

 uint16_t channels=1;
 uint16_t bits=16;
 uint16_t format=1;
 uint32_t dataSize=static_cast<uint32_t>(data.size()*2);
 uint32_t fileSize=36+dataSize;
 uint32_t byteRate=static_cast<uint32_t>(fs*channels*bits/8);
 uint16_t blockAlign=static_cast<uint16_t>(channels*bits/8);
 uint32_t fmtSize=16;
 file.write("RIFF",4);
 file.write(reinterpret_cast<const char*>(&fileSize),4);
 file.write("WAVE",4);
 file.write("fmt ",4);
 file.write(reinterpret_cast<const char*>(&fmtSize),4);
 file.write(reinterpret_cast<const char*>(&format),2);
 file.write(reinterpret_cast<const char*>(&channels),2);
 uint32_t sampleRate=static_cast<uint32_t>(fs);
 file.write(reinterpret_cast<const char*>(&sampleRate),4);
 file.write(reinterpret_cast<const char*>(&byteRate),4);
 file.write(reinterpret_cast<const char*>(&blockAlign),2);
 file.write(reinterpret_cast<const char*>(&bits),2);
 file.write("data",4);
 file.write(reinterpret_cast<const char*>(&dataSize),4);
 for(int i=0;i<static_cast<int>(data.size());i++)
 {
  double value=data[i];
  if(value>1.0)
  {
   value=1.0;
  }
  if(value<-1.0)
  {
   value=-1.0;
  }
  int16_t sample=static_cast<int16_t>(value*32767.0);
  file.write(reinterpret_cast<const char*>(&sample),2);
 }
 file.close();
}