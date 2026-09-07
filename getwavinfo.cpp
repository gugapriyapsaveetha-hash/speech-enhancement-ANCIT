#include<string>
#include<cstdint>
#include<fstream>
struct WavInfo
{
 int sampleRate;
 int numChannels;
 int bitsPerSample;
 int numSamples;
};
WavInfo getWavInfo(const std::string&filename)
{
 WavInfo info={0,0,0,0};
 std::ifstream file(filename,std::ios::binary);
 if(!file)
 {
  return info;
 }
 char riff[4];
 char wave[4];
 file.read(riff,4);
 file.seekg(4,std::ios::cur);
 file.read(wave,4);
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
   uint16_t audioFormat=0;
   uint16_t channels=0;
   uint32_t sampleRate=0;
   uint32_t byteRate=0;
   uint16_t blockAlign=0;
   uint16_t bits=0;
   file.read(reinterpret_cast<char*>(&audioFormat),2);
   file.read(reinterpret_cast<char*>(&channels),2);
   file.read(reinterpret_cast<char*>(&sampleRate),4);
   file.read(reinterpret_cast<char*>(&byteRate),4);
   file.read(reinterpret_cast<char*>(&blockAlign),2);
   file.read(reinterpret_cast<char*>(&bits),2);
   info.sampleRate=sampleRate;
   info.numChannels=channels;
   info.bitsPerSample=bits;
   if(size>16)
   {
    file.seekg(size-16,std::ios::cur);
   }
  }
  else if(id[0]=='d'&&id[1]=='a'&&id[2]=='t'&&id[3]=='a')
  {
   int bytesPerSample=info.bitsPerSample/8;
   if(bytesPerSample>0&&info.numChannels>0)
   {
    info.numSamples=size/(bytesPerSample*info.numChannels);
   }
   break;
  }
  else
  {
   file.seekg(size,std::ios::cur);
  }
 }
 file.close();
 return info;
}