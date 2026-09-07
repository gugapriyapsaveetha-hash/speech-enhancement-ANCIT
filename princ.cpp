#include<vector>
#include<cmath>
double princ(double angle)
{
 double phi=std::fmod(angle,2.0*std::acos(-1.0));
 if(phi<0)
 {
  phi=phi+2.0*std::acos(-1.0);
 }
 if(phi<std::acos(-1.0))
 {
  return phi;
 }
 else
 {
  return phi-2.0*std::acos(-1.0);
 }
}
std::vector<double>princ(const std::vector<double>&angle)
{
 std::vector<double>result(angle.size());
 for(size_t i=0;i<angle.size();i++)
 {
  result[i]=princ(angle[i]);
 }
 return result;
}