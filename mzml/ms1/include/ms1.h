#ifndef MS1_H
#define MS1_H

#include <memory>

class Ms1
{
public:
    Ms1();
    Ms1(float mz , float intensity , float rt);
    Ms1(double mz , double intensity , double rt);
    Ms1(double mz , double intensity , float rt);
public:
    void SetMz(float mz);
    void SetIntensity(float intensity);
    void SetRt(float rt);
    float GetMz();
    float GetIntensity();
    float GetRt();
private:
    float m_mz = 0;
    float m_intensity = 0;
    float m_rt = 0;
};

#endif // MS1_H
