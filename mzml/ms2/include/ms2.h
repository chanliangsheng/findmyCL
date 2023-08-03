#ifndef MS2_H
#define MS2_H


#include <vector>
#include <algorithm>
#include <QDebug>
#include <queue>
#include <cmath>
#include <headgroup.h>
#include <pa.h>
class Ms2
{
public:
    Ms2();
    Ms2(float precursor_ion_mz, float precursor_ion_intensity , float rt , std::vector<float> &fragment_ion_mz , std::vector<float> &fragment_ion_intensity);
    Ms2(float precursor_ion_mz, float precursor_ion_intensity , float rt , std::vector<double> &fragment_ion_mz , std::vector<double> &fragment_ion_intensity);
public:
    void SetPrecuisorIonMz(float precuisor_ion_mz);
    void SetPrecuisorIonIntensity(float precursor_ion_intensity);
    void SetRt(float rt);
    void SetFragmentIonMz(std::vector<float> fragment_ion_mz);
    void SetFragmentIonIntensity(std::vector<float> fragment_ion_intensity);
    void SetHeadgroup(std::vector<Headgroup> headgroup);

    float GetPrecuisorIonMz();
    float GetPrecuisorIonIntensity();
    float GetRt();
    std::vector<float> GetFragmentIonMz();
    std::vector<float> GetFragmentIonIntensity();
    float GetTotalIntensity();
    void CalculateTotalIntensity();
    std::vector<Headgroup> GetHeadgroup();

    void ClearHeadgroup();
    void ClearPaInfo();
    void ClearFaInfo();
    float GetMaxFragmentIntensity();
    float GetMinFragmentIntensity();
    void EmplaceBackHeadgroup(Headgroup& headgroup);
    void EmplaceBackPa(Pa pa);
    void EmplaceBackFa(Fa fa);
    int GetPaCount();
    int GetFaCount();
    std::vector<Fa>* GetFaVectorPtr();
    std::vector<Pa>* GetPaVectorPtr();
    void SortPaByChainLength();
    bool m_pa_sort_by_chain_length = 0;
    void SortFaByChainLength();
    bool m_fa_sort_by_chain_length = 0;
public:
    void DeleteLowIntensityFragment(float radio);//删除intensity最低的radio那么多个
private:
    //前体离子的信息
    float m_precursor_ion_mz = 0;
    float m_precursor_ion_intensity = 0;
    float m_total_intensity = 0;
    float m_rt = 0;
    //碎片离子的信息
    std::vector<float> m_fragment_ion_mz;
    std::vector<float> m_fragment_ion_intensity;
    //头基信息
    std::vector<Headgroup> m_headgroup;
    //Fa和Pa的信息
    std::vector<Pa> m_pa;
    std::vector<Fa> m_fa;
};


#endif // MS2_H
