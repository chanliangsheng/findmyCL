#ifndef FRAGMENTFINDER_H
#define FRAGMENTFINDER_H

#include <cmath>
#include <HeadgroupFinder.h>

class FragmentFinder:public Workflow
{
public:
    FragmentFinder();
    FragmentFinder(HeadgroupFinder& headgroup_finder);
    FragmentFinder(HeadgroupFinder& headgroup_finder , float ppm , float ppm_with_half_score , float mz_score_weight);
public:
    void CopyInfoFromMs2WithHeadgroup(HeadgroupFinder& headgroup_finder);
public:
    static float m_ppm;//ppm，初始化为30
    static float m_ppm_with_half_score;//初始化为5
    static float m_mz_score_weight;//初始化为0.5
    static float m_k;//m_ppm_with_half_score对应的常数项，由m_ppm_with_half_score得出
public:
    void FindMs2PaAndFa(Database& database);
};

#endif // MS1MATCH_H
