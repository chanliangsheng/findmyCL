#ifndef HEADGROUPFINDER_H
#define HEADGROUPFINDER_H

#include <workflow.h>
#include <MsLevelMatcher.h>
#include <cmath>

class HeadgroupFinder:public Workflow
{

public:
    static float m_ppm;//ppm，初始化为30
    static float m_ppm_with_half_score;//初始化为30
    static float m_mz_score_weight;//初始化为0.5
    static float m_k;//m_ppm_with_half_score对应的常数项，由m_ppm_with_half_score得出
public:
    HeadgroupFinder();
    HeadgroupFinder(MsLevelMatcher& ms_level_matcher);
    HeadgroupFinder(MsLevelMatcher& ms_level_matcher , float ppm);
    HeadgroupFinder(MsLevelMatcher& ms_level_matcher , float ppm , float ppm_with_half_score , float mz_score_weight);
public:
    void CheckMs2Headgroup();
    void CopyInfoFromCardiolipinMatchWithMs2(MsLevelMatcher& ms_level_matcher);
    void CopyInfoFromCardiolipinMatchWithMs2(MsLevelMatcher& ms_level_matcher , float ppm , float ppm_with_half_score , float mz_score_weight);
};

#endif // MS1MATCH_H
