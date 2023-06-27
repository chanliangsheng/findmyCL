#include <HeadgroupFinder.h>


float HeadgroupFinder::m_ppm = 30;
float HeadgroupFinder::m_ppm_with_half_score = 30;
float HeadgroupFinder::m_mz_score_weight = 0.5;
float HeadgroupFinder::m_k = -((1.17741/pow((HeadgroupFinder::m_ppm_with_half_score)/static_cast<double>(1000000), 2))/2);

using namespace std;
HeadgroupFinder::HeadgroupFinder()
{

}

HeadgroupFinder::HeadgroupFinder(MsLevelMatcher& ms_level_matcher)
{
    //获得一级找二级结果的复制
    this->m_cl_vector = ms_level_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms_level_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms_level_matcher.GetCopyDlclPairVector();
}

HeadgroupFinder::HeadgroupFinder(MsLevelMatcher& ms_level_matcher, float ppm)
{
    //获得一级找二级结果的复制
    this->m_cl_vector = ms_level_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms_level_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms_level_matcher.GetCopyDlclPairVector();
    this->m_ppm = ppm;
}

HeadgroupFinder::HeadgroupFinder(MsLevelMatcher& ms_level_matcher, float ppm, float ppm_with_half_score, float mz_score_weight)
{
    this->m_cl_vector = ms_level_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms_level_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms_level_matcher.GetCopyDlclPairVector();
    this->m_ppm = ppm;
    this->m_ppm_with_half_score = ppm_with_half_score;
    this->m_mz_score_weight = mz_score_weight;
    this->m_k = -((1.17741/pow((this->m_ppm_with_half_score)/static_cast<double>(1000000), 2))/2);

}

void HeadgroupFinder::CheckMs2Headgroup()
{
    //更新m_k
    this->m_k = -((1.17741/pow((this->m_ppm_with_half_score)/static_cast<double>(1000000), 2))/2);
    //cl
    for(auto itr = this->m_cl_vector.begin() ; itr != this->m_cl_vector.end();){
        itr->first.CheckHeadgroup(this->m_ppm , this->m_mz_score_weight , this->m_k);
        itr->second.CheckHeadgroup(this->m_ppm , this->m_mz_score_weight , this->m_k);
        //如果两个都没有找到头基，或者一个原本是没有的，另一个没有找到头基，则将其删除
        if(itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            itr = this->m_cl_vector.erase(itr);
        }
        else{
            ++itr;
        }
    }
    //mlcl
    for(auto itr = this->m_mlcl_vector.begin() ; itr != this->m_mlcl_vector.end();){
        itr->first.CheckHeadgroup(this->m_ppm , this->m_mz_score_weight , this->m_k);
        itr->second.CheckHeadgroup(this->m_ppm , this->m_mz_score_weight , this->m_k);
        //如果两个都没有找到头基，或者一个原本是没有的，另一个没有找到头基，则将其删除
        if(itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            itr = this->m_mlcl_vector.erase(itr);
        }
        else{
            ++itr;
        }
    }
    //dlcl
    for(auto itr = this->m_dlcl_vector.begin() ; itr != this->m_dlcl_vector.end();){
        itr->first.CheckHeadgroup(this->m_ppm , this->m_mz_score_weight , this->m_k);
        itr->second.CheckHeadgroup(this->m_ppm , this->m_mz_score_weight , this->m_k);
        //如果两个都没有找到头基，或者一个原本是没有的，另一个没有找到头基，则将其删除
        if(itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            itr = this->m_dlcl_vector.erase(itr);
        }
        else{
            ++itr;
        }
    }

    this->DeleteRedundantPair();//删除冗余心磷脂
}

void HeadgroupFinder::CopyInfoFromCardiolipinMatchWithMs2(MsLevelMatcher& ms_level_matcher)
{
    //获得一级找二级结果的复制
    this->m_cl_vector = ms_level_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms_level_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms_level_matcher.GetCopyDlclPairVector();
}

void HeadgroupFinder::CopyInfoFromCardiolipinMatchWithMs2(MsLevelMatcher& ms_level_matcher, float ppm, float ppm_with_half_score, float mz_score_weight)
{
    this->m_cl_vector = ms_level_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms_level_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms_level_matcher.GetCopyDlclPairVector();
    this->m_ppm = ppm;
    this->m_ppm_with_half_score = ppm_with_half_score;
    this->m_mz_score_weight = mz_score_weight;
    this->m_k = -((1.17741/pow((this->m_ppm_with_half_score)/static_cast<double>(1000000), 2))/2);
}
