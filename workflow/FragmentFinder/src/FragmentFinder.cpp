#include <FragmentFinder.h>

float FragmentFinder::m_ppm = 30;
float FragmentFinder::m_ppm_with_half_score = 30;
float FragmentFinder::m_mz_score_weight = 0.5;
float FragmentFinder::m_k = -((1.17741/pow((FragmentFinder::m_ppm_with_half_score)/static_cast<double>(1000000), 2))/2);

FragmentFinder::FragmentFinder()
{

}

FragmentFinder::FragmentFinder(HeadgroupFinder& headgroup_finder)
{
    this->m_cl_vector = headgroup_finder.GetCopyClPairVector();
    this->m_mlcl_vector = headgroup_finder.GetCopyMlclPairVector();
    this->m_dlcl_vector = headgroup_finder.GetCopyDlclPairVector();
}

FragmentFinder::FragmentFinder(HeadgroupFinder& headgroup_finder, float ppm, float ppm_with_half_score, float mz_score_weight)
{
    this->m_cl_vector = headgroup_finder.GetCopyClPairVector();
    this->m_mlcl_vector = headgroup_finder.GetCopyMlclPairVector();
    this->m_dlcl_vector = headgroup_finder.GetCopyDlclPairVector();
    this->m_ppm = ppm;
    this->m_ppm_with_half_score = ppm_with_half_score;
    this->m_mz_score_weight = mz_score_weight;
    this->m_k = -((1.17741/pow((ppm_with_half_score)/static_cast<double>(1000000), 2))/2);;
}

void FragmentFinder::CopyInfoFromMs2WithHeadgroup(HeadgroupFinder& headgroup_finder)
{
    this->m_cl_vector = headgroup_finder.GetCopyClPairVector();
    this->m_mlcl_vector = headgroup_finder.GetCopyMlclPairVector();
    this->m_dlcl_vector = headgroup_finder.GetCopyDlclPairVector();
}

void FragmentFinder::FindMs2PaAndFa(Database& database)
{
    //更新m_k
    this->m_k = -((1.17741/pow((this->m_ppm_with_half_score)/static_cast<double>(1000000), 2))/2);
    //cl
    for(auto itr = this->m_cl_vector.begin() ; itr != this->m_cl_vector.end();){
        itr->first.FindPaAndFa(this->m_ppm , this->m_mz_score_weight , this->m_k , database);
        itr->second.FindPaAndFa(this->m_ppm , this->m_mz_score_weight , this->m_k , database);

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
        itr->first.FindPaAndFa(this->m_ppm , this->m_mz_score_weight , this->m_k , database);
        itr->second.FindPaAndFa(this->m_ppm , this->m_mz_score_weight , this->m_k , database);
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
        itr->first.FindPaAndFa(this->m_ppm , this->m_mz_score_weight , this->m_k , database);
        itr->second.FindPaAndFa(this->m_ppm , this->m_mz_score_weight , this->m_k , database);
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


