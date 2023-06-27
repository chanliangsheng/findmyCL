#include "cardiolipin.h"


using namespace std;
Cardiolipin::Cardiolipin()
{
    this->m_Ms1_ptr = nullptr;//初始化为空指针
    this->m_Ms1_DatabaseRecord_ptr = nullptr;//初始化为空指针
    this->m_Ms2_vector_ptr = vector<Ms2*>();//初始化为空
}


Cardiolipin::Cardiolipin(Ms1 *Ms1_ptr, DatabaseRecord *DatabaseRecord_ptr)
{
    this->m_Ms1_ptr = Ms1_ptr;
    this->m_Ms1_DatabaseRecord_ptr = DatabaseRecord_ptr;
    this->m_Ms2_vector_ptr = vector<Ms2*>();//初始化为空
}

void Cardiolipin::SetMs1Ptr(Ms1 *ms1_ptr)
{
    this->m_Ms1_ptr = ms1_ptr;
}

void Cardiolipin::SetMs1DatabaseRecordPtr(DatabaseRecord *ms1_databaseRecord_ptr)
{
    this->m_Ms1_DatabaseRecord_ptr = ms1_databaseRecord_ptr;
}

void Cardiolipin::SetMs1MatchingScosre(float ms1_matching_score)
{
    this->m_Ms1_Matching_Score = ms1_matching_score;
}

void Cardiolipin::SetMs2VectorPtr(std::vector<Ms2 *> ms2_vector_ptr)
{
    this->m_Ms2_vector_ptr = ms2_vector_ptr;
}

Ms1 *Cardiolipin::GetMs1Ptr()
{
    return this->m_Ms1_ptr;
}

DatabaseRecord *Cardiolipin::GetMs1DatabaseRecordPtr()
{
    return this->m_Ms1_DatabaseRecord_ptr;
}

float Cardiolipin::GetMs1MatchingScore()
{
    return this->m_Ms1_Matching_Score;
}

std::vector<Ms2 *> Cardiolipin::GetMs2VectorPtr()
{
    return this->m_Ms2_vector_ptr;
}

unsigned int Cardiolipin::GetChainLength()
{
    return this->m_Ms1_DatabaseRecord_ptr->GetChainLength();
}

unsigned int Cardiolipin::GetUnsaturation()
{
    return this->m_Ms1_DatabaseRecord_ptr->GetUnsaturation();
}

unsigned int Cardiolipin::GetOxygen()
{
    return this->m_Ms1_DatabaseRecord_ptr->GetOxygen();
}

void Cardiolipin::ClearMs2Info()
{
    vector<Ms2*>().swap(this->m_Ms2_vector_ptr);//清空vector中的元素
}

QString Cardiolipin::ShowInfo()
{
    //如果是空对象，则返回空字符串
    if(this->CheckEmptyObject()){
        return "";
    }

    QString message;
    message = message + "MS1 informations->" + "MS1-mz:" + QString::number(this->m_Ms1_ptr->GetMz()) + ";" + "MS1-intensity:" + QString::number(this->m_Ms1_ptr->GetIntensity()) + ";" + "MS1-rt:" + QString::number(this->m_Ms1_ptr->GetRt() / 60) + " min" + "\n";
    message = message + "MS1 matching->" + "chain length:" + QString::number(this->GetChainLength()) + ";" + "unsaturation:" + QString::number(this->GetUnsaturation()) + ";" + "oxygen:" + QString::number(this->GetOxygen()) + "\n";

    //输出M-H二级信息
    message = message + "MS2:" + "\n";
    for(auto ms2_itr = this->m_Ms2_vector_ptr.begin() ; ms2_itr != this->m_Ms2_vector_ptr.end() ; ms2_itr++){
        float max_intensity = (*ms2_itr)->GetMaxFragmentIntensity();//最大的intensity
        float min_intensity = (*ms2_itr)->GetMinFragmentIntensity();//最小的intensity
        float precuisor_ion_mz = (*ms2_itr)->GetPrecuisorIonMz();
        message = message + "MS2 informations->" + "precuisor ion mz:" + QString::number(precuisor_ion_mz) + "max intensity:" + QString::number(max_intensity) + ";" + "min intensity" + QString::number(min_intensity) + ";" + "rt:" + QString::number((*ms2_itr)->GetRt() / 60) + " min"+ "\n";
    }
    return message;
}

std::array<unsigned int, 3> Cardiolipin::GetMs1DatabaseRecordComponent()
{
    //如果是空对象，则返回3个0
    if(this->CheckEmptyObject()){
        return array<unsigned int , 3>({0,0,0});
    }
    //否则，返回真实值
    else{
        return array<unsigned int , 3>({this->GetChainLength() , this->GetUnsaturation() , this->GetOxygen()});;
    }
}

QString Cardiolipin::GetMs1DatabaseRecordComponentWithQString()
{
    //如果是空对象，则返回3个0
    if(this->CheckEmptyObject()){
        return QString("0:0:0");
    }
    //否则，返回真实值
    else{
        QString result = QString::number(this->GetChainLength()) + ":" + QString::number(this->GetUnsaturation()) +  ":" + QString::number(this->GetOxygen());
        return result;
    }
}

float Cardiolipin::GetSampleMz()
{
    //如果是空对象，则返回0
    if(this->CheckEmptyObject()){
        return 0;
    }
    //否则，返回真实值
    else{
        return this->m_Ms1_ptr->GetMz();
    }
}

float Cardiolipin::GetSampleRt()
{
    //如果是空对象，则返回0
    if(this->CheckEmptyObject()){
        return 0;
    }
    //否则，返回真实值
    else{
        return this->m_Ms1_ptr->GetRt();
    }
}

float Cardiolipin::GetSampleIntensity()
{
    //如果是空对象，则返回0
    if(this->CheckEmptyObject()){
        return 0;
    }
    //否则，返回真实值
    else{
        return this->m_Ms1_ptr->GetIntensity();
    }
}

void Cardiolipin::ScoreMs1(float &k)
{
    if((this->m_Ms1_ptr != nullptr) && (this->m_Ms1_DatabaseRecord_ptr != nullptr)){
        float real_ppm = (this->m_Ms1_ptr->GetMz() - this->m_Ms1_DatabaseRecord_ptr->GetMz()) / (this->m_Ms1_DatabaseRecord_ptr->GetMz());//计算真实的ppm
        this->m_Ms1_Matching_Score = exp(k * pow(real_ppm , 2));//返回分数
    }
}

void Cardiolipin::MatchCardiolipinWithMs2(std::vector<Ms2> &ms2_vector, float &dalton, float &tolerance_rt_scope)
{
    //如果这个Cardiolipin是空的Cardiolipin，则直接返回
    if(this->CheckEmptyObject()){
        return;
    }

    // vector<Ms2>需要已经根据precursor_ion_mz排序
    int left = 0;
    int right = ms2_vector.size() - 1;
    while(left <= right){
        int mid = (left + right) / 2;
        float min_mz = ms2_vector[mid].GetPrecuisorIonMz() - dalton;//最小mz匹配值
        float max_mz = ms2_vector[mid].GetPrecuisorIonMz()  + dalton;//最大mz匹配值
        if((this->m_Ms1_ptr->GetMz() >= min_mz) && (this->m_Ms1_ptr->GetMz() <= max_mz)){
            float min_rt = ms2_vector[mid].GetRt() - tolerance_rt_scope;//最小rt匹配值
            float max_rt = ms2_vector[mid].GetRt() + tolerance_rt_scope;//最大rt匹配值
            if((this->m_Ms1_ptr->GetRt() >= min_rt) && (this->m_Ms1_ptr->GetRt() <= max_rt)){
                //如果rt也符合要求，则把MS2的地址放入m_Ms2_vector中
                this->m_Ms2_vector_ptr.emplace_back(&ms2_vector[mid]);
            }
            //mid的左边一位和右边一位的各自的符合区间
            int left_t = mid - 1;
            int right_t = mid + 1;

            float left_t_min_mz;
            float left_t_max_mz;
            float right_t_min_mz;
            float right_t_max_mz;
            while(left_t >= 0){
                left_t_min_mz = ms2_vector[left_t].GetPrecuisorIonMz() - dalton;
                left_t_max_mz = ms2_vector[left_t].GetPrecuisorIonMz() + dalton;
                float left_t_min_rt = ms2_vector[left_t].GetRt() - tolerance_rt_scope;//最小rt匹配值
                float left_t_max_rt = ms2_vector[left_t].GetRt() + tolerance_rt_scope;//最大rt匹配值
                if((this->m_Ms1_ptr->GetMz() >= left_t_min_mz) && (this->m_Ms1_ptr->GetMz() <= left_t_max_mz)){
                    //如果mz是符合的，那么判断rt是否符合；如果mz符合，rt不符合，继续向左寻找，直到mz不符合；如果mz符合，rt符合，则加入ms2_vector，继续向左寻找
                    if((this->m_Ms1_ptr->GetRt() >= left_t_min_rt) && (this->m_Ms1_ptr->GetRt() <= left_t_max_rt)){
                        this->m_Ms2_vector_ptr.emplace_back(&ms2_vector[left_t]);
                    }
                    left_t--;//向左继续找
                }
                else{
                    break;
                }
            }
            while(right_t <= ms2_vector.size() - 1){
                right_t_min_mz = ms2_vector[right_t].GetPrecuisorIonMz() - dalton;
                right_t_max_mz = ms2_vector[right_t].GetPrecuisorIonMz() + dalton;
                float right_t_min_rt = ms2_vector[right_t].GetRt() - tolerance_rt_scope;//最小rt匹配值
                float right_t_max_rt = ms2_vector[right_t].GetRt() + tolerance_rt_scope;//最大rt匹配值
                if((this->m_Ms1_ptr->GetMz() >= right_t_min_mz) && (this->m_Ms1_ptr->GetMz() <= right_t_max_mz)){
                    //如果mz是符合的，那么判断rt是否符合；如果mz符合，rt不符合，继续向右寻找，直到mz不符合；如果mz符合，rt符合，则加入ms2_vector，继续向右寻找
                    if((this->m_Ms1_ptr->GetRt() >= right_t_min_rt) && (this->m_Ms1_ptr->GetRt() <= right_t_max_rt)){
                        this->m_Ms2_vector_ptr.emplace_back(&ms2_vector[right_t]);
                    }
                    right_t++;//向右继续寻找
                }
                else{
                    break;
                }
            }
            break;
        }
        else if(this->m_Ms1_ptr->GetMz() > min_mz){
            left = mid + 1;
        }
        else if(this->m_Ms1_ptr->GetMz() < max_mz){
            right = mid - 1;
        }
    }
}

bool Cardiolipin::CheckMs2Exist()
{
    if(this->m_Ms2_vector_ptr.size() == 0){
        return 0;
    }
    else{
        return 1;
    }
}

bool Cardiolipin::CheckEmptyObject()
{
    if((this->m_Ms1_ptr == nullptr) && (this->m_Ms1_DatabaseRecord_ptr == nullptr)){
        return 1;
    }
    else{
        return 0;
    }
}

void Cardiolipin::EmptyObject()
{
    if(this->m_Ms1_ptr != nullptr){
        this->m_Ms1_ptr = nullptr;
    }
    if(this->m_Ms1_DatabaseRecord_ptr != nullptr){
        this->m_Ms1_DatabaseRecord_ptr = nullptr;
    }

    this->m_Ms1_Matching_Score = 0;//一级配对分数设置为0
    vector<Ms2*>().swap(this->m_Ms2_vector_ptr);//清空vector中的元素
}

void Cardiolipin::CheckHeadgroup(float& ppm , float& mz_score_weight , float& k)
{
    //如果是空对象，则直接返回
    if(this->CheckEmptyObject()){
        return;
    }

    //配对的最小值和最大值
    float min_mz = 152.9963485799 - 152.9963485799 * ppm/1000000;
    float max_mz = 152.9963485799 + 152.9963485799 * ppm/1000000;

    for(auto ms2_itr = this->m_Ms2_vector_ptr.begin() ; ms2_itr != this->m_Ms2_vector_ptr.end();){
        (*ms2_itr)->ClearHeadgroup();//清空头基信息
        bool find_head_group = 0;//判断是否找到了头基
        float max_intensity = (*ms2_itr)->GetMaxFragmentIntensity();//最大的intensity

        vector<float> fragment_ion_mz = (*ms2_itr)->GetFragmentIonMz();
        vector<float> fragment_ion_intensity = (*ms2_itr)->GetFragmentIonIntensity();
        for(auto fragment_mz_itr = fragment_ion_mz.begin() ; fragment_mz_itr != fragment_ion_mz.end() ; fragment_mz_itr++){
            //如果找到了头基
            if((*fragment_mz_itr >= min_mz) && (*fragment_mz_itr <= max_mz)){
                find_head_group = 1;//设置为真
                float real_ppm = (*fragment_mz_itr - 152.9963485799) / (152.9963485799);//计算真实的ppm
                float mz_score = exp(k * pow(real_ppm , 2));//mz分数
                int intensity_pos = distance(fragment_ion_mz.begin() , fragment_mz_itr);
                float intensity = fragment_ion_intensity[intensity_pos];
                float intensity_score = intensity / max_intensity;//intensity分数
                float total_score = mz_score_weight * mz_score + (1 - mz_score_weight)*intensity_score;//最终分数
                Headgroup new_headgroup = Headgroup(*fragment_mz_itr , intensity , total_score);//新建头基类
                (*ms2_itr)->EmplaceBackHeadgroup(new_headgroup);
            }
        }
        //如果这个二级中找到头基
        if(find_head_group){
            ++ms2_itr;
        }
        //如果没有找到头基，则将这个二级删除
        else{
            ms2_itr = this->m_Ms2_vector_ptr.erase(ms2_itr);
        }
    }

    //如果这个心磷脂的所有二级都没有找到头基，则将心磷脂的信息全部置为空
    if(this->m_Ms2_vector_ptr.size() == 0){
        this->EmptyObject();
    }
}

void Cardiolipin::FindPaAndFa(float &ppm, float &mz_score_weight, float &k, Database &database)
{
    //如果是空对象，则直接返回
    if(this->CheckEmptyObject()){
        return;
    }

    //取出PA和FA数据库，都已经根据mz进行排序
    vector<DatabaseRecord>* pa_match = database.GetLocalPA();
    vector<DatabaseRecord>* fa_match = database.GetLocalFA();
    //PA
    for(auto ms2_itr = this->m_Ms2_vector_ptr.begin() ; ms2_itr != this->m_Ms2_vector_ptr.end();ms2_itr++){
        (*ms2_itr)->ClearPaInfo();//把上次PA的搜索信息清空
        float max_intensity = (*ms2_itr)->GetMaxFragmentIntensity();//最大的intensity

        vector<float> fragment_ion_mz = (*ms2_itr)->GetFragmentIonMz();
        vector<float> fragment_ion_intensity = (*ms2_itr)->GetFragmentIonIntensity();
        for(auto fragment_mz_itr = fragment_ion_mz.begin() ; fragment_mz_itr != fragment_ion_mz.end() ; fragment_mz_itr ++){
            int left = 0;
            int right = pa_match->size() - 1;
            while(left <= right){
                int mid = (left + right)/2;
                //基于ppm所给定的范围，需要样本中的mz位于min_mz和max_mz中，才认为是符合的
                float min_mz = pa_match->at(mid).GetMz() - pa_match->at(mid).GetMz() * ppm/1000000;
                float max_mz = pa_match->at(mid).GetMz() + pa_match->at(mid).GetMz() * ppm/1000000;
                if((*fragment_mz_itr >= min_mz) && (*fragment_mz_itr <= max_mz)){
                    unsigned int oxygen_limit = ceil((pa_match->at(mid).GetChainLength() - (pa_match->at(mid).GetUnsaturation() * 2 - 2))/2);//允许的氧的最大个数
                    int intensity_pos = distance(fragment_ion_mz.begin() , fragment_mz_itr);//某个intensity对应的位置
                    float intensity = fragment_ion_intensity[intensity_pos];//取出intensity
                    float intensity_score = intensity / max_intensity;//intensity分数

                    //如果这个PA的氧个数超过了允许的最大个数，则说明不是PA，则退出
                    if(oxygen_limit < pa_match->at(mid).GetOxygen()){

                    }
                    else{
                        float real_ppm = (*fragment_mz_itr - pa_match->at(mid).GetMz())/pa_match->at(mid).GetMz();
                        float mz_score =  exp(k * pow(real_ppm , 2));//mz分数
                        float total_score = mz_score_weight * mz_score + (1 - mz_score_weight)*intensity_score;//最终分数
                        (*ms2_itr)->EmplaceBackPa(Pa(*fragment_mz_itr , intensity , total_score , &pa_match->at(mid)));//结果加入二级中
                    }
                    //mid的左边一位和右边一位的各自的符合区间
                    int left_t = mid - 1;
                    int right_t = mid + 1;
                    float left_t_min_mz;
                    float left_t_max_mz;
                    float right_t_min_mz;
                    float right_t_max_mz;
                    while(left_t >= 0){
                        //更新区间
                        left_t_min_mz = pa_match->at(left_t).GetMz() - pa_match->at(left_t).GetMz() * ppm/1000000;
                        left_t_max_mz = pa_match->at(left_t).GetMz() + pa_match->at(left_t).GetMz() * ppm/1000000;
                        if((*fragment_mz_itr >= left_t_min_mz) && (*fragment_mz_itr <= left_t_max_mz)){
                            unsigned int oxygen_limit = ceil((pa_match->at(left_t).GetChainLength() - (pa_match->at(left_t).GetUnsaturation() * 2 - 2))/2);//允许的氧的最大个数
                            if(oxygen_limit >= pa_match->at(left_t).GetOxygen()){
                                float real_ppm = (*fragment_mz_itr - pa_match->at(left_t).GetMz())/pa_match->at(left_t).GetMz();
                                float mz_score =  exp(k * pow(real_ppm , 2));//mz分数
                                float total_score = mz_score_weight * mz_score + (1 - mz_score_weight)*intensity_score;//最终分数，这里的intensity可以用上面算过的
                                (*ms2_itr)->EmplaceBackPa(Pa(*fragment_mz_itr , intensity , total_score , &pa_match->at(left_t)));//结果加入二级中
                            }
                            left_t--;//向左继续找
                        }
                        else{
                            break;
                        }
                    }
                    while(right_t <= pa_match->size() - 1){
                        //更新区间
                        right_t_min_mz = pa_match->at(right_t).GetMz() - pa_match->at(right_t).GetMz() * ppm/1000000;
                        right_t_max_mz = pa_match->at(right_t).GetMz() + pa_match->at(right_t).GetMz() * ppm/1000000;
                        if((*fragment_mz_itr >= right_t_min_mz) && (*fragment_mz_itr <= right_t_max_mz)){
                            unsigned int oxygen_limit = ceil((pa_match->at(right_t).GetChainLength() - (pa_match->at(right_t).GetUnsaturation() * 2 - 2))/2);//允许的氧的最大个数
                            if(oxygen_limit >= pa_match->at(right_t).GetOxygen()){
                                float real_ppm = (*fragment_mz_itr - pa_match->at(right_t).GetMz())/pa_match->at(right_t).GetMz();
                                float mz_score =  exp(k * pow(real_ppm , 2));//mz分数
                                float total_score = mz_score_weight * mz_score + (1 - mz_score_weight)*intensity_score;//最终分数，这里的intensity可以用上面算过的
                                (*ms2_itr)->EmplaceBackPa(Pa(*fragment_mz_itr , intensity , total_score , &pa_match->at(right_t)));//结果加入二级中
                            }
                            right_t++;//向右继续寻找
                        }
                        else{
                            break;
                        }
                    }
                    break;
                }
                else if(*fragment_mz_itr >= max_mz){
                    left = mid + 1;
                }
                else if(*fragment_mz_itr <= min_mz){
                    right = mid - 1;
                }
            }
        }
    }

    //Fa
    for(auto ms2_itr = this->m_Ms2_vector_ptr.begin() ; ms2_itr != this->m_Ms2_vector_ptr.end();ms2_itr++){
        (*ms2_itr)->ClearFaInfo();//把上次的FA的搜索信息清空
        float max_intensity = (*ms2_itr)->GetMaxFragmentIntensity();//最大的intensity

        vector<float> fragment_ion_mz = (*ms2_itr)->GetFragmentIonMz();
        vector<float> fragment_ion_intensity = (*ms2_itr)->GetFragmentIonIntensity();
        for(auto fragment_mz_itr = fragment_ion_mz.begin() ; fragment_mz_itr != fragment_ion_mz.end() ; fragment_mz_itr ++){
            int left = 0;
            int right = fa_match->size() - 1;
            while(left <= right){
                int mid = (left + right)/2;
                //基于ppm所给定的范围，需要样本中的mz位于min_mz和max_mz中，才认为是符合的
                float min_mz = fa_match->at(mid).GetMz() - fa_match->at(mid).GetMz()*ppm/1000000;
                float max_mz = fa_match->at(mid).GetMz() + fa_match->at(mid).GetMz()*ppm/1000000;

                if((*fragment_mz_itr >= min_mz) && (*fragment_mz_itr <= max_mz)){
                    unsigned int oxygen_limit = ceil((fa_match->at(mid).GetChainLength() - (fa_match->at(mid).GetUnsaturation() * 2 - 2))/2);//允许的氧的最大个数
                    int intensity_pos = distance(fragment_ion_mz.begin() , fragment_mz_itr);//某个intensity对应的位置
                    float intensity = fragment_ion_intensity[intensity_pos];//取出intensity
                    float intensity_score = intensity / max_intensity;//intensity分数

                    //如果这个PA的氧个数超过了允许的最大个数，则说明不是PA，则退出
                    if(oxygen_limit < fa_match->at(mid).GetOxygen()){

                    }
                    else{
                        float real_ppm = (*fragment_mz_itr - fa_match->at(mid).GetMz())/fa_match->at(mid).GetMz();
                        float mz_score =  exp(k * pow(real_ppm , 2));//mz分数
                        float total_score = mz_score_weight * mz_score + (1 - mz_score_weight)*intensity_score;//最终分数
                        (*ms2_itr)->EmplaceBackFa(Fa(*fragment_mz_itr , intensity , total_score , &fa_match->at(mid)));
                    }
                    //mid的左边一位和右边一位的各自的符合区间
                    int left_t = mid - 1;
                    int right_t = mid + 1;
                    float left_t_min_mz;
                    float left_t_max_mz;
                    float right_t_min_mz;
                    float right_t_max_mz;
                    while(left_t >= 0){
                        //更新区间
                        left_t_min_mz = fa_match->at(left_t).GetMz() - fa_match->at(left_t).GetMz()*ppm/1000000;
                        left_t_max_mz = fa_match->at(left_t).GetMz() + fa_match->at(left_t).GetMz()*ppm/1000000;
                        if((*fragment_mz_itr >= left_t_min_mz) && (*fragment_mz_itr <= left_t_max_mz)){
                            unsigned int oxygen_limit = ceil((fa_match->at(left_t).GetChainLength() - (fa_match->at(left_t).GetUnsaturation() * 2 - 2))/2);//允许的氧的最大个数
                            if(oxygen_limit >= fa_match->at(left_t).GetOxygen()){
                                float real_ppm = (*fragment_mz_itr - fa_match->at(left_t).GetMz())/fa_match->at(left_t).GetMz();
                                float mz_score =  exp(k * pow(real_ppm , 2));//mz分数
                                float total_score = mz_score_weight * mz_score + (1 - mz_score_weight)*intensity_score;//最终分数，这里的intensity可以用上面算过的
                                (*ms2_itr)->EmplaceBackFa(Fa(*fragment_mz_itr , intensity , total_score , &fa_match->at(left_t)));
                            }
                            left_t--;//向左继续找
                        }
                        else{
                            break;
                        }
                    }
                    while(right_t <= fa_match->size() - 1){
                        //更新区间
                        right_t_min_mz = fa_match->at(right_t).GetMz() - fa_match->at(right_t).GetMz()*ppm/1000000;
                        right_t_max_mz = fa_match->at(right_t).GetMz() + fa_match->at(right_t).GetMz()*ppm/1000000;
                        if((*fragment_mz_itr >= right_t_min_mz) && (*fragment_mz_itr <= right_t_max_mz)){
                            unsigned int oxygen_limit = ceil((fa_match->at(right_t).GetChainLength() - (fa_match->at(right_t).GetUnsaturation() * 2 - 2))/2);//允许的氧的最大个数
                            if(oxygen_limit >= fa_match->at(right_t).GetOxygen()){
                                float real_ppm = (*fragment_mz_itr - fa_match->at(right_t).GetMz())/fa_match->at(right_t).GetMz();
                                float mz_score =  exp(k * pow(real_ppm , 2));//mz分数
                                float total_score = mz_score_weight * mz_score + (1 - mz_score_weight)*intensity_score;//最终分数，这里的intensity可以用上面算过的
                                (*ms2_itr)->EmplaceBackFa(Fa(*fragment_mz_itr , intensity , total_score , &fa_match->at(right_t)));
                            }
                            right_t++;//向右继续寻找
                        }
                        else{
                            break;
                        }
                    }
                    break;
                }
                else if(*fragment_mz_itr >= max_mz){
                    left = mid + 1;
                }
                else if(*fragment_mz_itr <= min_mz){
                    right = mid - 1;
                }
            }
        }
    }

    //删除无PA和FA的二级，只是取消引用而已，没有把真正的二级从内存中抹除
    for(auto ms2_itr = this->m_Ms2_vector_ptr.begin() ; ms2_itr != this->m_Ms2_vector_ptr.end();){
        //如果既没有PA又没有FA，则把对二级的引用去除
        if(((*ms2_itr)->GetPaCount() == 0) && ((*ms2_itr)->GetFaCount() == 0)){
            ms2_itr = this->m_Ms2_vector_ptr.erase(ms2_itr);
        }
        else{
            ms2_itr++;
        }
    }


    //如果所有二级都没有PA和FA，则清空对象
    if(this->m_Ms2_vector_ptr.size() == 0){
        this->EmptyObject();
    }
}


void Cardiolipin::splice()
{
    qDebug() << "Class Cardiolipin call this 'splice' function!!!";
}
