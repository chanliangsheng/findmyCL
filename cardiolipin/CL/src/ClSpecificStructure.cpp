#include <ClSpecificStructure.h>

float ClSpecificStructure::m_fragment_score_weight = 1;
float ClSpecificStructure::m_fa_consistency_score_weight = 1;
float ClSpecificStructure::m_pa_exist_score_weight = 1;
float ClSpecificStructure::m_fa_intensity_variance_score_weight = 1;

using namespace std;

ClSpecificStructure::ClSpecificStructure()
{
    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    PaNode* right_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;
    this->m_right_pa_ptr = right_pa_node;

    this->m_pa_exist = 0;
    this->m_fa_exist = 0;
    this->m_left_pa_ptr = nullptr;
    this->m_right_pa_ptr = nullptr;
    this->m_left_pa_ptr->m_left_fa_ptr = nullptr;
    this->m_left_pa_ptr->m_right_fa_ptr = nullptr;
    this->m_right_pa_ptr->m_left_fa_ptr = nullptr;
    this->m_right_pa_ptr->m_right_fa_ptr = nullptr;
    this->m_ms2 = nullptr;

    this->m_score = 0;
    this->m_total_intensity = 0;
}

ClSpecificStructure::ClSpecificStructure(Fa* fa_1_ptr, Fa* fa_2_ptr, Fa* fa_3_ptr, Fa* fa_4_ptr , Ms2* ms2_ptr)
{
    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    PaNode* right_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;
    this->m_right_pa_ptr = right_pa_node;

    //把4个FA的chain，unsaturation，oxygen的信息存储到m_fa_info中
    this->m_fa_info.insert({fa_1_ptr->GetChainLength() , fa_1_ptr->GetUnsaturation() , fa_1_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_2_ptr->GetChainLength() , fa_2_ptr->GetUnsaturation() , fa_2_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_3_ptr->GetChainLength() , fa_3_ptr->GetUnsaturation() , fa_3_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_4_ptr->GetChainLength() , fa_4_ptr->GetUnsaturation() , fa_4_ptr->GetOxygen()});

    this->m_ms2 = ms2_ptr;
    this->m_left_pa_ptr->m_left_fa_ptr = fa_1_ptr;
    this->m_left_pa_ptr->m_right_fa_ptr = fa_2_ptr;
    this->m_right_pa_ptr->m_left_fa_ptr = fa_3_ptr;
    this->m_right_pa_ptr->m_right_fa_ptr = fa_4_ptr;

    this->m_left_pa_ptr->m_pa_ptr = nullptr;
    this->m_right_pa_ptr->m_pa_ptr = nullptr;

    //pa或者fa是否存在
    this->m_pa_exist = 0;
    this->m_fa_exist = 1;

    //进行该拼接的打分
    this->score();

    //计算总强度
    this->CalculateTotalIntensity();
}

ClSpecificStructure::ClSpecificStructure(Pa *pa_1_ptr, Pa *pa_2_ptr, Ms2* ms2_ptr)
{
    this->m_ms2 = ms2_ptr;

    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    PaNode* right_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;
    this->m_right_pa_ptr = right_pa_node;

    this->m_left_pa_ptr->m_pa_ptr = pa_1_ptr;
    this->m_right_pa_ptr->m_pa_ptr = pa_2_ptr;
    this->m_left_pa_ptr->m_left_fa_ptr = nullptr;
    this->m_left_pa_ptr->m_right_fa_ptr = nullptr;
    this->m_right_pa_ptr->m_left_fa_ptr = nullptr;
    this->m_right_pa_ptr->m_right_fa_ptr = nullptr;

    //pa或者fa是否存在
    this->m_pa_exist = 1;
    this->m_fa_exist = 0;

    //把2个PA的chain，unsaturation，oxygen的信息存储到m_pa_info中
    this->m_pa_info.insert({pa_1_ptr->GetChainLength() , pa_1_ptr->GetUnsaturation() , pa_1_ptr->GetOxygen()});
    this->m_pa_info.insert({pa_2_ptr->GetChainLength() , pa_2_ptr->GetUnsaturation() , pa_2_ptr->GetOxygen()});

    //进行该拼接的打分
    this->score();
    //计算总强度
    this->CalculateTotalIntensity();
}

ClSpecificStructure::ClSpecificStructure(Pa *pa_1_ptr, Fa *fa_1_ptr, Fa *fa_2_ptr, Pa *pa_2_ptr, Fa *fa_3_ptr, Fa *fa_4_ptr, Ms2* ms2_ptr)
{
    this->m_ms2 = ms2_ptr;

    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    PaNode* right_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;
    this->m_right_pa_ptr = right_pa_node;

    this->m_left_pa_ptr->m_pa_ptr = pa_1_ptr;
    this->m_left_pa_ptr->m_left_fa_ptr = fa_1_ptr;
    this->m_left_pa_ptr->m_right_fa_ptr = fa_2_ptr;

    this->m_right_pa_ptr->m_pa_ptr = pa_2_ptr;
    this->m_right_pa_ptr->m_left_fa_ptr = fa_3_ptr;
    this->m_right_pa_ptr->m_right_fa_ptr = fa_4_ptr;

    //pa或者fa是否存在
    this->m_pa_exist = 1;
    this->m_fa_exist = 1;

    //把2个PA，4个FA的chain，unsaturation，oxygen的信息存储到m_pa_info中
    this->m_pa_info.insert({pa_1_ptr->GetChainLength() , pa_1_ptr->GetUnsaturation() , pa_1_ptr->GetOxygen()});
    this->m_pa_info.insert({pa_2_ptr->GetChainLength() , pa_2_ptr->GetUnsaturation() , pa_2_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_1_ptr->GetChainLength() , fa_1_ptr->GetUnsaturation() , fa_1_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_2_ptr->GetChainLength() , fa_2_ptr->GetUnsaturation() , fa_2_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_3_ptr->GetChainLength() , fa_3_ptr->GetUnsaturation() , fa_3_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_4_ptr->GetChainLength() , fa_4_ptr->GetUnsaturation() , fa_4_ptr->GetOxygen()});

    //进行该拼接的打分
    this->score();
    //计算总强度
    this->CalculateTotalIntensity();
}

float ClSpecificStructure::GetScore()
{
    return this->m_score;
}

void ClSpecificStructure::score()
{
    float fragment_score = 0;
    float fa_consistency_score = 0;


//    vector<int> match_info = {66,2,0};
//    if(this->GetTotalInfo() == match_info){
//        qDebug() << "stop";
//    }

    //如果这个组合的FA不存在，但是PA存在
    if(!this->m_fa_exist && this->m_pa_exist){
        //最后总分就是权重*1，因为pa存在就是1
        this->m_score = this->m_pa_exist_score_weight;
    }
    //如果FA存在，PA不存在
    else if(this->m_fa_exist && !this->m_pa_exist){
        //把4个FA的指针加入到哈希表中，这样可以把一样的FA合并，合并的结果是key：FA，value：FA出现的次数
        unordered_map<Fa*,float> fa_hashmap;
        fa_hashmap.insert({this->m_left_pa_ptr->m_left_fa_ptr , 1});

        auto find_itr = fa_hashmap.find(this->m_left_pa_ptr->m_right_fa_ptr);
        if(find_itr == fa_hashmap.end()){
            fa_hashmap.insert({this->m_left_pa_ptr->m_right_fa_ptr,1});
        }
        else{
            find_itr->second++;
        }

        find_itr = fa_hashmap.find(this->m_right_pa_ptr->m_left_fa_ptr);
        if(find_itr == fa_hashmap.end()){
            fa_hashmap.insert({this->m_right_pa_ptr->m_left_fa_ptr,1});
        }
        else{
            find_itr->second++;
        }

        find_itr = fa_hashmap.find(this->m_right_pa_ptr->m_right_fa_ptr);
        if(find_itr == fa_hashmap.end()){
            fa_hashmap.insert({this->m_right_pa_ptr->m_right_fa_ptr,1});
        }
        else{
            find_itr->second++;
        }


        //得到碎片分数
        for(auto itr = fa_hashmap.begin() ; itr != fa_hashmap.end() ; itr++){
            fragment_score += (*itr).first->GetScore();
        }


        //碎片分数是总分的平均值
        fragment_score = fragment_score / 4;

        //fa强度标准差分数
        vector<float> fa_score_vector;
        for(auto itr = fa_hashmap.begin() ; itr != fa_hashmap.end() ; itr++){
            for(float i = 0 ; i < itr->second ; i++){
                fa_score_vector.push_back(itr->first->GetScore() / itr->second);
            }
        }
        //求平方和
        float sum = std::accumulate(fa_score_vector.begin(), fa_score_vector.end(), 0.0, [fragment_score](float total, float value) {
            return total + std::pow(value - fragment_score, 2);
        });
        float fa_intensity_variance_score = 1 - sqrt(sum/4);//分数为1-标准差，方差越大，分数越小

        //fa一致性分数
        fa_consistency_score = float(this->m_fa_info.size()) / 4;


        //总分是碎片分数加上一致性分数
        this->m_score = this->m_fragment_score_weight*fragment_score + this->m_fa_consistency_score_weight*fa_consistency_score + this->m_fa_intensity_variance_score_weight*fa_intensity_variance_score;
    }
    //如果FA和PA都存在
    else if(this->m_fa_exist && this->m_pa_exist){
        //把4个FA的指针加入到哈希表中，这样可以把一样的FA合并，合并的结果是key：FA，value：FA出现的次数
        unordered_map<Fa*,float> fa_hashmap;
        fa_hashmap.insert({this->m_left_pa_ptr->m_left_fa_ptr , 1});

        auto find_itr = fa_hashmap.find(this->m_left_pa_ptr->m_right_fa_ptr);
        if(find_itr == fa_hashmap.end()){
            fa_hashmap.insert({this->m_left_pa_ptr->m_right_fa_ptr,1});
        }
        else{
            find_itr->second++;
        }

        find_itr = fa_hashmap.find(this->m_right_pa_ptr->m_left_fa_ptr);
        if(find_itr == fa_hashmap.end()){
            fa_hashmap.insert({this->m_right_pa_ptr->m_left_fa_ptr,1});
        }
        else{
            find_itr->second++;
        }

        find_itr = fa_hashmap.find(this->m_right_pa_ptr->m_right_fa_ptr);
        if(find_itr == fa_hashmap.end()){
            fa_hashmap.insert({this->m_right_pa_ptr->m_right_fa_ptr,1});
        }
        else{
            find_itr->second++;
        }


        //得到碎片分数
        for(auto itr = fa_hashmap.begin() ; itr != fa_hashmap.end() ; itr++){
            fragment_score += (*itr).first->GetScore();
        }

        //碎片分数是总分的平均值
        fragment_score = fragment_score / 4;

        //fa强度标准差分数
        vector<float> fa_score_vector;
        for(auto itr = fa_hashmap.begin() ; itr != fa_hashmap.end() ; itr++){
            for(float i = 0 ; i < itr->second ; i++){
                fa_score_vector.push_back(itr->first->GetScore() / itr->second);
            }
        }
        //求平方和
        float sum = std::accumulate(fa_score_vector.begin(), fa_score_vector.end(), 0.0, [fragment_score](float total, float value) {
            return total + std::pow(value - fragment_score, 2);
        });
        float fa_intensity_variance_score = 1 - sqrt(sum/4);//分数为1-标准差，方差越大，分数越小

        //fa一致性分数
        fa_consistency_score = float(this->m_fa_info.size()) / 4;
        //总分是碎片分数+一致性分数+pa存在性权重
        this->m_score = this->m_fragment_score_weight*fragment_score + this->m_fa_consistency_score_weight*fa_consistency_score + this->m_pa_exist_score_weight + this->m_fa_intensity_variance_score_weight*fa_intensity_variance_score;
    }
}

float ClSpecificStructure::GetTotalIntensity()
{
    return this->m_total_intensity;
}

void ClSpecificStructure::CalculateTotalIntensity()
{
    if(!this->m_fa_exist && this->m_pa_exist){
        //把2个PA指针加入set中，因为可以把一样的PA合并
        set<Pa*> pa_set;
        pa_set.insert(this->m_left_pa_ptr->m_pa_ptr);
        pa_set.insert(this->m_right_pa_ptr->m_pa_ptr);
        for(auto itr = pa_set.begin() ; itr != pa_set.end() ; itr++){
            this->m_total_intensity = this->m_total_intensity + (*itr)->GetIntensity();
        }
    }
    else if(this->m_fa_exist && !this->m_pa_exist){
        //把4个FA指针加入set中，因为可以把一样的FA合并
        set<Fa*> fa_set;
        fa_set.insert(this->m_left_pa_ptr->m_left_fa_ptr);
        fa_set.insert(this->m_left_pa_ptr->m_right_fa_ptr);
        fa_set.insert(this->m_right_pa_ptr->m_left_fa_ptr);
        fa_set.insert(this->m_right_pa_ptr->m_right_fa_ptr);

        //加分数
        for(auto itr = fa_set.begin();itr != fa_set.end();itr++){
            this->m_total_intensity += (*itr)->GetIntensity();
        }
    }
    else{
        set<Pa*> pa_set;
        set<Fa*> fa_set;
        pa_set.insert(this->m_left_pa_ptr->m_pa_ptr);
        pa_set.insert(this->m_right_pa_ptr->m_pa_ptr);
        fa_set.insert(this->m_left_pa_ptr->m_left_fa_ptr);
        fa_set.insert(this->m_left_pa_ptr->m_right_fa_ptr);
        fa_set.insert(this->m_right_pa_ptr->m_left_fa_ptr);
        fa_set.insert(this->m_right_pa_ptr->m_right_fa_ptr);
        for(auto itr = pa_set.begin() ; itr != pa_set.end() ; itr++){
            this->m_total_intensity = this->m_total_intensity + (*itr)->GetIntensity();
        }
        for(auto itr = fa_set.begin();itr != fa_set.end();itr++){
            this->m_total_intensity += (*itr)->GetIntensity();
        }
    }
}

float ClSpecificStructure::GetMs2TotalIntensity()
{
    return this->m_ms2->GetTotalIntensity();
}

bool ClSpecificStructure::operator==(const ClSpecificStructure &other)
{
//    //如果两个Cl来源于同一个二级，则不可能相同，这是由拼接的过程决定的
//    if(this->m_ms2 == other.m_ms2){
//        return false;
//    }
    //如果两者的FA都不存在，则比较PA是否相同
    if(!this->m_fa_exist && !other.m_fa_exist){
        //如果两者的PA信息相同，则把this的分数更新为最大值
        if(this->m_pa_info == other.m_pa_info){
            this->m_score = max(this->m_score , other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果两者的PA都不存在，则比较FA是否相同
    else if(!this->m_pa_exist && !other.m_pa_exist){
        //如果两者的FA信息相同，则把this的分数更新为最大值
        if(this->m_fa_info == other.m_fa_info){
            this->m_score = max(this->m_score , other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果两者的PA，FA都存在
    else if(this->m_fa_exist && other.m_fa_exist && this->m_pa_exist && other.m_pa_exist){
        //如果两者的PA和FA信息相同，则把this的分数更新为最大值
        if((this->m_pa_info == other.m_pa_info) && (this->m_fa_info == other.m_fa_info)){
            this->m_score = max(this->m_score , other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果一个只有PA，另一个既有PA，又有FA
    else if((!this->m_fa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_fa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的PA信息相同，把this的内容换成新的内容
        if(this->m_pa_info == other.m_pa_info){
            if(!this->m_fa_exist){
                *this = other;
                return true;
            }
            else if(!other.m_fa_exist){
                return true;
            }
        }
        else{
            return false;
        }
    }
    //如果一个只有FA，另一个既有PA，又有FA
    else if((!this->m_pa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_pa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的FA信息相同，把this的内容换成新的内容，并且把分数换为最大值
        if(this->m_fa_info == other.m_fa_info){
            float max_score = max(this->m_score , other.m_score);
            if(!this->m_pa_exist){
                *this = other;
                this->m_score = max_score;
                return true;
            }
            else if(!other.m_pa_exist){
                this->m_score = max_score;
                return true;
            }
        }
        else{
            return false;
        }
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other.m_fa_exist) || (!other.m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        ClSpecificStructure* object_only_fa;
        ClSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            ClSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            ClSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        //判定object_only_fa的4个FA是否能判定成object_only_pa的2个PA
        //用vector存储4个FA
        vector<Fa*> fa_ptr_vector;
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_right_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_pa_ptr->m_left_fa_ptr);
        //哈系法解两数之和
        bool splice_success = 0;
        multimap<unsigned int,Fa*> hash_map;//以链长为键，对应的Fa为值构建多重哈希表
        for(auto fa_ptr_vector_itr = fa_ptr_vector.begin() ; fa_ptr_vector_itr!= fa_ptr_vector.end() ; fa_ptr_vector_itr++){
            //寻找是否有2个FA可以组成第一个PA
            auto fist_pa_search_pair_itr = hash_map.equal_range(object_only_pa->m_left_pa_ptr->m_pa_ptr->GetChainLength() - (*fa_ptr_vector_itr)->GetChainLength());
            //如果没有找到
            if(fist_pa_search_pair_itr.first == fist_pa_search_pair_itr.second){

            }
            else{
                for(auto itr = fist_pa_search_pair_itr.first ; itr!=fist_pa_search_pair_itr.second ; itr++){
                  if(((*fa_ptr_vector_itr)->GetUnsaturation() + itr->second->GetUnsaturation() == object_only_pa->m_left_pa_ptr->m_pa_ptr->GetUnsaturation())&&
                     ((*fa_ptr_vector_itr)->GetOxygen() + itr->second->GetOxygen() == object_only_pa->m_left_pa_ptr->m_pa_ptr->GetOxygen())){

                    splice_success = 1;//设定为真
                    //将配对上第1个PA的2个FA的信息加入到只有PA的对象中
                    object_only_pa->m_left_pa_ptr->m_left_fa_ptr = *fa_ptr_vector_itr;
                    object_only_pa->m_left_pa_ptr->m_right_fa_ptr = itr->second;
                    //去除那两个拼接成第1个PA的两个FA
                    fa_ptr_vector.erase(fa_ptr_vector_itr);
                    for(auto vector_itr = fa_ptr_vector.begin();vector_itr!=fa_ptr_vector.end();){
                      if(*vector_itr == itr->second){
                        vector_itr = fa_ptr_vector.erase(vector_itr);
                      }
                      else{
                        vector_itr++;
                      }
                    }
                    //赋值给第2个PA
                    object_only_pa->m_right_pa_ptr->m_left_fa_ptr = fa_ptr_vector[0];
                    object_only_pa->m_right_pa_ptr->m_right_fa_ptr = fa_ptr_vector[1];

                    if(this == object_only_fa){
                      *this = *object_only_pa;
                    }
                    this->update();//更新分数等信息
                    return true;
                  }
                }
            }
            hash_map.insert({(*fa_ptr_vector_itr)->GetChainLength() , *fa_ptr_vector_itr});//把链长作为键，Fa的指针作为值加入到哈希表
        }
        //如果无法找到4个FA可以拼接成这2个PA，则返回false
        if(!splice_success){
            return false;
        }
    }
    return false;
}

shared_ptr<ClSpecificStructure> ClSpecificStructure::merge(ClSpecificStructure* other , bool merge_m_h_m_2h)
{
//    //如果不是用于合并M-H和M-2H
//    if(!merge_m_h_m_2h){
//        //如果两个Cl来源于同一个二级，则不可能相同，这是由拼接的过程决定的
//        if(this->m_ms2 == other->m_ms2){
//            return nullptr;
//        }
//    }
    //如果两者的FA都不存在，则比较PA是否相同
    if(!this->m_fa_exist && !other->m_fa_exist){
        //如果两者的PA信息相同，则返回一个新的对象的共享指针
        if(this->m_pa_info == other->m_pa_info){
            shared_ptr<ClSpecificStructure> new_object = make_shared<ClSpecificStructure>(*this);
            new_object->m_score = max(this->m_score , other->m_score);
            return new_object;
        }
        else{
            return nullptr;
        }
    }
    //如果两者的PA都不存在，则比较FA是否相同
    else if(!this->m_pa_exist && !other->m_pa_exist){
        //如果两者的FA信息相同，则返回一个新的对象的指针
        if(this->m_fa_info == other->m_fa_info){
            shared_ptr<ClSpecificStructure> new_object = make_shared<ClSpecificStructure>(*this);
            new_object->m_score = max(this->m_score , other->m_score);
            return new_object;
        }
        else{
            return nullptr;
        }
    }
    //如果两者的PA，FA都存在
    else if(this->m_fa_exist && other->m_fa_exist && this->m_pa_exist && other->m_pa_exist){
        //如果两者的PA和FA信息相同，则返回一个新的对象的共享指针
        if((this->m_pa_info == other->m_pa_info) && (this->m_fa_info == other->m_fa_info)){
            shared_ptr<ClSpecificStructure> new_object = make_shared<ClSpecificStructure>(*this);
            new_object->m_score = max(this->m_score , other->m_score);
            return new_object;
        }
        else{
            return nullptr;
        }
    }
    //如果一个只有PA，另一个既有PA，又有FA
    else if((!this->m_fa_exist && other->m_fa_exist && other->m_pa_exist) || (!other->m_fa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的PA信息相同，把this的内容换成新的内容
        if(this->m_pa_info == other->m_pa_info){
            if(!this->m_fa_exist){
                shared_ptr<ClSpecificStructure> new_object = make_shared<ClSpecificStructure>(*this);
                return new_object;
            }
            else if(!other->m_fa_exist){
                shared_ptr<ClSpecificStructure> new_object = make_shared<ClSpecificStructure>(*this);
                return new_object;
            }
        }
        else{
            return nullptr;
        }
    }
    //如果一个只有FA，另一个既有PA，又有FA
    else if((!this->m_pa_exist && other->m_fa_exist && other->m_pa_exist) || (!other->m_pa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的FA信息相同，把this的内容换成新的内容，并且把分数换为最大值
        if(this->m_fa_info == other->m_fa_info){
            float max_score = max(this->m_score , other->m_score);
            if(!this->m_pa_exist){
                shared_ptr<ClSpecificStructure> new_object = make_shared<ClSpecificStructure>(*other);
                new_object->m_score = max_score;
                return new_object;
            }
            else if(!other->m_pa_exist){
                shared_ptr<ClSpecificStructure> new_object = make_shared<ClSpecificStructure>(*this);
                new_object->m_score = max_score;
                return new_object;
            }
        }
        else{
            return nullptr;
        }
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other->m_fa_exist) || (!other->m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        ClSpecificStructure* object_only_fa;
        shared_ptr<ClSpecificStructure> object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            object_only_pa = make_shared<ClSpecificStructure>(*other);
        }
        else{
            object_only_fa = other;
            object_only_pa = make_shared<ClSpecificStructure>(*this);
        }

        //判定object_only_fa的4个FA是否能判定成object_only_pa的2个PA
        //用vector存储4个FA
        vector<Fa*> fa_ptr_vector;
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_right_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_pa_ptr->m_left_fa_ptr);
        //哈系法解两数之和
        bool splice_success = 0;
        multimap<unsigned int,Fa*> hash_map;//以链长为键，对应的Fa为值构建多重哈希表
        for(auto fa_ptr_vector_itr = fa_ptr_vector.begin() ; fa_ptr_vector_itr!= fa_ptr_vector.end() ; fa_ptr_vector_itr++){
            //寻找是否有2个FA可以组成第一个PA
            auto fist_pa_search_pair_itr = hash_map.equal_range(object_only_pa->m_left_pa_ptr->m_pa_ptr->GetChainLength() - (*fa_ptr_vector_itr)->GetChainLength());
            //如果没有找到
            if(fist_pa_search_pair_itr.first == fist_pa_search_pair_itr.second){

            }
            else{
                //遍历可能配对上的结果，如果配对上了，multimap可能没有配对上，它会取d大于母包键的范围
                for(auto itr = fist_pa_search_pair_itr.first ; itr!=fist_pa_search_pair_itr.second ; itr++){
                    if(((*fa_ptr_vector_itr)->GetUnsaturation() + itr->second->GetUnsaturation() == object_only_pa->m_left_pa_ptr->m_pa_ptr->GetUnsaturation())&&
                       ((*fa_ptr_vector_itr)->GetOxygen() + itr->second->GetOxygen() == object_only_pa->m_left_pa_ptr->m_pa_ptr->GetOxygen())){

                        splice_success = 1;//设定为真
                        //将配对上第1个PA的2个FA的信息加入到只有PA的对象中
                        object_only_pa->m_left_pa_ptr->m_left_fa_ptr = *fa_ptr_vector_itr;
                        object_only_pa->m_left_pa_ptr->m_right_fa_ptr = itr->second;
                        //去除那两个拼接成第1个PA的两个FA
                        fa_ptr_vector.erase(fa_ptr_vector_itr);
                        for(auto vector_itr = fa_ptr_vector.begin();vector_itr!=fa_ptr_vector.end();){
                            if(*vector_itr == itr->second){
                                vector_itr = fa_ptr_vector.erase(vector_itr);
                            }
                            else{
                                vector_itr++;
                            }
                        }
                        //赋值给第2个PA
                        object_only_pa->m_right_pa_ptr->m_left_fa_ptr = fa_ptr_vector[0];
                        object_only_pa->m_right_pa_ptr->m_right_fa_ptr = fa_ptr_vector[1];
                        object_only_pa->update();//更新分数等信息
                        return object_only_pa;
                    }
                }
            }
            hash_map.insert({(*fa_ptr_vector_itr)->GetChainLength() , *fa_ptr_vector_itr});//把链长作为键，Fa的指针作为值加入到哈希表
        }
        //如果无法找到4个FA可以拼接成这2个PA，则返回false
        if(!splice_success){
            return nullptr;
        }
    }
    return nullptr;
}

void ClSpecificStructure::update()
{
    //清空FA和PA信息
    set<array<unsigned int,3>>().swap(this->m_pa_info);
    set<array<unsigned int,3>>().swap(this->m_fa_info);

    //更新PA信息
    if((this->m_left_pa_ptr->m_pa_ptr != nullptr) && (this->m_right_pa_ptr->m_pa_ptr != nullptr)){
        this->m_pa_exist = 1;
        this->m_pa_info.insert({this->m_left_pa_ptr->m_pa_ptr->GetChainLength() , this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation() , this->m_left_pa_ptr->m_pa_ptr->GetOxygen()});
        this->m_pa_info.insert({this->m_right_pa_ptr->m_pa_ptr->GetChainLength() , this->m_right_pa_ptr->m_pa_ptr->GetUnsaturation() , this->m_right_pa_ptr->m_pa_ptr->GetOxygen()});
    }
    //更新FA信息
    if((this->m_left_pa_ptr->m_left_fa_ptr != nullptr) && (this->m_right_pa_ptr->m_right_fa_ptr != nullptr) && (this->m_left_pa_ptr->m_right_fa_ptr != nullptr) && (this->m_right_pa_ptr->m_left_fa_ptr != nullptr)){
        this->m_fa_exist = 1;
        this->m_fa_info.insert({this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() , this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() , this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()});
        this->m_fa_info.insert({this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength() , this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation() , this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()});
        this->m_fa_info.insert({this->m_right_pa_ptr->m_left_fa_ptr->GetChainLength() , this->m_right_pa_ptr->m_left_fa_ptr->GetUnsaturation() , this->m_right_pa_ptr->m_left_fa_ptr->GetOxygen()});
        this->m_fa_info.insert({this->m_right_pa_ptr->m_right_fa_ptr->GetChainLength() , this->m_right_pa_ptr->m_right_fa_ptr->GetUnsaturation() , this->m_right_pa_ptr->m_right_fa_ptr->GetOxygen()});
    }

    //更新得分信息
    this->score();
}

bool ClSpecificStructure::StrictMerge(ClSpecificStructure &other)
{
    //如果两者的FA都不存在，则比较PA是否相同
    if(!this->m_fa_exist && !other.m_fa_exist){
        //如果两者的PA信息相同，则把this的分数更新为最大值
        if(this->m_pa_info == other.m_pa_info){
            this->m_score = max(this->m_score , other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果两者的PA都不存在，则比较FA是否相同
    else if(!this->m_pa_exist && !other.m_pa_exist){
        //如果两者的FA信息相同，则把this的分数更新为新值
        if(this->m_fa_info == other.m_fa_info){
            this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果两者的PA，FA都存在
    else if(this->m_fa_exist && other.m_fa_exist && this->m_pa_exist && other.m_pa_exist){
        //如果两者的PA和FA信息相同，则把this的分数更新为最大值
        if((this->m_pa_info == other.m_pa_info) && (this->m_fa_info == other.m_fa_info)){
            this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果一个只有PA，另一个既有PA，又有FA
    else if((!this->m_fa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_fa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的PA信息相同，把this的内容换成新的内容
        if(this->m_pa_info == other.m_pa_info){
            if(!this->m_fa_exist){
                *this = other;
                return true;
            }
            else if(!other.m_fa_exist){
                return true;
            }
        }
        else{
            return false;
        }
    }
    //如果一个只有FA，另一个既有PA，又有FA
    else if((!this->m_pa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_pa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的FA信息相同，把this的内容换成新的内容，并且把分数换为最大值
        if(this->m_fa_info == other.m_fa_info){
            if(!this->m_pa_exist){
                *this = other;
                this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
                return true;
            }
            else if(!other.m_pa_exist){
                this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
                return true;
            }
        }
        else{
            return false;
        }
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other.m_fa_exist) || (!other.m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        ClSpecificStructure* object_only_fa;
        ClSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            ClSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            ClSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        //判定object_only_fa的4个FA是否能判定成object_only_pa的2个PA
        //用vector存储4个FA
        vector<Fa*> fa_ptr_vector;
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_right_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_pa_ptr->m_left_fa_ptr);
        //哈系法解两数之和
        bool splice_success = 0;
        multimap<unsigned int,Fa*> hash_map;//以链长为键，对应的Fa为值构建多重哈希表
        for(auto fa_ptr_vector_itr = fa_ptr_vector.begin() ; fa_ptr_vector_itr!= fa_ptr_vector.end() ; fa_ptr_vector_itr++){
            //寻找是否有2个FA可以组成第一个PA
            auto fist_pa_search_pair_itr = hash_map.equal_range(object_only_pa->m_left_pa_ptr->m_pa_ptr->GetChainLength() - (*fa_ptr_vector_itr)->GetChainLength());
            //如果没有找到
            if(fist_pa_search_pair_itr.first == fist_pa_search_pair_itr.second){

            }
            else{
                for(auto itr = fist_pa_search_pair_itr.first ; itr!=fist_pa_search_pair_itr.second ; itr++){
                  if(((*fa_ptr_vector_itr)->GetUnsaturation() + itr->second->GetUnsaturation() == object_only_pa->m_left_pa_ptr->m_pa_ptr->GetUnsaturation())&&
                     ((*fa_ptr_vector_itr)->GetOxygen() + itr->second->GetOxygen() == object_only_pa->m_left_pa_ptr->m_pa_ptr->GetOxygen())){

                    splice_success = 1;//设定为真
                    //将配对上第1个PA的2个FA的信息加入到只有PA的对象中
                    object_only_pa->m_left_pa_ptr->m_left_fa_ptr = *fa_ptr_vector_itr;
                    object_only_pa->m_left_pa_ptr->m_right_fa_ptr = itr->second;
                    //去除那两个拼接成第1个PA的两个FA
                    fa_ptr_vector.erase(fa_ptr_vector_itr);
                    for(auto vector_itr = fa_ptr_vector.begin();vector_itr!=fa_ptr_vector.end();){
                      if(*vector_itr == itr->second){
                        vector_itr = fa_ptr_vector.erase(vector_itr);
                      }
                      else{
                        vector_itr++;
                      }
                    }
                    //赋值给第2个PA
                    object_only_pa->m_right_pa_ptr->m_left_fa_ptr = fa_ptr_vector[0];
                    object_only_pa->m_right_pa_ptr->m_right_fa_ptr = fa_ptr_vector[1];

                    if(this == object_only_fa){
                      *this = *object_only_pa;
                    }
                    this->update();//更新分数等信息
                    return true;
                  }
                }
            }
            hash_map.insert({(*fa_ptr_vector_itr)->GetChainLength() , *fa_ptr_vector_itr});//把链长作为键，Fa的指针作为值加入到哈希表
        }
        //如果无法找到4个FA可以拼接成这2个PA，则返回false
        if(!splice_success){
            return false;
        }
    }
    return false;
}

bool ClSpecificStructure::FlexibleMerge(ClSpecificStructure &other)
{
//    //如果两个Cl来源于同一个二级，则不可能相同，这是由拼接的过程决定的
//    if(this->m_ms2 == other.m_ms2){
//        return false;
//    }
    //如果两者的FA都不存在，则比较PA是否相同
    if(!this->m_fa_exist && !other.m_fa_exist){
        //如果两者的PA信息相同，则把this的分数更新为最大值
        if(this->m_pa_info == other.m_pa_info){
            this->m_score = max(this->m_score , other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果两者的PA都不存在，则比较FA是否相同
    else if(!this->m_pa_exist && !other.m_pa_exist){
        //如果两者的FA信息相同，则把this的分数更新为新值
        if(this->m_fa_info == other.m_fa_info){
            this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果两者的PA，FA都存在
    else if(this->m_fa_exist && other.m_fa_exist && this->m_pa_exist && other.m_pa_exist){
        //如果两者的PA和FA信息相同，则把this的分数更新为最大值
//        if((this->m_pa_info == other.m_pa_info) && (this->m_fa_info == other.m_fa_info)){
//            this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
//            return true;
//        }
        //如果FA的信息相同
        if(this->m_fa_info == other.m_fa_info){
            //如果PA的信息也相同
            if(this->m_pa_info == other.m_pa_info){
                this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
                return true;
            }
            //如果PA的信息不相同
            else{
                this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
                this->m_flexible_mode_dif_pa_merge.emplace_back(&other);//把跟this的PA不同的Cl的信息存储到m_flexible_mode_dif_pa_merge中
                return true;
            }
        }
        else{
            return false;
        }
    }
    //如果一个只有PA，另一个既有PA，又有FA
    else if((!this->m_fa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_fa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的PA信息相同，把this的内容换成新的内容
        if(this->m_pa_info == other.m_pa_info){
            if(!this->m_fa_exist){
                *this = other;
                return true;
            }
            else if(!other.m_fa_exist){
                return true;
            }
        }
        else{
            return false;
        }
    }
    //如果一个只有FA，另一个既有PA，又有FA
    else if((!this->m_pa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_pa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的FA信息相同，把this的内容换成新的内容，并且把分数换为最大值
        if(this->m_fa_info == other.m_fa_info){
            if(!this->m_pa_exist){
                *this = other;
                this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
                return true;
            }
            else if(!other.m_pa_exist){
                this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
                return true;
            }
        }
        else{
            return false;
        }
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other.m_fa_exist) || (!other.m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        ClSpecificStructure* object_only_fa;
        ClSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            ClSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            ClSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        //判定object_only_fa的4个FA是否能判定成object_only_pa的2个PA
        //用vector存储4个FA
        vector<Fa*> fa_ptr_vector;
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_right_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_pa_ptr->m_left_fa_ptr);
        //哈系法解两数之和
        bool splice_success = 0;
        multimap<unsigned int,Fa*> hash_map;//以链长为键，对应的Fa为值构建多重哈希表
        for(auto fa_ptr_vector_itr = fa_ptr_vector.begin() ; fa_ptr_vector_itr!= fa_ptr_vector.end() ; fa_ptr_vector_itr++){
            //寻找是否有2个FA可以组成第一个PA
            auto fist_pa_search_pair_itr = hash_map.equal_range(object_only_pa->m_left_pa_ptr->m_pa_ptr->GetChainLength() - (*fa_ptr_vector_itr)->GetChainLength());
            //如果没有找到
            if(fist_pa_search_pair_itr.first == fist_pa_search_pair_itr.second){

            }
            else{
                for(auto itr = fist_pa_search_pair_itr.first ; itr!=fist_pa_search_pair_itr.second ; itr++){
                  if(((*fa_ptr_vector_itr)->GetUnsaturation() + itr->second->GetUnsaturation() == object_only_pa->m_left_pa_ptr->m_pa_ptr->GetUnsaturation())&&
                     ((*fa_ptr_vector_itr)->GetOxygen() + itr->second->GetOxygen() == object_only_pa->m_left_pa_ptr->m_pa_ptr->GetOxygen())){

                    splice_success = 1;//设定为真
                    //将配对上第1个PA的2个FA的信息加入到只有PA的对象中
                    object_only_pa->m_left_pa_ptr->m_left_fa_ptr = *fa_ptr_vector_itr;
                    object_only_pa->m_left_pa_ptr->m_right_fa_ptr = itr->second;
                    //去除那两个拼接成第1个PA的两个FA
                    fa_ptr_vector.erase(fa_ptr_vector_itr);
                    for(auto vector_itr = fa_ptr_vector.begin();vector_itr!=fa_ptr_vector.end();){
                      if(*vector_itr == itr->second){
                        vector_itr = fa_ptr_vector.erase(vector_itr);
                      }
                      else{
                        vector_itr++;
                      }
                    }
                    //赋值给第2个PA
                    object_only_pa->m_right_pa_ptr->m_left_fa_ptr = fa_ptr_vector[0];
                    object_only_pa->m_right_pa_ptr->m_right_fa_ptr = fa_ptr_vector[1];

                    if(this == object_only_fa){
                      *this = *object_only_pa;
                    }
                    this->update();//更新分数等信息
                    return true;
                  }
                }
            }
            hash_map.insert({(*fa_ptr_vector_itr)->GetChainLength() , *fa_ptr_vector_itr});//把链长作为键，Fa的指针作为值加入到哈希表
        }
        //如果无法找到4个FA可以拼接成这2个PA，则返回false
        if(!splice_success){
            return false;
        }
    }
    return false;
}

QString ClSpecificStructure::ShowInfo()
{
    QString message;
    //如果PA和FA都存在
    if(this->m_fa_exist && this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_pa_ptr->GetAdditiveForm()) + ");";
        QString Pa2_message = QString::number(this->m_right_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_right_pa_ptr->m_pa_ptr->GetAdditiveForm()) + ")";
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_left_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_right_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_3_message = QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_right_pa_ptr->m_left_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_4_message = QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_right_pa_ptr->m_right_fa_ptr->GetAdditiveForm()) + ");";
        message =  message + "This is Cl : have both PA and FA :" + "PA1:" + Pa1_message + " FA1:" +  Fa_1_message  + "  FA2 :" + Fa_2_message + "   PA2" + Pa2_message + "  FA3 :" + Fa_3_message + "   FA4 :" + Fa_4_message + "   score :" + QString::number(this->m_score);
    }
    //如果只有FA存在
    else if(this->m_fa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_left_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_right_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_3_message = QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_right_pa_ptr->m_left_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_4_message = QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetScore())  + "add_form:" + QString::fromStdString(this->m_right_pa_ptr->m_right_fa_ptr->GetAdditiveForm()) + ");";
        message =  message + "This is Cl : have only FA :" "FA1:" +  Fa_1_message  + "  FA2 :" + Fa_2_message + "   FA3 :" + Fa_3_message + "   FA4 :" + Fa_4_message + "   score :" + QString::number(this->m_score);
    }
    //如果只有PA存在
    else if(this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetScore()) + QString::fromStdString(this->m_left_pa_ptr->m_pa_ptr->GetAdditiveForm()) + ");";
        QString Pa2_message = QString::number(this->m_right_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_right_pa_ptr->m_pa_ptr->GetAdditiveForm()) + ")";
        message =  message + "This is Cl : have only PA :" + "PA1:" + Pa1_message + "   PA2" + Pa2_message + "  score :" + QString::number(this->m_score);
    }

    return message;
}

QString ClSpecificStructure::ShowSimpleInfo()
{
    QString message;
    //如果PA和FA都存在
    if(this->m_fa_exist && this->m_pa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen());
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen());
        QString Fa_3_message = QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetOxygen());
        QString Fa_4_message = QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetOxygen());
        message =  message +  Fa_1_message + "/" + Fa_2_message + "/" + Fa_3_message + "/" + Fa_4_message;
    }
    //如果只有FA存在
    else if(this->m_fa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen());
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen());
        QString Fa_3_message = QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_pa_ptr->m_left_fa_ptr->GetOxygen());
        QString Fa_4_message = QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_right_pa_ptr->m_right_fa_ptr->GetOxygen());
        message =  message +  Fa_1_message + "/" + Fa_2_message + "/" + Fa_3_message + "/" + Fa_4_message;
    }
    //如果只有PA存在
    else if(this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen());
        QString Pa2_message = QString::number(this->m_right_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_right_pa_ptr->m_pa_ptr->GetOxygen());
        message =  message + Pa1_message + "/" + Pa2_message;
    }

    return message;
}

std::vector<int> ClSpecificStructure::GetTotalInfo()
{
    return {this->GetTotalChainLength() , this->GetTotalUnsaturation() , this->GetTotalOxygen()};
}

int ClSpecificStructure::GetTotalChainLength()
{
    if(this->m_pa_exist){
        return this->m_left_pa_ptr->m_pa_ptr->GetChainLength() + this->m_right_pa_ptr->m_pa_ptr->GetChainLength();
    }
    else{
        return this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() + this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength() + this->m_right_pa_ptr->m_left_fa_ptr->GetChainLength() + this->m_right_pa_ptr->m_right_fa_ptr->GetChainLength();
    }
}

int ClSpecificStructure::GetTotalUnsaturation()
{
    if(this->m_pa_exist){
        return this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation() + this->m_right_pa_ptr->m_pa_ptr->GetUnsaturation();
    }
    else{
        return this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() + this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation() + this->m_right_pa_ptr->m_left_fa_ptr->GetUnsaturation() + this->m_right_pa_ptr->m_right_fa_ptr->GetUnsaturation();
    }
}

int ClSpecificStructure::GetTotalOxygen()
{
    if(this->m_pa_exist){
        return this->m_left_pa_ptr->m_pa_ptr->GetOxygen() + this->m_right_pa_ptr->m_pa_ptr->GetOxygen();
    }
    else{
        return this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen() + this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen() + this->m_right_pa_ptr->m_left_fa_ptr->GetOxygen() + this->m_right_pa_ptr->m_right_fa_ptr->GetOxygen();
    }
}
