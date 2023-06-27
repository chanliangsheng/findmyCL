#include <MlclSpecificStructure.h>

float MlclSpecificStructure::m_fragment_score_weight = 1;
float MlclSpecificStructure::m_fa_consistency_score_weight = 1;
float MlclSpecificStructure::m_pa_exist_score_weight = 1;

using namespace std;

MlclSpecificStructure::MlclSpecificStructure()
{
    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;

    this->m_pa_exist = 0;
    this->m_fa_exist = 0;
    this->m_left_pa_ptr = nullptr;
    this->m_left_pa_ptr->m_left_fa_ptr = nullptr;
    this->m_left_pa_ptr->m_right_fa_ptr = nullptr;
    this->m_right_fa_ptr = nullptr;
    this->m_ms2 = nullptr;
    this->m_score = 0;
}

MlclSpecificStructure::MlclSpecificStructure(Fa *fa_1_ptr, Fa *fa_2_ptr, Fa *fa_3_ptr, Ms2 *ms2_ptr)
{
    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;


    //把3个FA的chain，unsaturation，oxygen的信息存储到m_fa_info中
    this->m_fa_info.insert({fa_1_ptr->GetChainLength() , fa_1_ptr->GetUnsaturation() , fa_1_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_2_ptr->GetChainLength() , fa_2_ptr->GetUnsaturation() , fa_2_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_3_ptr->GetChainLength() , fa_3_ptr->GetUnsaturation() , fa_3_ptr->GetOxygen()});

    this->m_ms2 = ms2_ptr;
    this->m_left_pa_ptr->m_left_fa_ptr = fa_1_ptr;
    this->m_left_pa_ptr->m_right_fa_ptr = fa_2_ptr;
    this->m_right_fa_ptr = fa_3_ptr;
    this->m_left_pa_ptr->m_pa_ptr = nullptr;

    //pa或者fa是否存在
    this->m_pa_exist = 0;
    this->m_fa_exist = 1;

    //进行打分
    this->score();
}

MlclSpecificStructure::MlclSpecificStructure(Pa *pa_1_ptr, Fa *fa_1_ptr, Fa *fa_2_ptr, Fa *fa_3_ptr, Ms2 *ms2_ptr)
{
    this->m_ms2 = ms2_ptr;

    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;

    this->m_left_pa_ptr->m_pa_ptr = pa_1_ptr;
    this->m_left_pa_ptr->m_left_fa_ptr = fa_1_ptr;
    this->m_left_pa_ptr->m_right_fa_ptr = fa_2_ptr;

    this->m_right_fa_ptr = fa_3_ptr;

    //pa或者fa是否存在
    this->m_pa_exist = 1;
    this->m_fa_exist = 1;
    //把1个PA，3个FA的chain，unsaturation，oxygen的信息存储到m_pa_info中
    this->m_pa_info.insert({pa_1_ptr->GetChainLength() , pa_1_ptr->GetUnsaturation() , pa_1_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_1_ptr->GetChainLength() , fa_1_ptr->GetUnsaturation() , fa_1_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_2_ptr->GetChainLength() , fa_2_ptr->GetUnsaturation() , fa_2_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_3_ptr->GetChainLength() , fa_3_ptr->GetUnsaturation() , fa_3_ptr->GetOxygen()});

    //进行打分
    this->score();
}

MlclSpecificStructure::MlclSpecificStructure(Pa *pa_1_ptr, Fa *fa_3_ptr, Ms2 *ms2_ptr)
{
    this->m_ms2 = ms2_ptr;
    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;
    this->m_left_pa_ptr->m_left_fa_ptr = nullptr;
    this->m_left_pa_ptr->m_right_fa_ptr = nullptr;

    this->m_left_pa_ptr->m_pa_ptr = pa_1_ptr;
    this->m_right_fa_ptr = fa_3_ptr;

    //pa或者fa是否存在
    this->m_pa_exist = 1;
    this->m_fa_exist = 0;

    //把1个PA的chain，unsaturation，oxygen的信息存储到m_pa_info中
    this->m_pa_info.insert({pa_1_ptr->GetChainLength() , pa_1_ptr->GetUnsaturation() , pa_1_ptr->GetOxygen()});

    this->m_fa_info.insert({fa_3_ptr->GetChainLength() , fa_3_ptr->GetUnsaturation() , fa_3_ptr->GetOxygen()});

    //进行打分
    this->score();
}

void MlclSpecificStructure::score()
{
    float fragment_score = 0;
//    float pa_exist_score = 0;
    float fa_consistency_score = 0;

    //如果这个组合的FA不存在，但是PA存在
    if(!this->m_fa_exist && this->m_pa_exist){
        //碎片分数
        fragment_score = this->m_right_fa_ptr->GetScore();
        //总分是碎片分数加上一致性分数
        this->m_score = this->m_pa_exist_score_weight + this->m_fragment_score_weight*fragment_score;
    }
    //如果FA存在，PA不存在
    else if(this->m_fa_exist && !this->m_pa_exist){
        fragment_score = (this->m_left_pa_ptr->m_left_fa_ptr->GetScore() + this->m_left_pa_ptr->m_right_fa_ptr->GetScore() + this->m_right_fa_ptr->GetScore()) / 3;//碎片分数是3个FA分数的平均值
        fa_consistency_score = float(this->m_fa_info.size()) / 3;
        //总分是碎片分数加上一致性分数
        this->m_score = this->m_fragment_score_weight*fragment_score + this->m_fa_consistency_score_weight*fa_consistency_score;
    }
    //如果FA和PA都存在
    else if(this->m_fa_exist && this->m_pa_exist){
        fragment_score = (this->m_left_pa_ptr->m_left_fa_ptr->GetScore() + this->m_left_pa_ptr->m_right_fa_ptr->GetScore() + this->m_right_fa_ptr->GetScore()) / 3;//碎片分数是4个FA分数的平均值
        fa_consistency_score = float(this->m_fa_info.size()) / 3;
        //总分是碎片分数+一致性分数+pa存在性权重
        this->m_score = this->m_fragment_score_weight*fragment_score + this->m_fa_consistency_score_weight*fa_consistency_score + this->m_pa_exist_score_weight;
    }
}

bool MlclSpecificStructure::operator==(const MlclSpecificStructure &other)
{
//    //如果两个Cl来源于同一个二级，则不可能相同，这是由拼接的过程决定的
//    if(this->m_ms2 == other.m_ms2){
//        return false;
//    }
    //如果两者的FA都不存在，指的是组成PA的2个FA不存在，则比较PA是否相同
    if(!this->m_fa_exist && !other.m_fa_exist){
        //如果两者的PA信息相同，则把this的分数更新为最大值
        if(this->m_pa_info == other.m_pa_info && this->m_fa_info == other.m_fa_info){
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
    //如果一个的PA没有FA，另一个的PA有FA
    else if((!this->m_fa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_fa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的PA信息相同，把this的内容换成新的内容；如果PA相同，则右侧的FA也一定相同
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
    //如果一个只有FA，但是没有PA，另一个既有PA，PA有对应的2个FA，又有FA
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
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other.m_fa_exist) || (!other.m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        MlclSpecificStructure* object_only_fa;
        MlclSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            MlclSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            MlclSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        //判定object_only_fa的3个FA是否能判定成object_only_pa的1个PA
        //用vector存储3个FA
        vector<Fa*> fa_ptr_vector;
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_right_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_fa_ptr);

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
                        //赋值给右侧的那个FA
                        object_only_pa->m_right_fa_ptr = fa_ptr_vector[0];
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
        //如果无法找到2个FA可以拼接成这1个PA，则返回false
        if(!splice_success){
            return false;
        }
    }
    return false;
}

std::shared_ptr<MlclSpecificStructure> MlclSpecificStructure::merge(MlclSpecificStructure *other , bool merge_m_h_m_2h)
{
//    //如果不是用于合并M-H和M-2H
//    if(!merge_m_h_m_2h){
//        //如果两个Cl来源于同一个二级，则不可能相同，这是由拼接的过程决定的
//        if(this->m_ms2 == other->m_ms2){
//            return nullptr;
//        }
//    }
    //如果两者的FA都不存在，指的是组成PA的2个FA不存在，则比较PA是否相同
     if(!this->m_fa_exist && !other->m_fa_exist){
        //如果两者的PA信息相同，则把this的分数更新为最大值
        if(this->m_pa_info == other->m_pa_info && this->m_fa_info == other->m_fa_info){
            shared_ptr<MlclSpecificStructure> new_object = make_shared<MlclSpecificStructure>(*this);
            new_object->m_score = max(this->m_score , other->m_score);
            return new_object;
        }
        else{
            return nullptr;
        }
    }
    //如果两者的PA都不存在，则比较FA是否相同
    else if(!this->m_pa_exist && !other->m_pa_exist){
        //如果两者的FA信息相同，则把this的分数更新为最大值
        if(this->m_fa_info == other->m_fa_info){
            shared_ptr<MlclSpecificStructure> new_object = make_shared<MlclSpecificStructure>(*this);
            new_object->m_score = max(this->m_score , other->m_score);
            return new_object;
        }
        else{
            return nullptr;
        }
    }
    //如果两者的PA，FA都存在
    else if(this->m_fa_exist && other->m_fa_exist && this->m_pa_exist && other->m_pa_exist){
        //如果两者的PA和FA信息相同，则把this的分数更新为最大值
        if((this->m_pa_info == other->m_pa_info) && (this->m_fa_info == other->m_fa_info)){
            shared_ptr<MlclSpecificStructure> new_object = make_shared<MlclSpecificStructure>(*this);
            new_object->m_score = max(this->m_score , other->m_score);
            return new_object;
        }
        else{
            return nullptr;
        }
    }
    //如果一个的PA没有FA，另一个的PA有FA
    else if((!this->m_fa_exist && other->m_fa_exist && other->m_pa_exist) || (!other->m_fa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的PA信息相同，把this的内容换成新的内容；如果PA相同，则右侧的FA也一定相同
        if(this->m_pa_info == other->m_pa_info){
            if(!this->m_fa_exist){
                shared_ptr<MlclSpecificStructure> new_object = make_shared<MlclSpecificStructure>(*other);
                return new_object;
            }
            else if(!other->m_fa_exist){
                shared_ptr<MlclSpecificStructure> new_object = make_shared<MlclSpecificStructure>(*this);
                return new_object;
            }
        }
        else{
            return nullptr;
        }
    }
    //如果一个只有FA，但是没有PA，另一个既有PA，PA有对应的2个FA，又有FA
    else if((!this->m_pa_exist && other->m_fa_exist && other->m_pa_exist) || (!other->m_pa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的FA信息相同，把this的内容换成新的内容，并且把分数换为最大值
        if(this->m_fa_info == other->m_fa_info){
            float max_score = max(this->m_score , other->m_score);
            if(!this->m_pa_exist){
                shared_ptr<MlclSpecificStructure> new_object = make_shared<MlclSpecificStructure>(*other);
                new_object->m_score = max_score;
                return new_object;
            }
            else if(!other->m_pa_exist){
                shared_ptr<MlclSpecificStructure> new_object = make_shared<MlclSpecificStructure>(*this);
                new_object->m_score = max_score;
                return new_object;
            }
        }
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other->m_fa_exist) || (!other->m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        MlclSpecificStructure* object_only_fa;
        shared_ptr<MlclSpecificStructure> object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            object_only_pa = make_shared<MlclSpecificStructure>(*other);
        }
        else{
            object_only_fa = other;
            object_only_pa = make_shared<MlclSpecificStructure>(*this);
        }

        //判定object_only_fa的3个FA是否能判定成object_only_pa的1个PA
        //用vector存储3个FA
        vector<Fa*> fa_ptr_vector;
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_right_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_fa_ptr);

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
                        //赋值给右侧的那个FA
                        object_only_pa->m_right_fa_ptr = fa_ptr_vector[0];
                        object_only_pa->update();//更新分数等信息
                        return object_only_pa;
                    }
                }

            }
            hash_map.insert({(*fa_ptr_vector_itr)->GetChainLength() , *fa_ptr_vector_itr});//把链长作为键，Fa的指针作为值加入到哈希表
        }
        //如果无法找到2个FA可以拼接成这1个PA，则返回false
        if(!splice_success){
            return nullptr;
        }
    }
    return nullptr;
}

void MlclSpecificStructure::update()
{
    //清空FA和PA信息
    set<array<unsigned int,3>>().swap(this->m_pa_info);
    set<array<unsigned int,3>>().swap(this->m_fa_info);
    //更新PA信息
    if(this->m_left_pa_ptr->m_pa_ptr != nullptr){
        this->m_pa_exist = 1;
        this->m_pa_info.insert({this->m_left_pa_ptr->m_pa_ptr->GetChainLength() , this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation() , this->m_left_pa_ptr->m_pa_ptr->GetOxygen()});
    }
    //更新FA信息
    if((this->m_left_pa_ptr->m_left_fa_ptr != nullptr) && (this->m_left_pa_ptr->m_right_fa_ptr != nullptr) && (this->m_right_fa_ptr != nullptr)){
        this->m_fa_exist = 1;
        this->m_fa_info.insert({this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() , this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() , this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()});
        this->m_fa_info.insert({this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength() , this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation() , this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()});
        this->m_fa_info.insert({this->m_right_fa_ptr->GetChainLength() , this->m_right_fa_ptr->GetUnsaturation() , this->m_right_fa_ptr->GetOxygen()});
    }

    //更新得分信息
    this->score();
}

bool MlclSpecificStructure::StrictMerge(MlclSpecificStructure &other)
{
    //如果两者的FA都不存在，指的是组成PA的2个FA不存在，则比较PA是否相同
    if(!this->m_fa_exist && !other.m_fa_exist){
        //如果两者的PA信息相同，则把this的分数更新为最大值
        if(this->m_pa_info == other.m_pa_info && this->m_fa_info == other.m_fa_info){
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
    //如果一个的PA没有FA，另一个的PA有FA
    else if((!this->m_fa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_fa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的PA信息相同，把this的内容换成新的内容；如果PA相同，则右侧的FA也一定相同
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
    //如果一个只有FA，但是没有PA，另一个既有PA，PA有对应的2个FA，又有FA
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
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other.m_fa_exist) || (!other.m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        MlclSpecificStructure* object_only_fa;
        MlclSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            MlclSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            MlclSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        //判定object_only_fa的3个FA是否能判定成object_only_pa的1个PA
        //用vector存储3个FA
        vector<Fa*> fa_ptr_vector;
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_right_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_fa_ptr);

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
                        //赋值给右侧的那个FA
                        object_only_pa->m_right_fa_ptr = fa_ptr_vector[0];
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
        //如果无法找到2个FA可以拼接成这1个PA，则返回false
        if(!splice_success){
            return false;
        }
    }
    return false;
}

bool MlclSpecificStructure::FlexibleMerge(MlclSpecificStructure &other)
{
    //如果两者的FA都不存在，指的是组成PA的2个FA不存在，则比较PA是否相同
    if(!this->m_fa_exist && !other.m_fa_exist){
        //如果两者的PA信息相同，则把this的分数更新为最大值
        if(this->m_pa_info == other.m_pa_info && this->m_fa_info == other.m_fa_info){
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
            this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
            return true;
        }
        else{
            return false;
        }
    }
    //如果两者的PA，FA都存在
    else if(this->m_fa_exist && other.m_fa_exist && this->m_pa_exist && other.m_pa_exist){
        //如果两者FA信息相同，则把this的分数更新为最大值
        if(this->m_fa_info == other.m_fa_info){
            //如果两者的PA信息也相同
            if(this->m_pa_info == other.m_pa_info){
                this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
                return true;
            }
            //如果两者的PA信息不同
            else{
                this->m_score = 1 - (1 - this->m_score)*(1 - other.m_score);
                this->m_flexible_mode_dif_pa_merge.emplace_back(&other);//把跟this的PA不同的MLCL的信息存储到m_flexible_mode_dif_pa_merge中
                return true;
            }
        }
        else{
            return false;
        }
    }
    //如果一个的PA没有FA，另一个的PA有FA
    else if((!this->m_fa_exist && other.m_fa_exist && other.m_pa_exist) || (!other.m_fa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的PA信息相同，把this的内容换成新的内容；如果PA相同，则右侧的FA也一定相同
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
    //如果一个只有FA，但是没有PA，另一个既有PA，PA有对应的2个FA，又有FA
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
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other.m_fa_exist) || (!other.m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        MlclSpecificStructure* object_only_fa;
        MlclSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            MlclSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            MlclSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        //判定object_only_fa的3个FA是否能判定成object_only_pa的1个PA
        //用vector存储3个FA
        vector<Fa*> fa_ptr_vector;
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_left_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_left_pa_ptr->m_right_fa_ptr);
        fa_ptr_vector.push_back(object_only_fa->m_right_fa_ptr);

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
                        //赋值给右侧的那个FA
                        object_only_pa->m_right_fa_ptr = fa_ptr_vector[0];
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
        //如果无法找到2个FA可以拼接成这1个PA，则返回false
        if(!splice_success){
            return false;
        }
    }
    return false;
}

QString MlclSpecificStructure::ShowInfo()
{
    QString message;
    //如果PA和FA都存在
    if(this->m_fa_exist && this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_pa_ptr->GetAdditiveForm()) + ");";
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_left_fa_ptr->GetAdditiveForm())+ ");";
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_right_fa_ptr->GetAdditiveForm())+ ");";
        QString Fa_3_message = QString::number(this->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_right_fa_ptr->GetAdditiveForm()) + ");";
        message =  message + "This is MLCl : have both PA and FA :" + "PA1:" + Pa1_message + "  FA1:" +  Fa_1_message  + "  FA2 :" + Fa_2_message + "   FA3 :" + Fa_3_message + "score :" + QString::number(this->m_score);
    }
    //如果只有FA存在
    else if(this->m_fa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_left_fa_ptr->GetAdditiveForm())+ ");";
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_right_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_3_message = QString::number(this->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_right_fa_ptr->GetAdditiveForm()) + ");";

        message =  message + "This is MLCl : have only FA :" "FA1:" +  Fa_1_message  + "    FA2 :" + Fa_2_message + "   FA3 :" + Fa_3_message + "   score :" + QString::number(this->m_score);
    }
    else if(this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_pa_ptr->GetAdditiveForm()) + ");";
        QString Fa_3_message = QString::number(this->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_right_fa_ptr->GetAdditiveForm()) + ");";
        message =  message + "This is MLCl : have only PA :" + "PA1:" + Pa1_message + " FA3" + Fa_3_message + " score :" + QString::number(this->m_score);
    }

    return message;
}

QString MlclSpecificStructure::ShowSimpleInfo()
{
    QString message;
    //如果PA和FA都存在
    if(this->m_fa_exist && this->m_pa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen());
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen());
        QString Fa_3_message = QString::number(this->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_fa_ptr->GetOxygen());
        message =  message +  Fa_1_message + "/" + Fa_2_message + "/" + Fa_3_message;
    }
    //如果只有FA存在
    else if(this->m_fa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen());
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen());
        QString Fa_3_message = QString::number(this->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_fa_ptr->GetOxygen());
        message =  message +  Fa_1_message + "/" + Fa_2_message + "/" + Fa_3_message;
    }
    //如果只有PA存在
    else if(this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen());
        QString Fa_3_message = QString::number(this->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_right_fa_ptr->GetOxygen());
        message =  message + Pa1_message + "/" + Fa_3_message;
    }

    return message;
}
