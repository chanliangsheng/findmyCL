#include <DlclSpecificStructure.h>

float DlclSpecificStructure::m_fragment_score_weight = 1;
float DlclSpecificStructure::m_fa_consistency_score_weight = 1;
float DlclSpecificStructure::m_pa_exist_score_weight = 1;
float DlclSpecificStructure::m_fa_intensity_variance_score_weight = 1;

using namespace std;

DlclSpecificStructure::DlclSpecificStructure()
{
    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;

    this->m_pa_exist = 0;
    this->m_fa_exist = 0;
    this->m_left_pa_ptr = nullptr;
    this->m_left_pa_ptr->m_left_fa_ptr = nullptr;
    this->m_left_pa_ptr->m_right_fa_ptr = nullptr;
    this->m_ms2 = nullptr;

    this->m_score = 0;
}

DlclSpecificStructure::DlclSpecificStructure(Fa *fa_1_ptr, Fa *fa_2_ptr, Ms2 *ms2_ptr)
{
    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;

    //把3个FA的chain，unsaturation，oxygen的信息存储到m_fa_info中
    this->m_fa_info.insert({fa_1_ptr->GetChainLength() , fa_1_ptr->GetUnsaturation() , fa_1_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_2_ptr->GetChainLength() , fa_2_ptr->GetUnsaturation() , fa_2_ptr->GetOxygen()});

    this->m_ms2 = ms2_ptr;
    this->m_left_pa_ptr->m_left_fa_ptr = fa_1_ptr;
    this->m_left_pa_ptr->m_right_fa_ptr = fa_2_ptr;
    this->m_left_pa_ptr->m_pa_ptr = nullptr;

    //pa或者fa是否存在
    this->m_pa_exist = 0;
    this->m_fa_exist = 1;

    //进行打分
    this->score();
    //计算总强度
    this->CalculateTotalIntensity();
}

DlclSpecificStructure::DlclSpecificStructure(Pa *pa_1_ptr, Fa *fa_1_ptr, Fa *fa_2_ptr, Ms2 *ms2_ptr)
{
    this->m_ms2 = ms2_ptr;

    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;

    this->m_left_pa_ptr->m_pa_ptr = pa_1_ptr;
    this->m_left_pa_ptr->m_left_fa_ptr = fa_1_ptr;
    this->m_left_pa_ptr->m_right_fa_ptr = fa_2_ptr;

    //pa或者fa是否存在
    this->m_pa_exist = 1;
    this->m_fa_exist = 1;

    //把1个PA，3个FA的chain，unsaturation，oxygen的信息存储到m_pa_info中
    this->m_pa_info.insert({pa_1_ptr->GetChainLength() , pa_1_ptr->GetUnsaturation() , pa_1_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_1_ptr->GetChainLength() , fa_1_ptr->GetUnsaturation() , fa_1_ptr->GetOxygen()});
    this->m_fa_info.insert({fa_2_ptr->GetChainLength() , fa_2_ptr->GetUnsaturation() , fa_2_ptr->GetOxygen()});

    //进行打分
    this->score();
    //计算总强度
    this->CalculateTotalIntensity();
}

DlclSpecificStructure::DlclSpecificStructure(Pa *pa_1_ptr, Ms2 *ms2_ptr)
{
    this->m_ms2 = ms2_ptr;
    //初始化两个PA节点，方式this的两个Pa节点为空
    PaNode* left_pa_node = new PaNode;
    this->m_left_pa_ptr = left_pa_node;
    this->m_left_pa_ptr->m_left_fa_ptr = nullptr;
    this->m_left_pa_ptr->m_right_fa_ptr = nullptr;
    this->m_left_pa_ptr->m_pa_ptr = pa_1_ptr;

    //pa或者fa是否存在
    this->m_pa_exist = 1;
    this->m_fa_exist = 0;

    //把1个PA的chain，unsaturation，oxygen的信息存储到m_pa_info中
    this->m_pa_info.insert({pa_1_ptr->GetChainLength() , pa_1_ptr->GetUnsaturation() , pa_1_ptr->GetOxygen()});

    //进行打分
    this->score();
    //计算总强度
    this->CalculateTotalIntensity();
}

void DlclSpecificStructure::score()
{
    float fragment_score = 0;
    float fa_consistency_score = 0;

    //如果这个组合的FA不存在，但是PA存在
    if(!this->m_fa_exist && this->m_pa_exist){
        //总分是PA存在的权重
        this->m_score = this->m_pa_exist_score_weight;
    }
    //如果FA存在，PA不存在
    else if(this->m_fa_exist && !this->m_pa_exist){
        //把2个FA的指针加入到哈希表中，这样可以把一样的FA合并，合并的结果是key：FA，value：FA出现的次数
        unordered_map<Fa*,float> fa_hashmap;
        fa_hashmap.insert({this->m_left_pa_ptr->m_left_fa_ptr , 1});

        auto find_itr = fa_hashmap.find(this->m_left_pa_ptr->m_right_fa_ptr);
        if(find_itr == fa_hashmap.end()){
            fa_hashmap.insert({this->m_left_pa_ptr->m_right_fa_ptr,1});
        }
        else{
            find_itr->second++;
        }


        //得到碎片分数
        for(auto itr = fa_hashmap.begin() ; itr != fa_hashmap.end() ; itr++){
            fragment_score += (*itr).first->GetScore();
        }

        //碎片分数是总分的平均值
        fragment_score = fragment_score / 2;

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
        float fa_intensity_variance_score = 1 - sqrt(sum/2);//分数为1-标准差，方差越大，分数越小

        //fa一致性分数
        fa_consistency_score = float(this->m_fa_info.size()) / 2;
        //总分是碎片分数加上一致性分数
        this->m_score = this->m_fragment_score_weight*fragment_score + this->m_fa_consistency_score_weight*fa_consistency_score + this->m_fa_intensity_variance_score_weight*fa_intensity_variance_score;
    }
    //如果FA和PA都存在
    else if(this->m_fa_exist && this->m_pa_exist){
        //把2个FA的指针加入到哈希表中，这样可以把一样的FA合并，合并的结果是key：FA，value：FA出现的次数
        unordered_map<Fa*,float> fa_hashmap;
        fa_hashmap.insert({this->m_left_pa_ptr->m_left_fa_ptr , 1});

        auto find_itr = fa_hashmap.find(this->m_left_pa_ptr->m_right_fa_ptr);
        if(find_itr == fa_hashmap.end()){
            fa_hashmap.insert({this->m_left_pa_ptr->m_right_fa_ptr,1});
        }
        else{
            find_itr->second++;
        }


        //得到碎片分数
        for(auto itr = fa_hashmap.begin() ; itr != fa_hashmap.end() ; itr++){
            fragment_score += (*itr).first->GetScore();
        }

        //碎片分数是总分的平均值
        fragment_score = fragment_score / 2;

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
        float fa_intensity_variance_score = 1 - sqrt(sum/2);//分数为1-标准差，方差越大，分数越小

        //fa一致性分数
        fa_consistency_score = float(this->m_fa_info.size()) / 2;
        //总分是碎片分数加上一致性分数
        this->m_score = this->m_fragment_score_weight*fragment_score + this->m_fa_consistency_score_weight*fa_consistency_score + this->m_pa_exist_score_weight + this->m_fa_intensity_variance_score_weight*fa_intensity_variance_score;
    }
}

float DlclSpecificStructure::GetTotalIntensity()
{
    return this->m_total_intensity;
}

void DlclSpecificStructure::CalculateTotalIntensity()
{
    //如果这个组合的FA不存在，但是PA存在
    if(!this->m_fa_exist && this->m_pa_exist){
        this->m_total_intensity += this->m_left_pa_ptr->m_pa_ptr->GetIntensity();
    }
    //如果FA存在，PA不存在
    else if(this->m_fa_exist && !this->m_pa_exist){
        //把2个FA指针加入set中，因为可以把一样的FA合并
        set<Fa*> fa_set;
        fa_set.insert(this->m_left_pa_ptr->m_left_fa_ptr);
        fa_set.insert(this->m_left_pa_ptr->m_right_fa_ptr);
        for(auto itr = fa_set.begin();itr != fa_set.end();itr++){
            this->m_total_intensity += (*itr)->GetIntensity();
        }
    }
    //如果FA和PA都存在
    else if(this->m_fa_exist && this->m_pa_exist){
        this->m_total_intensity += this->m_left_pa_ptr->m_pa_ptr->GetIntensity();
        //把2个FA指针加入set中，因为可以把一样的FA合并
        set<Fa*> fa_set;
        fa_set.insert(this->m_left_pa_ptr->m_left_fa_ptr);
        fa_set.insert(this->m_left_pa_ptr->m_right_fa_ptr);
        for(auto itr = fa_set.begin();itr != fa_set.end();itr++){
            this->m_total_intensity += (*itr)->GetIntensity();
        }

    }
}

float DlclSpecificStructure::GetMs2TotalIntensity()
{
    return this->m_ms2->GetTotalIntensity();
}

bool DlclSpecificStructure::operator==(const DlclSpecificStructure &other)
{
//    //如果两个Cl来源于同一个二级，则不可能相同，这是由拼接的过程决定的
//    if(this->m_ms2 == other.m_ms2){
//        return false;
//    }
    //如果两者的FA都不存在，指的是组成PA的2个FA不存在，则比较PA是否相同
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
    //如果一个只有FA，另一个既有PA，PA有对应的2个FA
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
        DlclSpecificStructure* object_only_fa;
        DlclSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            DlclSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            DlclSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        bool splice_success = 0;
        if((object_only_pa->m_left_pa_ptr->m_pa_ptr->GetChainLength() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) &&
           (object_only_pa->m_left_pa_ptr->m_pa_ptr->GetUnsaturation() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) &&
           (object_only_pa->m_left_pa_ptr->m_pa_ptr->GetOxygen() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetOxygen() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetOxygen())){
            splice_success = 1;
            object_only_pa->m_left_pa_ptr->m_left_fa_ptr = object_only_fa->m_left_pa_ptr->m_left_fa_ptr;
            object_only_pa->m_left_pa_ptr->m_right_fa_ptr = object_only_fa->m_left_pa_ptr->m_right_fa_ptr;
            if(this == object_only_fa){
              *this = *object_only_pa;
            }
            this->update();//更新分数等信息
            return true;
        }
        splice_success = 1;
        if(!splice_success){
            return false;
        }
    }
    return false;
}

std::shared_ptr<DlclSpecificStructure> DlclSpecificStructure::merge(DlclSpecificStructure *other , bool merge_m_h_m_2h)
{
    //如果不是用于合并M-H和M-2H
    if(!merge_m_h_m_2h){
        //如果两个Cl来源于同一个二级，则不可能相同，这是由拼接的过程决定的
        if(this->m_ms2 == other->m_ms2){
            return nullptr;
        }
    }
    //如果两者的FA都不存在，指的是组成PA的2个FA不存在，则比较PA是否相同
    else if(!this->m_fa_exist && !other->m_fa_exist){
        //如果两者的PA信息相同，则把this的分数更新为最大值
        if(this->m_pa_info == other->m_pa_info){
            shared_ptr<DlclSpecificStructure> new_object = make_shared<DlclSpecificStructure>(*this);
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
            shared_ptr<DlclSpecificStructure> new_object = make_shared<DlclSpecificStructure >(*this);
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
            shared_ptr<DlclSpecificStructure> new_object = make_shared<DlclSpecificStructure>(*this);
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
                shared_ptr<DlclSpecificStructure> new_object = make_shared<DlclSpecificStructure>(*other);
                return new_object;
            }
            else if(!other->m_fa_exist){
                shared_ptr<DlclSpecificStructure> new_object = make_shared<DlclSpecificStructure>(*this);
                return new_object;
            }
        }
        else{
            return nullptr;
        }
    }
    //如果一个只有FA，另一个既有PA，PA有对应的2个FA
    else if((!this->m_pa_exist && other->m_fa_exist && other->m_pa_exist) || (!other->m_pa_exist && this->m_fa_exist && this->m_pa_exist)){
        //如果两者的FA信息相同，把this的内容换成新的内容，并且把分数换为最大值
        if(this->m_fa_info == other->m_fa_info){
            float max_score = max(this->m_score , other->m_score);
            if(!this->m_pa_exist){
                shared_ptr<DlclSpecificStructure> new_object = make_shared<DlclSpecificStructure>(*other);
                new_object->m_score = max_score;
                return new_object;
            }
            else if(!other->m_pa_exist){
                shared_ptr<DlclSpecificStructure> new_object = make_shared<DlclSpecificStructure>(*this);
                new_object->m_score = max_score;
                return new_object;
            }
        }
    }
    //如果一个只有FA，另一个只有PA
    else if((!this->m_pa_exist && !other->m_fa_exist) || (!other->m_pa_exist && !this->m_pa_exist)){
        //抽象出两个指针，一个指向只有FA的，另一个指向只有PA的；交换指针
        DlclSpecificStructure* object_only_fa;
        shared_ptr<DlclSpecificStructure> object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            object_only_pa = make_shared<DlclSpecificStructure>(*other);
        }
        else{
            object_only_fa = other;
            object_only_pa = make_shared<DlclSpecificStructure>(*this);
        }

        bool splice_success = 0;
        if((object_only_pa->m_left_pa_ptr->m_pa_ptr->GetChainLength() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) &&
           (object_only_pa->m_left_pa_ptr->m_pa_ptr->GetUnsaturation() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) &&
           (object_only_pa->m_left_pa_ptr->m_pa_ptr->GetOxygen() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetOxygen() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetOxygen())){
            splice_success = 1;
            object_only_pa->m_left_pa_ptr->m_left_fa_ptr = object_only_fa->m_left_pa_ptr->m_left_fa_ptr;
            object_only_pa->m_left_pa_ptr->m_right_fa_ptr = object_only_fa->m_left_pa_ptr->m_right_fa_ptr;
            object_only_pa->update();//更新分数等信息
            return object_only_pa;
        }
        splice_success = 1;
        if(!splice_success){
            return nullptr;
        }
    }
    return nullptr;
}

void DlclSpecificStructure::update()
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
    if((this->m_left_pa_ptr->m_left_fa_ptr != nullptr) && (this->m_left_pa_ptr->m_right_fa_ptr != nullptr)){
        this->m_fa_exist = 1;
        this->m_fa_info.insert({this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() , this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() , this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()});
        this->m_fa_info.insert({this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength() , this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation() , this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()});
    }
    //更新得分信息
    this->score();
}

bool DlclSpecificStructure::StrictMerge(DlclSpecificStructure &other)
{
    //如果两者的FA都不存在，指的是组成PA的2个FA不存在，则比较PA是否相同
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
    //如果一个只有FA，另一个既有PA，PA有对应的2个FA
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
        DlclSpecificStructure* object_only_fa;
        DlclSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            DlclSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            DlclSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        bool splice_success = 0;
        if((object_only_pa->m_left_pa_ptr->m_pa_ptr->GetChainLength() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) &&
           (object_only_pa->m_left_pa_ptr->m_pa_ptr->GetUnsaturation() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) &&
           (object_only_pa->m_left_pa_ptr->m_pa_ptr->GetOxygen() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetOxygen() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetOxygen())){
            splice_success = 1;
            object_only_pa->m_left_pa_ptr->m_left_fa_ptr = object_only_fa->m_left_pa_ptr->m_left_fa_ptr;
            object_only_pa->m_left_pa_ptr->m_right_fa_ptr = object_only_fa->m_left_pa_ptr->m_right_fa_ptr;
            if(this == object_only_fa){
              *this = *object_only_pa;
            }
            this->update();//更新分数等信息
            return true;
        }
        splice_success = 1;
        if(!splice_success){
            return false;
        }
    }
    return false;
}

bool DlclSpecificStructure::FlexibleMerge(DlclSpecificStructure &other)
{
    //如果两者的FA都不存在，指的是组成PA的2个FA不存在，则比较PA是否相同
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
    //如果一个只有FA，另一个既有PA，PA有对应的2个FA
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
        DlclSpecificStructure* object_only_fa;
        DlclSpecificStructure* object_only_pa;
        if(!this->m_pa_exist){
            object_only_fa = this;
            DlclSpecificStructure tem_var = other;
            object_only_pa = &tem_var;
        }
        else{
            DlclSpecificStructure tem_var = other;
            object_only_fa = &tem_var;
            object_only_pa = this;
        }

        bool splice_success = 0;
        if((object_only_pa->m_left_pa_ptr->m_pa_ptr->GetChainLength() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) &&
           (object_only_pa->m_left_pa_ptr->m_pa_ptr->GetUnsaturation() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) &&
           (object_only_pa->m_left_pa_ptr->m_pa_ptr->GetOxygen() == object_only_fa->m_left_pa_ptr->m_left_fa_ptr->GetOxygen() + object_only_fa->m_left_pa_ptr->m_right_fa_ptr->GetOxygen())){
            splice_success = 1;
            object_only_pa->m_left_pa_ptr->m_left_fa_ptr = object_only_fa->m_left_pa_ptr->m_left_fa_ptr;
            object_only_pa->m_left_pa_ptr->m_right_fa_ptr = object_only_fa->m_left_pa_ptr->m_right_fa_ptr;
            if(this == object_only_fa){
              *this = *object_only_pa;
            }
            this->update();//更新分数等信息
            return true;
        }
        splice_success = 1;
        if(!splice_success){
            return false;
        }
    }
    return false;
}

QString DlclSpecificStructure::ShowInfo()
{
    QString message;
    //如果PA和FA都存在
    if(this->m_fa_exist && this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_pa_ptr->GetAdditiveForm()) + ");";
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_left_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_right_fa_ptr->GetAdditiveForm()) + ");";
        message =  message + "This is DLCl : have both PA and FA :" + "PA1:" + Pa1_message + "  FA1:" +  Fa_1_message  + "  FA2 :" + Fa_2_message + "   score :" + QString::number(this->m_score);
    }
    //如果只有FA存在
    else if(this->m_fa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_left_fa_ptr->GetAdditiveForm()) + ");";
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_right_fa_ptr->GetAdditiveForm()) + ");";
        message =  message + "This is DLCl : have only FA :" "FA1:" +  Fa_1_message  + "    FA2 :" + Fa_2_message + "   score :" + QString::number(this->m_score);
    }
    else if(this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen()) + "(mz:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetMz()) +"intensity:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetIntensity()) + "score:" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetScore()) + "add_form:" + QString::fromStdString(this->m_left_pa_ptr->m_pa_ptr->GetAdditiveForm()) + ");";
        message =  message + "This is DLCl : have only PA :" + "PA1:" + Pa1_message + " score :" + QString::number(this->m_score);
    }
    return message;
}

QString DlclSpecificStructure::ShowSimpleInfo()
{
    QString message;
    //如果PA和FA都存在
    if(this->m_fa_exist && this->m_pa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen());
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen());
        message =  message +  Fa_1_message + "/" + Fa_2_message;
    }
    //如果只有FA存在
    else if(this->m_fa_exist){
        QString Fa_1_message = QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen());
        QString Fa_2_message = QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation()) + ":" +QString::number(this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen());
        message =  message +  Fa_1_message + "/" + Fa_2_message;
    }
    //如果只有PA存在
    else if(this->m_pa_exist){
        QString Pa1_message = QString::number(this->m_left_pa_ptr->m_pa_ptr->GetChainLength()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation()) + ":" + QString::number(this->m_left_pa_ptr->m_pa_ptr->GetOxygen());
        message =  message + Pa1_message;
    }

    return message;
}

std::vector<int> DlclSpecificStructure::GetTotalInfo()
{
    return {this->GetTotalChainLength() , this->GetTotalUnsaturation() , this->GetTotalOxygen()};
}

int DlclSpecificStructure::GetTotalChainLength()
{
    if(this->m_pa_exist){
        return this->m_left_pa_ptr->m_pa_ptr->GetChainLength();
    }
    else{
        return this->m_left_pa_ptr->m_left_fa_ptr->GetChainLength() + this->m_left_pa_ptr->m_right_fa_ptr->GetChainLength();
    }
}

int DlclSpecificStructure::GetTotalUnsaturation()
{
    if(this->m_pa_exist){
        return this->m_left_pa_ptr->m_pa_ptr->GetUnsaturation();
    }
    else{
        return this->m_left_pa_ptr->m_left_fa_ptr->GetUnsaturation() + this->m_left_pa_ptr->m_right_fa_ptr->GetUnsaturation();
    }
}

int DlclSpecificStructure::GetTotalOxygen()
{
    if(this->m_pa_exist){
        return this->m_left_pa_ptr->m_pa_ptr->GetOxygen();
    }
    else{
        return this->m_left_pa_ptr->m_left_fa_ptr->GetOxygen() + this->m_left_pa_ptr->m_right_fa_ptr->GetOxygen();
    }
}
