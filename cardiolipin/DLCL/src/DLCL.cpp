#include <DLCL.h>

float DlclSpecificStructure::m_fragment_score_weight = 1;
float DlclSpecificStructure::m_fa_consistency_score_weight = 1;
float DlclSpecificStructure::m_pa_exist_score_weight = 1;


using namespace std;
Dlcl::Dlcl()
{

}

Dlcl::Dlcl(Cardiolipin &superclass)
{
    if(superclass.GetMs1Ptr() != nullptr){
        this->SetMs1Ptr(superclass.GetMs1Ptr());
    }
    else{
        this->SetMs1Ptr(nullptr);
    }
    if(superclass.GetMs1DatabaseRecordPtr() != nullptr){
        this->SetMs1DatabaseRecordPtr(superclass.GetMs1DatabaseRecordPtr());
    }
    else{
        this->SetMs1DatabaseRecordPtr(nullptr);
    }

    this->SetMs1MatchingScosre(superclass.GetMs1MatchingScore());
    this->SetMs2VectorPtr(superclass.GetMs2VectorPtr());
}

void Dlcl::splice()
{
    //如果是空对象，则不进行操作
    if(this->CheckEmptyObject()){
        return;
    }

    //对不同的情况进行拼接操作
    std::vector<Ms2*> ms2_ptr_vector = this->GetMs2VectorPtr();
    //对不同的情况进行拼接操作
    for(auto ms2_itr = ms2_ptr_vector.begin() ; ms2_itr != ms2_ptr_vector.end();ms2_itr++){
        //如果这个二级没有PA，尝试用2个FA拼接
        if((*ms2_itr)->GetPaCount() == 0){
            this->TwoFaSpliceDlcl(*ms2_itr);
        }
        //如果这个二级即有FA，又有PA
        else if((*ms2_itr)->GetPaCount() != 0 && (*ms2_itr)->GetFaCount() != 0)
        {
            this->OnePaTwoFaSpliceDlcl(*ms2_itr);
        }
    }

    //如果没有拼接成功，则清空对象
    if(this->m_dlcl_specific_structure_vector.size() == 0){
        this->EmptyObject();
        return;
    }

    //合并不同二级之间的拼接结果
    this->MergeSplice();
}

void Dlcl::MergeSplice()
{
//    //如果只有一个二级，则不需要合并，因为一个二级中的所有拼接结果都是不重复的
//    if(this->m_Ms2_vector_ptr.size() == 1){
//        return;
//    }
    list<DlclSpecificStructure> non_repeating_only_pa_list;//存储只有PA和FA的对象，用list来存储
    list<DlclSpecificStructure> non_repeating_only_fa_list;//存储只有PA和FA的对象，用list来存储
    list<DlclSpecificStructure> non_repeating_have_both_pa_fa_list;//存储既有PA和又有FA的对象，用list来存储

    list<DlclSpecificStructure> compare_to_have_both_pa_fa_list;

    //先把只有PA的跟只有PA的合并，只有FA的跟只有FA的合并，两者都有的更两者都有的合并；分别存储到object_only_pa_only_fa_list和object_with_pa_fa_list中
    for(auto first_itr = this->m_dlcl_specific_structure_vector.begin() ; first_itr != this->m_dlcl_specific_structure_vector.end();first_itr++){
        if(!first_itr->m_fa_exist){
            for(auto second_itr = next(first_itr) ; second_itr != this->m_dlcl_specific_structure_vector.end();){
                if(!second_itr->m_fa_exist){
                    //如果这两个Cl相同，进行这个操作的时候，first_itr对应的Cl内部会发生改变
                    if(*first_itr == *second_itr){
                        //如果这两个Cl相同，那么从m_cl_specific_structure_vector中去除第二个Cl
                        second_itr =  this->m_dlcl_specific_structure_vector.erase(second_itr);
                    }
                    else{
                        second_itr++;
                    }
                }
                else{
                    second_itr++;
                }
            }
            non_repeating_only_pa_list.push_back(*first_itr);//加入到object_only_pa_only_fa_list中
        }
        else if(!first_itr->m_pa_exist){
            for(auto second_itr = next(first_itr) ; second_itr != this->m_dlcl_specific_structure_vector.end();){
                if(!second_itr->m_pa_exist){
                    //如果这两个Cl相同，进行这个操作的时候，first_itr对应的Cl内部会发生改变
                    if(*first_itr == *second_itr){
                        second_itr =  this->m_dlcl_specific_structure_vector.erase(second_itr);
                    }
                    else{
                        second_itr++;
                    }
                }
                else{
                    second_itr++;
                }
            }
            non_repeating_only_fa_list.push_back(*first_itr);//加入到object_only_pa_only_fa_list中
        }
        else if(first_itr->m_fa_exist && first_itr->m_pa_exist){
            for(auto second_itr = next(first_itr) ; second_itr != this->m_dlcl_specific_structure_vector.end();){
                if(second_itr->m_fa_exist && second_itr->m_pa_exist){
                    //如果这两个Cl相同，进行这个操作的时候，first_itr对应的Cl内部会发生改变
                    if(*first_itr == *second_itr){
                        second_itr =  this->m_dlcl_specific_structure_vector.erase(second_itr);
                    }
                    else{
                        second_itr++;
                    }
                }
                else{
                    second_itr++;
                }
            }
            non_repeating_have_both_pa_fa_list.push_back(*first_itr);//加入到object_only_pa_only_fa_list中
        }
    }

    //进行仅有PA的和仅有FA的CL的重复合并，最终结果追加到only_pa_fa_merge_to_both_pafa_list和only_pa_fa_not_merge_to_both_pafa_list中
    set<DlclSpecificStructure*> fa_can_be_spliced_with_pa;//判断哪些FA被合并了
    //合并仅有PA的和仅有FA的,non_repeating_only_pa_list和non_repeating_only_fa_list进行合并
    for(auto first_itr = non_repeating_only_pa_list.begin() ; first_itr != non_repeating_only_pa_list.end();first_itr++){
        bool first_itr_success = 0;//判定这个Cl是否能在后面找到成功拼接上的
        for(auto second_itr = non_repeating_only_fa_list.begin() ; second_itr != non_repeating_only_fa_list.end();second_itr++){
            //如果这两个Cl相同，返回新的Cl的共享指针，如果不同，返回nullptr
            shared_ptr<DlclSpecificStructure> merge_result = first_itr->merge(&*second_itr);
            if(merge_result != nullptr){
                first_itr_success = 1;//设置成成功
                compare_to_have_both_pa_fa_list.push_back(*merge_result);//追加到结果中
                fa_can_be_spliced_with_pa.insert(&*second_itr);
            }
            else{

            }
        }
        //如果没有找到成功拼接的情况，则把这个仅有PA或者仅有FA的Cl加入到结果中
        if(!first_itr_success){
            compare_to_have_both_pa_fa_list.push_back(*first_itr);//追加到结果中
        }
    }
    //判断non_repeating_only_fa_list中的哪些FA被使用了
    for(auto itr = non_repeating_only_fa_list.begin() ; itr != non_repeating_only_fa_list.end() ;){
        if(fa_can_be_spliced_with_pa.find(&*itr) != fa_can_be_spliced_with_pa.end()){
            itr = non_repeating_only_fa_list.erase(itr);
        }
        else{
            itr++;
        }
    }
    //PA和FA合并产生了两个新结果，compare_to_have_both_pa_fa_list和non_repeating_only_fa_list的减少，把两者连起来，里面的元素都是互相之间无法拼接的
    compare_to_have_both_pa_fa_list.splice(compare_to_have_both_pa_fa_list.end() , non_repeating_only_fa_list);


    list<DlclSpecificStructure>().swap(this->m_dlcl_specific_structure_vector);//清空m_cl_specific_structure_vector的内容
    //最终合并
    for(auto first_itr = non_repeating_have_both_pa_fa_list.begin() ; first_itr!=non_repeating_have_both_pa_fa_list.end() ; first_itr++){
        for(auto second_itr = compare_to_have_both_pa_fa_list.begin() ; second_itr != compare_to_have_both_pa_fa_list.end();){
            //如果这两个Cl相同，进行这个操作的时候，first_itr对应的Cl内部会发生改变
            if(*first_itr == *second_itr){
                second_itr =  compare_to_have_both_pa_fa_list.erase(second_itr);
            }
            else{
                second_itr++;
            }
        }
        this->m_dlcl_specific_structure_vector.push_back(*first_itr);
    }
    this->m_dlcl_specific_structure_vector.splice(this->m_dlcl_specific_structure_vector.end() , compare_to_have_both_pa_fa_list);
}

void Dlcl::TwoFaSpliceDlcl(Ms2 *ms2_ptr)
{
    vector<Fa>* fa_vector_ptr = ms2_ptr->GetFaVectorPtr();

    //哈系法解两数之和
    multimap<unsigned int,Fa*> hash_map;//以链长为键，对应的Fa为值构建多重哈希表
    for(auto fa_vector_itr = fa_vector_ptr->begin();fa_vector_itr != fa_vector_ptr->end();fa_vector_itr++){
        hash_map.insert({fa_vector_itr->GetChainLength() , &(*fa_vector_itr)});//插入数据
        auto search_pair_itr = hash_map.equal_range(this->GetChainLength() - fa_vector_itr->GetChainLength());//寻找剩余链长的键的范围
        //如果没有找到
        if(search_pair_itr.first == search_pair_itr.second){

        }
        else{
            //遍历搜索到的范围
            for(auto match_itr = search_pair_itr.first ; match_itr != search_pair_itr.second ; match_itr++){
                //如果链长相加起来等于目标链长的情况下，如果不饱和度和氧个数加起来都等于目标的这些属性，加入到结果中
                if(((fa_vector_itr->GetUnsaturation() + match_itr->second->GetUnsaturation()) == this->GetUnsaturation()) &&
                   ((fa_vector_itr->GetOxygen() + match_itr->second->GetOxygen()) == this->GetOxygen())){
                    this->m_dlcl_specific_structure_vector.emplace_back(DlclSpecificStructure(&*fa_vector_itr , match_itr->second , ms2_ptr));//加入到结果中
                }
            }
        }
    }

}

void Dlcl::OnePaTwoFaSpliceDlcl(Ms2 *ms2_ptr)
{
    //用于存储哪2个FA拼接成了1个PA
    set<set<Fa*>> which_two_fa_splice_one_pa;

    //以FA的链长对fa_vector进行排序
    ms2_ptr->SortFaByChainLength();
    vector<Fa>* fa_vector_ptr = ms2_ptr->GetFaVectorPtr();
    vector<Pa>* pa_vector_ptr = ms2_ptr->GetPaVectorPtr();

    multimap<unsigned int,Fa*> fa_hash_map;//以链长为键，对应的fa为值构建多重哈希表
    //把FA的信息加入到表中
    for(auto fa_vector_itr = fa_vector_ptr->begin() ; fa_vector_itr != fa_vector_ptr->end();fa_vector_itr++){
        fa_hash_map.insert({fa_vector_itr->GetChainLength() , &*fa_vector_itr});
    }

    //寻找哪些PA可以构成这个DLCL，把可以组成的PA的2个FA加入到set中
    for(auto pa_vector_itr = pa_vector_ptr->begin();pa_vector_itr != pa_vector_ptr->end();pa_vector_itr++){
        if((pa_vector_itr->GetChainLength() == this->GetChainLength()) &&
           (pa_vector_itr->GetUnsaturation() == this->GetUnsaturation()) &&
           (pa_vector_itr->GetOxygen() == this->GetOxygen())){
            bool splice_pa_success = 0;
            //用2个FA拼接成这个PA
            for(auto fa_vector_itr = fa_vector_ptr->begin() ; fa_vector_itr != fa_vector_ptr->end();fa_vector_itr++){
                auto left_fa_search_itr = fa_hash_map.equal_range(pa_vector_itr->GetChainLength() - fa_vector_itr->GetChainLength());
                //如果没有找到
                if(left_fa_search_itr.first == left_fa_search_itr.second){

                }
                else{
                    for(auto itr = left_fa_search_itr.first ; itr!= left_fa_search_itr.second ; itr++){
                        if(((pa_vector_itr->GetUnsaturation() == fa_vector_itr->GetUnsaturation()) + itr->second->GetUnsaturation()) &&
                           ((pa_vector_itr->GetOxygen() == fa_vector_itr->GetOxygen()) + itr->second->GetOxygen())){
                            splice_pa_success = 1;
                            //2FA组成了这个PA
                            this->m_dlcl_specific_structure_vector.emplace_back(DlclSpecificStructure(&*pa_vector_itr , itr->second , &*fa_vector_itr , ms2_ptr));
                            which_two_fa_splice_one_pa.insert({itr->second , &*fa_vector_itr});
                        }
                    }
                }
            }
            //如果不成功，则加入1PA到结果中
            if(!splice_pa_success){
                this->m_dlcl_specific_structure_vector.emplace_back(DlclSpecificStructure(&*pa_vector_itr , ms2_ptr));
            }

        }
    }

     //用2个FA去拼接这个DLCL，如果这2个FA可以拼接成PA，则跳过这种组合
    for(auto fa_vector_itr = fa_vector_ptr->begin() ; fa_vector_itr != fa_vector_ptr->end();fa_vector_itr++){
        auto left_fa_search_itr = fa_hash_map.equal_range(this->GetChainLength() - fa_vector_itr->GetChainLength());
        //如果没有找到
        if(left_fa_search_itr.first == left_fa_search_itr.second){

        }
        else{
            for(auto itr = left_fa_search_itr.first ; itr!= left_fa_search_itr.second ; itr++){
                if((itr->second->GetUnsaturation() + fa_vector_itr->GetUnsaturation() == this->GetUnsaturation()) &&
                   (itr->second->GetOxygen() + fa_vector_itr->GetOxygen() == this->GetOxygen())){
                    //如果这2个FA的组合可以组成PA，则不考虑这种组合
                    if(which_two_fa_splice_one_pa.find({itr->second , &*fa_vector_itr}) == which_two_fa_splice_one_pa.end()){
                        this->m_dlcl_specific_structure_vector.emplace_back(DlclSpecificStructure(&*fa_vector_itr , itr->second , ms2_ptr));
                    }
                }
            }
        }
    }
}

void Dlcl::EmptyObject()
{
    if(this->GetMs1Ptr() != nullptr){
        this->SetMs1Ptr(nullptr);
    }
    if(this->GetMs1DatabaseRecordPtr() != nullptr){
        this->SetMs1DatabaseRecordPtr(nullptr);
    }

    this->SetMs1MatchingScosre(0);
    this->ClearMs2Info();//清空vector中的元素

    if(this->m_dlcl_specific_structure_vector.size() != 0){
        //去除每个DLCL详细结构中分配的PA内存
        for(auto itr = this->m_dlcl_specific_structure_vector.begin() ; itr != this->m_dlcl_specific_structure_vector.end() ; itr++){
            if(itr->m_left_pa_ptr != nullptr){
                delete itr->m_left_pa_ptr;
            }
        }
        //把所有的MLCL详细结构清空
        list<DlclSpecificStructure>().swap(this->m_dlcl_specific_structure_vector);
    }
}

