#include <CL.h>

//float ClSpecificStructure::m_fragment_score_weight = 1;
//float ClSpecificStructure::m_fa_consistency_score_weight = 1;
//float ClSpecificStructure::m_pa_exist_score_weight = 1;

using namespace std;


Cl::Cl()
{

}

Cl::Cl(Cardiolipin &superclass)
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


void Cl::splice()
{
    //如果是空对象，则不进行操作
    if(this->CheckEmptyObject()){
        return;
    }

//        if(this->GetChainLength() == 70 && this->GetUnsaturation() == 6 && this->GetOxygen() == 0){
//            qDebug() << "stop";
//        }

    //对不同的情况进行拼接操作
    std::vector<Ms2*> ms2_ptr_vector = this->GetMs2VectorPtr();
    for(auto ms2_itr = ms2_ptr_vector.begin() ; ms2_itr != ms2_ptr_vector.end();ms2_itr++){
        //如果这个二级没有PA，尝试用4个FA拼接
        if((*ms2_itr)->GetPaCount() == 0){
            this->FourFaSpliceCl(*ms2_itr);
        }
        //如果这个二级没有FA，则尝试用2个PA拼接
        else if((*ms2_itr)->GetFaCount() == 0){
            this->TwoPaSpliceCl(*ms2_itr);
        }
        //如果这个二级即有FA，又有PA
        else{
            this->TwoPaFourFaSpliceCl(*ms2_itr);
        }
    }

    //去除强度低的拼接结果
    this->DeleteRedundantSpliceResult();

    //如果没有拼接成功，则清空对象
    if(this->m_cl_specific_structure_vector.size() == 0){
        this->EmptyObject();
        return;
    }

    //合并不同二级之间的拼接结果
    this->MergeSplice();

}


void Cl::MergeSplice()
{
    //    //如果只有一个二级，则不需要合并，因为一个二级中的所有拼接结果都是不重复的
    //    if(this->m_Ms2_vector_ptr.size() == 1){
    //        return;
    //    }

    list<ClSpecificStructure> non_repeating_only_pa_list;//存储只有PA和FA的对象，用list来存储
    list<ClSpecificStructure> non_repeating_only_fa_list;//存储只有PA和FA的对象，用list来存储
    list<ClSpecificStructure> non_repeating_have_both_pa_fa_list;//存储既有PA和又有FA的对象，用list来存储

    list<ClSpecificStructure> compare_to_have_both_pa_fa_list;


    //先把只有PA的跟只有PA的合并，只有FA的跟只有FA的合并，两者都有的更两者都有的合并；分别存储到non_repeating_only_pa_list和non_repeating_only_fa_list和non_repeating_have_both_pa_fa_list中
    for(auto first_itr = this->m_cl_specific_structure_vector.begin() ; first_itr != this->m_cl_specific_structure_vector.end();first_itr++){
        if(!first_itr->m_fa_exist){
            for(auto second_itr = next(first_itr) ; second_itr != this->m_cl_specific_structure_vector.end();){
                if(!second_itr->m_fa_exist){
                    //如果这两个Cl相同，进行这个操作的时候，first_itr对应的Cl内部会发生改变
                    if(*first_itr == *second_itr){
                        //如果这两个Cl相同，那么从m_cl_specific_structure_vector中去除第二个Cl
                        second_itr =  this->m_cl_specific_structure_vector.erase(second_itr);
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
            for(auto second_itr = next(first_itr) ; second_itr != this->m_cl_specific_structure_vector.end();){
                if(!second_itr->m_pa_exist){
                    //如果这两个Cl相同，进行这个操作的时候，first_itr对应的Cl内部会发生改变
                    if(*first_itr == *second_itr){
                        second_itr =  this->m_cl_specific_structure_vector.erase(second_itr);
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
            for(auto second_itr = next(first_itr) ; second_itr != this->m_cl_specific_structure_vector.end();){
                if(second_itr->m_fa_exist && second_itr->m_pa_exist){
                    //如果这两个Cl相同，进行这个操作的时候，first_itr对应的Cl内部会发生改变
                    if(*first_itr == *second_itr){
                        second_itr =  this->m_cl_specific_structure_vector.erase(second_itr);
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

    //进行仅有PA的和仅有FA的CL的重复合并，最终结果追加到compare_to_have_both_pa_fa_list中
    set<ClSpecificStructure*> fa_can_be_spliced_with_pa;//判断哪些FA被合并了
    //合并仅有PA的和仅有FA的,non_repeating_only_pa_list和non_repeating_only_fa_list进行合并
    for(auto first_itr = non_repeating_only_pa_list.begin() ; first_itr != non_repeating_only_pa_list.end();first_itr++){
        bool first_itr_success = 0;//判定这个Cl是否能在后面找到成功拼接上的
        for(auto second_itr = non_repeating_only_fa_list.begin() ; second_itr != non_repeating_only_fa_list.end();second_itr++){
            //如果这两个Cl相同，返回新的Cl的共享指针，如果不同，返回nullptr
            shared_ptr<ClSpecificStructure> merge_result = first_itr->merge(&*second_itr);
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


    list<ClSpecificStructure>().swap(this->m_cl_specific_structure_vector);//清空m_cl_specific_structure_vector的内容


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
        this->m_cl_specific_structure_vector.push_back(*first_itr);
    }
    this->m_cl_specific_structure_vector.splice(this->m_cl_specific_structure_vector.end() , compare_to_have_both_pa_fa_list);

    //    if(this->GetChainLength() == 70 && this->GetUnsaturation() == 4){
    //        for(auto itr = this->m_cl_specific_structure_vector.begin() ; itr != this->m_cl_specific_structure_vector.end() ; itr++){
    //            itr->Print();
    //        }
    //    }
    //问题应该不出在这里
}


void Cl::FourFaSpliceCl(Ms2* ms2_ptr)
{
    //以FA的链长对fa_vector进行排序
    ms2_ptr->SortFaByChainLength();
    vector<Fa>* fa_vector_ptr = ms2_ptr->GetFaVectorPtr();


    int fa_vector_size = fa_vector_ptr->size();//获取有多少个FA
    //四数之和解题
    for(int i = 0 ; i < fa_vector_size ; i++){
        for(int j = i ; j < fa_vector_size ; j++){
            int left = j;
            int right = fa_vector_size - 1;
            while(left <= right){
                if(((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) == (this->GetChainLength()) &&
                    ((fa_vector_ptr->at(i).GetUnsaturation() + fa_vector_ptr->at(j).GetUnsaturation() + fa_vector_ptr->at(left).GetUnsaturation() + fa_vector_ptr->at(right).GetUnsaturation()) == (this->GetUnsaturation()) &&
                     ((fa_vector_ptr->at(i).GetOxygen() + fa_vector_ptr->at(j).GetOxygen() + fa_vector_ptr->at(left).GetOxygen() + fa_vector_ptr->at(right).GetOxygen()) == (this->GetOxygen()))))){
                    this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(&fa_vector_ptr->at(i) , &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left) , &fa_vector_ptr->at(right) , ms2_ptr));//加入到结果中
                    int left_t = left;
                    while((left_t < fa_vector_size - 1) && (fa_vector_ptr->at(left_t).GetChainLength() == fa_vector_ptr->at(left_t + 1).GetChainLength()) && (fa_vector_ptr->at(left_t).GetUnsaturation() == fa_vector_ptr->at(left_t + 1).GetUnsaturation()) && (fa_vector_ptr->at(left_t).GetOxygen() == fa_vector_ptr->at(left_t + 1).GetOxygen())){
                        this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(&fa_vector_ptr->at(i) , &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left_t + 1) , &fa_vector_ptr->at(right) , ms2_ptr));//加入到结果中
                        left_t++;
                    }
                    right--;
                }
                //如果所有链长加起来大于目标链长，说明大了，right需要向左移
                else if((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) > (this->GetChainLength())){
                    right--;
                }
                //如果所有链长加起来小于目标链长，说明小了，right需要向右移
                else if((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) < (this->GetChainLength())){
                    left++;
                }
                //如果链长加起来等于目标链长，但是其他属性加起来不等于目标属性，则向left向右寻找是否有符合要求的，寻找完right--
                else if((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) == (this->GetChainLength())){
                    int left_t = left;
                    int right_t = right;
                    while(left_t + 1 < right_t){
                        if(((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left_t + 1).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) == (this->GetChainLength())) &&
                                ((fa_vector_ptr->at(i).GetUnsaturation() + fa_vector_ptr->at(j).GetUnsaturation() + fa_vector_ptr->at(left_t + 1).GetUnsaturation() + fa_vector_ptr->at(right).GetUnsaturation()) == (this->GetUnsaturation())) &&
                                ((fa_vector_ptr->at(i).GetOxygen() + fa_vector_ptr->at(j).GetOxygen() + fa_vector_ptr->at(left_t + 1).GetOxygen() + fa_vector_ptr->at(right).GetOxygen()) == (this->GetOxygen()))){
                            this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(&fa_vector_ptr->at(i) , &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left_t + 1) , &fa_vector_ptr->at(right) , ms2_ptr));//加入到结果中
                        }
                        left_t++;
                    }
                    right--;
                }
            }
        }
    }
}

void Cl::TwoPaSpliceCl(Ms2* ms2_ptr)
{
    //哈系法解两数之和
    multimap<unsigned int,Pa*> hash_map;//以链长为键，对应的Pa为值构建多重哈希表
    vector<Pa>* pa_vector_ptr = ms2_ptr->GetPaVectorPtr();

    for(auto pa_vector_itr = pa_vector_ptr->begin();pa_vector_itr != pa_vector_ptr->end();pa_vector_itr++){
        hash_map.insert({pa_vector_itr->GetChainLength() , &(*pa_vector_itr)});//插入数据
        auto search_pair_itr = hash_map.equal_range(this->GetChainLength() - pa_vector_itr->GetChainLength());//寻找剩余链长的键的范围
        //如果没有找到
        if(search_pair_itr.first == search_pair_itr.second){

        }
        else{
            //遍历搜索到的范围
            for(auto match_itr = search_pair_itr.first ; match_itr != search_pair_itr.second ; match_itr++){
                //如果链长相加起来等于目标链长的情况下，如果不饱和度和氧个数加起来都等于目标的这些属性，加入到结果中
                if(((pa_vector_itr->GetUnsaturation() + match_itr->second->GetUnsaturation()) == this->GetUnsaturation()) &&
                        ((pa_vector_itr->GetOxygen() + match_itr->second->GetOxygen()) == this->GetOxygen())){
                    this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(&*pa_vector_itr , match_itr->second , ms2_ptr));//加入到结果中
                }
            }
        }
    }
}

void Cl::TwoPaFourFaSpliceCl(Ms2* ms2_ptr)
{
    //用于存储哪4个FA拼接成了2个PA
    set<set<Fa*>> which_four_fa_splice_two_pa;
    //以FA的链长对fa_vector进行排序
    ms2_ptr->SortFaByChainLength();
    vector<Fa>* fa_vector_ptr = ms2_ptr->GetFaVectorPtr();
    vector<Pa>* pa_vector_ptr = ms2_ptr->GetPaVectorPtr();

    //先寻找有哪两个PA可以拼接成Cl
    multimap<unsigned int,Pa*> hash_map;//以链长为键，对应的Pa为值构建多重哈希表
    for(auto pa_vector_ptr_itr = pa_vector_ptr->begin();pa_vector_ptr_itr != pa_vector_ptr->end();pa_vector_ptr_itr++){
        hash_map.insert({pa_vector_ptr_itr->GetChainLength() , &*pa_vector_ptr_itr});
        auto search_pair_itr = hash_map.equal_range(this->GetChainLength() - pa_vector_ptr_itr->GetChainLength());//寻找剩余链长的键的范围
        //如果没有找到
        if(search_pair_itr.first == search_pair_itr.second){

        }
        else{
            //遍历搜索到的范围
            for(auto match_itr = search_pair_itr.first ; match_itr != search_pair_itr.second ; match_itr++){
                //如果链长相加起来等于目标链长的情况下，如果不饱和度和氧个数加起来都等于目标的这些属性，那么寻找是否有FA可以拼接成这些Pa
                if(((pa_vector_ptr_itr->GetUnsaturation() + match_itr->second->GetUnsaturation()) == this->GetUnsaturation()) &&
                        ((pa_vector_ptr_itr->GetOxygen() + match_itr->second->GetOxygen()) == this->GetOxygen())){
                    which_four_fa_splice_two_pa = this->FourFaSpliceTwoPa(&*pa_vector_ptr_itr , match_itr->second , fa_vector_ptr , which_four_fa_splice_two_pa , ms2_ptr);//四个FA拼接这两个PA，which_four_fa_splice_two_pa存储能拼接成PA的四个FA的信息
                }
            }
        }
    }

    int fa_vector_size = fa_vector_ptr->size();//获取有多少个FA
    //四数之和解题，需要把可以拼接成2个PA的那4个FA的组合排除，用到了set结构
    for(int i = 0 ; i < fa_vector_size ; i++){
        for(int j = i ; j < fa_vector_size ; j++){
            int left = j;
            int right = fa_vector_size - 1;
            while(left <= right){
                if(((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) == (this->GetChainLength())) &&
                        ((fa_vector_ptr->at(i).GetUnsaturation() + fa_vector_ptr->at(j).GetUnsaturation() + fa_vector_ptr->at(left).GetUnsaturation() + fa_vector_ptr->at(right).GetUnsaturation()) == (this->GetUnsaturation())) &&
                        ((fa_vector_ptr->at(i).GetOxygen() + fa_vector_ptr->at(j).GetOxygen() + fa_vector_ptr->at(left).GetOxygen() + fa_vector_ptr->at(right).GetOxygen()) == (this->GetOxygen()))){
                    //如果找不到，说明这四个FA无法拼接成那2个PA，把这个四个FA的组合情况加入到结果中
                    if(which_four_fa_splice_two_pa.find({&fa_vector_ptr->at(i) , &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left) , &fa_vector_ptr->at(right)}) == which_four_fa_splice_two_pa.end()){
                        this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(&fa_vector_ptr->at(i) , &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left) , &fa_vector_ptr->at(right) , ms2_ptr));//加入到结果中
                    }
                    int left_t = left;
                    while((left_t < fa_vector_size - 1)  && (fa_vector_ptr->at(left_t).GetChainLength() == fa_vector_ptr->at(left_t + 1).GetChainLength()) && (fa_vector_ptr->at(left_t).GetUnsaturation() == fa_vector_ptr->at(left_t + 1).GetUnsaturation()) && (fa_vector_ptr->at(left_t).GetOxygen() == fa_vector_ptr->at(left_t + 1).GetOxygen())){
                        //如果找不到，说明这四个FA无法拼接成那2个PA，把这个四个FA的组合情况加入到结果中
                        if(which_four_fa_splice_two_pa.find({&fa_vector_ptr->at(i) , &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left_t + 1) , &fa_vector_ptr->at(right)}) == which_four_fa_splice_two_pa.end()){

                            this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(&fa_vector_ptr->at(i), &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left + 1) , &fa_vector_ptr->at(right) , ms2_ptr));//加入到结果中
                            //                            if(this->GetChainLength() == 68 && this->GetUnsaturation() == 4){
                            //                                this->m_cl_specific_structure_vector.back().Print();

                            //                            }
                        }
                        left_t++;
                    }
                    right--;
                }
                //如果所有链长加起来大于目标链长，说明大了，right需要向左移
                else if((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) > (this->GetChainLength())){
                    right--;
                }
                //如果所有链长加起来小于目标链长，说明小了，right需要向右移
                else if((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) < (this->GetChainLength())){
                    left++;
                }
                //如果链长加起来等于目标链长，但是其他属性加起来不等于目标属性，则向left向右寻找是否有符合要求的，寻找完right--
                else if((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) == (this->GetChainLength())){
                    int left_t = left;
                    int right_t = right;
                    while(left_t + 1 < right_t){
                        if(((fa_vector_ptr->at(i).GetChainLength() + fa_vector_ptr->at(j).GetChainLength() + fa_vector_ptr->at(left_t + 1).GetChainLength() + fa_vector_ptr->at(right).GetChainLength()) == (this->GetChainLength())) &&
                                ((fa_vector_ptr->at(i).GetUnsaturation() + fa_vector_ptr->at(j).GetUnsaturation() + fa_vector_ptr->at(left_t + 1).GetUnsaturation() + fa_vector_ptr->at(right).GetUnsaturation()) == (this->GetUnsaturation())) &&
                                ((fa_vector_ptr->at(i).GetOxygen() + fa_vector_ptr->at(j).GetOxygen() + fa_vector_ptr->at(left_t + 1).GetOxygen() + fa_vector_ptr->at(right).GetOxygen()) == (this->GetOxygen()))){
                            //如果找不到，说明这四个FA无法拼接成那2个PA，把这个四个FA的组合情况加入到结果中
                            if(which_four_fa_splice_two_pa.find({&fa_vector_ptr->at(i) , &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left + 1) , &fa_vector_ptr->at(right)}) == which_four_fa_splice_two_pa.end()){
                                this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(&fa_vector_ptr->at(i) , &fa_vector_ptr->at(j) , &fa_vector_ptr->at(left_t + 1) , &fa_vector_ptr->at(right) , ms2_ptr));//加入到结果中
                            }
                        }
                        left_t++;
                    }
                    right--;
                }
            }
        }
    }
}

std::set<std::set<Fa*>>& Cl::FourFaSpliceTwoPa(Pa* pa_1_ptr , Pa* pa_2_ptr , std::vector<Fa>* fa_vector_ptr , std::set<std::set<Fa*>>& store , Ms2* ms2_ptr)
{
    //记录哪两个FA拼接成了这两个PA
    vector<pair<Fa*,Fa*>> tow_fa_splice_first_pa;

    multimap<unsigned int , Fa*> first_pa_hash_map;//以链长为键，对应的Fa为值构建多重哈希表，先拼接第一个PA
    for(auto fa_vector_itr = fa_vector_ptr->begin() ; fa_vector_itr != fa_vector_ptr->end();fa_vector_itr++){
        first_pa_hash_map.insert({fa_vector_itr->GetChainLength() , &(*fa_vector_itr)});//插入数据
        auto search_first_pa_pair_itr = first_pa_hash_map.equal_range(pa_1_ptr->GetChainLength() - fa_vector_itr->GetChainLength());//寻找剩余链长的键的范围
        //如果没有找到
        if(search_first_pa_pair_itr.first == search_first_pa_pair_itr.second){

        }
        else{
            //遍历搜索到的范围
            for(auto match_itr = search_first_pa_pair_itr.first ; match_itr != search_first_pa_pair_itr.second ; match_itr++){
                //如果链长相加起来等于目标链长的情况下，如果不饱和度和氧个数加起来都等于目标的这些属性，那么就找到了可以拼接成第一个Pa的Fa
                if(((fa_vector_itr->GetUnsaturation() + match_itr->second->GetUnsaturation()) == pa_1_ptr->GetUnsaturation()) &&
                        ((fa_vector_itr->GetOxygen() + match_itr->second->GetOxygen()) == pa_1_ptr->GetOxygen())){
                    tow_fa_splice_first_pa.emplace_back(pair<Fa*,Fa*>(&*fa_vector_itr , match_itr->second));//记录哪两个FA拼接成了第一个PA
                }
            }
        }
    }

    //如果无法拼接成第一个PA，则直接返回结果
    if(tow_fa_splice_first_pa.size() == 0){
        this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(pa_1_ptr , pa_2_ptr , ms2_ptr));//追加两个PA的拼接结果
        return store;
    }

    //判定是否有4个FA拼接成2个PA成功
    bool four_fa_splice_success = 0;

    //拼接第二个PA
    for(auto fa_vector_itr = fa_vector_ptr->begin() ; fa_vector_itr != fa_vector_ptr->end();fa_vector_itr++){
        //不需要再继续创建新的哈希表，直接用拼接第一个PA的哈希表即可
        auto search_second_pa_pair_itr = first_pa_hash_map.equal_range(pa_2_ptr->GetChainLength() - fa_vector_itr->GetChainLength());//寻找剩余链长的键的范围
        //如果没有找到
        if((search_second_pa_pair_itr.first == first_pa_hash_map.end()) && (search_second_pa_pair_itr.second == first_pa_hash_map.end())){

        }
        else{
            //设定为成功
            four_fa_splice_success = 1;
            //遍历搜索到的范围
            for(auto match_itr = search_second_pa_pair_itr.first ; match_itr != search_second_pa_pair_itr.second ; match_itr++){
                //如果链长相加起来等于目标链长的情况下，如果不饱和度和氧个数加起来都等于目标的这些属性，那么就找到了可以拼接成第二个Pa的Fa
                if(((fa_vector_itr->GetUnsaturation() + match_itr->second->GetUnsaturation()) == pa_2_ptr->GetUnsaturation()) &&
                        ((fa_vector_itr->GetOxygen() + match_itr->second->GetOxygen()) == pa_2_ptr->GetOxygen())){
                    //遍历拼接第一个PA的结果，把结果追加到m_cl_specific_structure_vector和store中
                    for(auto tow_fa_splice_first_pa_itr = tow_fa_splice_first_pa.begin() ; tow_fa_splice_first_pa_itr != tow_fa_splice_first_pa.end() ; tow_fa_splice_first_pa_itr++){

                        //                        set<set<Fa*>>::iterator search =  store.find({tow_fa_splice_first_pa_itr->first , tow_fa_splice_first_pa_itr->second , &*fa_vector_itr , match_itr->second});
                        store.insert({tow_fa_splice_first_pa_itr->first , tow_fa_splice_first_pa_itr->second , &*fa_vector_itr , match_itr->second});//追加这些可以拼接成PA的FA组合
                        this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(pa_1_ptr , tow_fa_splice_first_pa_itr->first , tow_fa_splice_first_pa_itr->second , pa_2_ptr , &*fa_vector_itr , match_itr->second , ms2_ptr));
                    }
                }
            }
        }
    }
    //store有可能数目会比this的m_cl_specific_structure_vector少，因为store是个set，没有重复项；如果有2个重复的FA，另外2个是互相不重复的，游有可能这样可以组成2种PA组合

    if(!four_fa_splice_success){
        //如果第一个成功拼接，但是第二个无法拼接，则把PA追加到结果中
        this->m_cl_specific_structure_vector.emplace_back(ClSpecificStructure(pa_1_ptr , pa_2_ptr , ms2_ptr));//追加两个PA的拼接结果
    }

    return store;
}

void Cl::EmptyObject()
{
    if(this->GetMs1Ptr() != nullptr){
        this->SetMs1Ptr(nullptr);
    }
    if(this->GetMs1DatabaseRecordPtr() != nullptr){
        this->SetMs1DatabaseRecordPtr(nullptr);
    }

    this->SetMs1MatchingScosre(0);
    this->ClearMs2Info();//清空vector中的元素

    if(this->m_cl_specific_structure_vector.size() != 0){
        //去除每个CL详细结构中分配的PA内存
        for(auto itr = this->m_cl_specific_structure_vector.begin() ; itr != this->m_cl_specific_structure_vector.end() ; itr++){
            if(itr->m_left_pa_ptr != nullptr){
                delete itr->m_left_pa_ptr;
            }
            if(itr->m_right_pa_ptr != nullptr){
                delete itr->m_right_pa_ptr;
            }
        }
        //把所有的CL详细结构清空
        this->ClearSpliceResult();
    }
}

void Cl::ClearSpliceResult()
{
    //把所有的CL详细结构清空
    list<ClSpecificStructure>().swap(this->m_cl_specific_structure_vector);
}

void Cl::DeleteRedundantSpliceResult()
{
    //删除拼接结果总强度小于二级总强度5%的拼接结果
    for(auto itr = this->m_cl_specific_structure_vector.begin() ; itr != this->m_cl_specific_structure_vector.end();){
        if((itr->GetTotalIntensity() / itr->GetMs2TotalIntensity()) <= this->m_delete_redundant_splice_result_radio){
            itr = this->m_cl_specific_structure_vector.erase(itr);
        }
        else{
            itr++;
        }
    }
}
