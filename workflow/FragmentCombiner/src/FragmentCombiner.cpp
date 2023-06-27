#include <FragmentCombiner.h>

using namespace std;
string FragmentCombiner::mode = "strict";

FragmentCombiner::FragmentCombiner()
{

}

FragmentCombiner::FragmentCombiner(FragmentFinder& fragment_finder , string mode)
{
    //把fragment_finder中的结果复制过来，不采用指针的方法是为了每一步在在执行之后，如果再执行，用的还是上一步的数据；所以workflow中没有使用指针
    this->m_cl_vector = fragment_finder.GetCopyClPairVector();
    this->m_mlcl_vector = fragment_finder.GetCopyMlclPairVector();
    this->m_dlcl_vector = fragment_finder.GetCopyDlclPairVector();

    this->mode = mode;
}

void FragmentCombiner::CopyInfoFromMs2WithPaAndFa(FragmentFinder& fragment_finder)
{
    //把fragment_finder中的结果复制过来，不采用指针的方法是为了每一步在在执行之后，如果再执行，用的还是上一步的数据；所以workflow中没有使用指针
    this->m_cl_vector = fragment_finder.GetCopyClPairVector();
    this->m_mlcl_vector = fragment_finder.GetCopyMlclPairVector();
    this->m_dlcl_vector = fragment_finder.GetCopyDlclPairVector();
}

void FragmentCombiner::splice()
{
    //cl
    for(auto itr = this->m_cl_vector.begin() ; itr != this->m_cl_vector.end();){

        if(itr->first.GetMs1DatabaseRecordComponent() == array<unsigned int , 3>({58,3,0})){
            qDebug() << "stop";
        }

        itr->first.splice();
        itr->second.splice();
        //如果两个都没有拼接成功；或者一个原本是没有的，另一个没有拼接成功，则将其删除
        if(itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            itr = this->m_cl_vector.erase(itr);
        }
        else{
            ++itr;
        }
    }

    //mlcl
    for(auto itr = this->m_mlcl_vector.begin() ; itr != this->m_mlcl_vector.end();){
        itr->first.splice();
        itr->second.splice();
        //如果两个都没有拼接成功；或者一个原本是没有的，另一个没有拼接成功，则将其删除
        if(itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            itr = this->m_mlcl_vector.erase(itr);
        }
        else{
            ++itr;
        }
    }

    //dlcl
    for(auto itr = this->m_dlcl_vector.begin() ; itr != this->m_dlcl_vector.end();){
        itr->first.splice();
        itr->second.splice();
        //如果两个都没有拼接成功；或者一个原本是没有的，另一个没有拼接成功，则将其删除
        if(itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            itr = this->m_dlcl_vector.erase(itr);
        }
        else{
            ++itr;
        }
    }

    this->DeleteRedundantPair();//删除冗余心磷脂

    //合并M-H和M-2H
    for(auto itr = this->m_cl_vector.begin() ; itr != this->m_cl_vector.end() ; itr++){
        this->MergeClPair(*itr);
    }

    //合并M-H和M-2H
    for(auto itr = this->m_mlcl_vector.begin() ; itr != this->m_mlcl_vector.end() ; itr++){
        this->MergeMlclPair(*itr);
    }

    //合并M-H和M-2H
    for(auto itr = this->m_dlcl_vector.begin() ; itr != this->m_dlcl_vector.end() ; itr++){
        this->MergeDlclPair(*itr);
    }


    //对哈希表中的每个拼接结果的分数从大到小进行排序
   for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end() ; itr++){
       itr->second->sort([](const ClSpecificStructure& lhs, const ClSpecificStructure& rhs){
           return lhs.m_score > rhs.m_score;
       });
   }
   for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end() ; itr++){
       itr->second->sort([](const MlclSpecificStructure& lhs, const MlclSpecificStructure& rhs){
           return lhs.m_score > rhs.m_score;
       });
   }
   for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end() ; itr++){
       itr->second->sort([](const DlclSpecificStructure& lhs, const DlclSpecificStructure& rhs){
           return lhs.m_score > rhs.m_score;
       });
   }



}

void FragmentCombiner::MergeClPair(std::pair<Cl, Cl> &cl_pair)
{
    //如果M-H或者M-2H某个是空的，则把结果加入到m_cl_merge_hash_table中，退出函数
    if(cl_pair.first.CheckEmptyObject()){
        this->m_cl_merge_hash_table.insert({&cl_pair , &cl_pair.second.m_cl_specific_structure_vector});
        return;
    }
    else if(cl_pair.second.CheckEmptyObject())
    {
        this->m_cl_merge_hash_table.insert({&cl_pair , &cl_pair.first.m_cl_specific_structure_vector});
        return;
    }


    auto m_h_cl_specific_structure_vector_copy = cl_pair.first.m_cl_specific_structure_vector;
    auto m_2h_cl_specific_structure_vector_copy = cl_pair.second.m_cl_specific_structure_vector;

    //新建merge_pair_cl_specific_structure_vector_ptr指针，加入M-H和M-2H的拼接结果
    list<ClSpecificStructure>* merge_pair_cl_specific_structure_vector_ptr = new list<ClSpecificStructure>;


    //如果是严格模式，则跟CL中的正常合并的流程相同，但是打分不同，所以合并用StrictMerge
    if(this->mode == "strict"){
        set<ClSpecificStructure*> m_2h_fa_merged_set;//存储M-2H中可以被用来合并的FA
        set<ClSpecificStructure*> m_2h_pa_merged_set;//存储M-2H中可以被用来合并的PA
        for(auto first_itr = m_h_cl_specific_structure_vector_copy.begin() ; first_itr != m_h_cl_specific_structure_vector_copy.end(); first_itr++){
            //如果遍历的M-H是两者都有
            if(first_itr->m_fa_exist && first_itr->m_pa_exist){
                for(auto second_itr = m_2h_cl_specific_structure_vector_copy.begin() ; second_itr != m_2h_cl_specific_structure_vector_copy.end();){
                    //如果成功拼接
                    if(first_itr->StrictMerge(*second_itr)){
                       //如果这个M-2H的这个CL既有PA又有FA，则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                       if(second_itr->m_fa_exist && second_itr->m_pa_exist){
                           second_itr = m_2h_cl_specific_structure_vector_copy.erase(second_itr);
                       }
                       //如果M-2H的这个CL只有PA，那么M-H之后的FA仍然有可能与之拼接成Cl
                       else if(!second_itr->m_fa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                       }
                       //如果M-2H的这个CL只有FA，那么M-H之后的PA仍然有可能与之拼接成Cl
                       else if(!second_itr->m_pa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                       }
                       break;//退出循环
                   }
                    //如果无法成功拼接
                   else{
                       second_itr++;
                   }
                }
                merge_pair_cl_specific_structure_vector_ptr->push_back(*first_itr);
            }

            //如果遍历的M-H只有PA
            else if(!first_itr->m_fa_exist){
                bool splice_success = 0;
                for(auto second_itr = m_2h_cl_specific_structure_vector_copy.begin() ; second_itr != m_2h_cl_specific_structure_vector_copy.end();){
                    shared_ptr<ClSpecificStructure> merge_result = first_itr->merge(&*second_itr , 1);
                    //如果成功拼接
                    if(merge_result != nullptr){
                        splice_success = 1;
                        merge_pair_cl_specific_structure_vector_ptr->push_back(*merge_result);
                        //如果M-2H的这个CL只有PA,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了；M-H后的FA不可能可以拼接成这个仅有PA的CL，因为如果有，那么这个FA肯定能和
                        //first_itr对应的PA早已经拼接了，所以没有
                        if(!second_itr->m_fa_exist){
                            second_itr = m_2h_cl_specific_structure_vector_copy.erase(second_itr);
                            break;
                        }
                        //如果M-2H的这个CL只有FA，M-H后面仍然存在只有PA的与之可以拼接，所以加入到集合
                        else if(!second_itr->m_pa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                        }
                        //如果M-2H的这个CL两者都有,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                        else if(second_itr->m_pa_exist && second_itr->m_fa_exist){
                            second_itr = m_2h_cl_specific_structure_vector_copy.erase(second_itr);
                        }
                    }
                    //如果无法成功拼接
                    else{
                        second_itr++;
                    }
                }
                //如果这个M-H的CL在M-2H总无法找到重复CL，则将这个CL加入到最终结果中
                if(!splice_success){
                    merge_pair_cl_specific_structure_vector_ptr->push_back(*first_itr);
                }
            }

            //如果遍历的M-H只有FA
            else if(!first_itr->m_pa_exist){
                bool splice_success = 0;
                for(auto second_itr = m_2h_cl_specific_structure_vector_copy.begin() ; second_itr != m_2h_cl_specific_structure_vector_copy.end();){
                    shared_ptr<ClSpecificStructure> merge_result = first_itr->merge(&*second_itr , 1);
                    //如果成功合并
                    if(merge_result != nullptr){
                        splice_success = 1;
                        merge_pair_cl_specific_structure_vector_ptr->push_back(*merge_result);
                        //如果M-2H的这个CL只有FA
                        if(!second_itr->m_pa_exist){
                            second_itr = m_2h_cl_specific_structure_vector_copy.erase(second_itr);
                            break;
                        }
                        //如果M-2H的这个CL只有PA，M-H后面仍然存在只有FA的与之可以拼接，所以加入到集合
                        else if(!second_itr->m_fa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                        }
                        //如果M-2H的这个CL两者都有,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                        else if(second_itr->m_pa_exist && second_itr->m_fa_exist){
                            second_itr = m_2h_cl_specific_structure_vector_copy.erase(second_itr);
                        }
                    }
                    //如果无法成功合并
                    else{
                        second_itr++;
                    }
                }
                //如果这个M-H的CL在M-2H总无法找到重复CL，则将这个CL加入到最终结果中
                if(!splice_success){
                    merge_pair_cl_specific_structure_vector_ptr->push_back(*first_itr);
                }
            }
        }

        //去除M-2H中被用于合并的仅有FA或仅有PA的CL
        for(auto itr = m_2h_cl_specific_structure_vector_copy.begin() ; itr != m_2h_cl_specific_structure_vector_copy.end();){
            //如果这个M-2H的CL只有PA
            if(!itr->m_fa_exist){
                if(m_2h_pa_merged_set.find(&*itr) == m_2h_pa_merged_set.end()){
                    itr = m_2h_cl_specific_structure_vector_copy.erase(itr);
                }
                else{
                    itr++;
                }
            }
            //如果这个M-2H的CL只有FA
            else if(!itr->m_pa_exist){
                if(m_2h_pa_merged_set.find(&*itr) == m_2h_pa_merged_set.end()){
                    itr = m_2h_cl_specific_structure_vector_copy.erase(itr);
                }
                else{
                    itr++;
                }
            }
            else{
                itr++;
            }
        }

        //把最终结果和M-2H中无法与M-H合并的部分，相连接
        merge_pair_cl_specific_structure_vector_ptr->splice(merge_pair_cl_specific_structure_vector_ptr->end() , m_2h_cl_specific_structure_vector_copy);

        //加入到哈希表中
        this->m_cl_merge_hash_table.insert({&cl_pair , merge_pair_cl_specific_structure_vector_ptr});
    }

    //如果是宽松模式
    else if(this->mode == "flexible"){
        set<ClSpecificStructure*> m_2h_merged_set;//存储M-2H中可以被用来合并的PA
        //M-H和M-2H合并，M-2H被合并过的Cl加入set中
        for(auto first_itr = m_h_cl_specific_structure_vector_copy.begin() ; first_itr != m_h_cl_specific_structure_vector_copy.end(); first_itr++){
            for(auto second_itr = m_2h_cl_specific_structure_vector_copy.begin() ; second_itr != m_2h_cl_specific_structure_vector_copy.end() ; second_itr++){
                //合并
                if(first_itr->FlexibleMerge(*second_itr)){
                    m_2h_merged_set.insert(&*second_itr);//如果成功被合并，则加入set中
                }
                else{

                }
            }
            merge_pair_cl_specific_structure_vector_ptr->push_back(*first_itr);//加入到最终结果中
        }

        //M-2H去除被合并过的Cl
        for(auto itr = m_2h_cl_specific_structure_vector_copy.begin() ; itr != m_2h_cl_specific_structure_vector_copy.end();){
            //如果这个Cl被合并过
            if(m_2h_merged_set.find(&*itr) != m_2h_merged_set.end()){
                itr = m_2h_cl_specific_structure_vector_copy.erase(itr);//删除
            }
            else{
                itr++;
            }
        }

        //把M-2H不能被合并的部分加入到最终结果中
        merge_pair_cl_specific_structure_vector_ptr->splice(merge_pair_cl_specific_structure_vector_ptr->end() , m_2h_cl_specific_structure_vector_copy);

        //加入到哈希表中
        this->m_cl_merge_hash_table.insert({&cl_pair , merge_pair_cl_specific_structure_vector_ptr});
    }
}

void FragmentCombiner::MergeMlclPair(std::pair<Mlcl, Mlcl> &mlcl_pair)
{
    //如果M-H或者M-2H某个是空的，则把结果加入到m_cl_merge_hash_table中，退出函数
    if(mlcl_pair.first.CheckEmptyObject()){
        this->m_mlcl_merge_hash_table.insert({&mlcl_pair , &mlcl_pair.second.m_mlcl_specific_structure_vector});
        return;
    }
    else if(mlcl_pair.second.CheckEmptyObject())
    {
        this->m_mlcl_merge_hash_table.insert({&mlcl_pair , &mlcl_pair.first.m_mlcl_specific_structure_vector});
        return;
    }


    auto m_h_mlcl_specific_structure_vector_copy = mlcl_pair.first.m_mlcl_specific_structure_vector;
    auto m_2h_mlcl_specific_structure_vector_copy = mlcl_pair.second.m_mlcl_specific_structure_vector;

    //新建merge_pair_cl_specific_structure_vector_ptr指针，加入M-H和M-2H的拼接结果
    list<MlclSpecificStructure>* merge_pair_mlcl_specific_structure_vector_ptr = new list<MlclSpecificStructure>;


    //如果是严格模式，则跟CL中的正常合并的流程相同，但是打分不同，所以合并用StrictMerge
    if(this->mode == "strict"){
        set<MlclSpecificStructure*> m_2h_fa_merged_set;//存储M-2H中可以被用来合并的FA
        set<MlclSpecificStructure*> m_2h_pa_merged_set;//存储M-2H中可以被用来合并的PA
        for(auto first_itr = m_h_mlcl_specific_structure_vector_copy.begin() ; first_itr != m_h_mlcl_specific_structure_vector_copy.end(); first_itr++){
            //如果遍历的M-H是两者都有
            if(first_itr->m_fa_exist && first_itr->m_pa_exist){
                for(auto second_itr = m_2h_mlcl_specific_structure_vector_copy.begin() ; second_itr != m_2h_mlcl_specific_structure_vector_copy.end();){
                    //如果成功拼接
                    if(first_itr->StrictMerge(*second_itr)){
                       //如果这个M-2H的这个CL既有PA又有FA，则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                       if(second_itr->m_fa_exist && second_itr->m_pa_exist){
                           second_itr = m_2h_mlcl_specific_structure_vector_copy.erase(second_itr);
                       }
                       //如果M-2H的这个CL只有PA，那么M-H之后的FA仍然有可能与之拼接成Cl
                       else if(!second_itr->m_fa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                       }
                       //如果M-2H的这个CL只有FA，那么M-H之后的PA仍然有可能与之拼接成Cl
                       else if(!second_itr->m_pa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                       }
                       break;//退出循环
                   }
                    //如果无法成功拼接
                   else{
                       second_itr++;
                   }
                }
                merge_pair_mlcl_specific_structure_vector_ptr->push_back(*first_itr);
            }

            //如果遍历的M-H只有PA
            else if(!first_itr->m_fa_exist){
                bool splice_success = 0;
                for(auto second_itr = m_2h_mlcl_specific_structure_vector_copy.begin() ; second_itr != m_2h_mlcl_specific_structure_vector_copy.end();){
                    shared_ptr<MlclSpecificStructure> merge_result = first_itr->merge(&*second_itr , 1);
                    //如果成功拼接
                    if(merge_result != nullptr){
                        splice_success = 1;
                        merge_pair_mlcl_specific_structure_vector_ptr->push_back(*merge_result);
                        //如果M-2H的这个CL只有PA,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了；M-H后的FA不可能可以拼接成这个仅有PA的CL，因为如果有，那么这个FA肯定能和
                        //first_itr对应的PA早已经拼接了，所以没有
                        if(!second_itr->m_fa_exist){
                            second_itr = m_2h_mlcl_specific_structure_vector_copy.erase(second_itr);
                            break;
                        }
                        //如果M-2H的这个CL只有FA，M-H后面仍然存在只有PA的与之可以拼接，所以加入到集合
                        else if(!second_itr->m_pa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                        }
                        //如果M-2H的这个CL两者都有,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                        else if(second_itr->m_pa_exist && second_itr->m_fa_exist){
                            second_itr = m_2h_mlcl_specific_structure_vector_copy.erase(second_itr);
                        }
                    }
                    //如果无法成功拼接
                    else{
                        second_itr++;
                    }
                }
                //如果这个M-H的CL在M-2H总无法找到重复CL，则将这个CL加入到最终结果中
                if(!splice_success){
                    merge_pair_mlcl_specific_structure_vector_ptr->push_back(*first_itr);
                }
            }

            //如果遍历的M-H只有FA
            else if(!first_itr->m_pa_exist){
                bool splice_success = 0;
                for(auto second_itr = m_2h_mlcl_specific_structure_vector_copy.begin() ; second_itr != m_2h_mlcl_specific_structure_vector_copy.end();){
                    shared_ptr<MlclSpecificStructure> merge_result = first_itr->merge(&*second_itr , 1);
                    //如果成功合并
                    if(merge_result != nullptr){
                        splice_success = 1;
                        merge_pair_mlcl_specific_structure_vector_ptr->push_back(*merge_result);
                        //如果M-2H的这个CL只有FA
                        if(!second_itr->m_pa_exist){
                            second_itr = m_2h_mlcl_specific_structure_vector_copy.erase(second_itr);
                            break;
                        }
                        //如果M-2H的这个CL只有PA，M-H后面仍然存在只有FA的与之可以拼接，所以加入到集合
                        else if(!second_itr->m_fa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                        }
                        //如果M-2H的这个CL两者都有,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                        else if(second_itr->m_pa_exist && second_itr->m_fa_exist){
                            second_itr = m_2h_mlcl_specific_structure_vector_copy.erase(second_itr);
                        }
                    }
                    //如果无法成功合并
                    else{
                        second_itr++;
                    }
                }
                //如果这个M-H的CL在M-2H总无法找到重复CL，则将这个CL加入到最终结果中
                if(!splice_success){
                    merge_pair_mlcl_specific_structure_vector_ptr->push_back(*first_itr);
                }
            }
        }

        //去除M-2H中被用于合并的仅有FA或仅有PA的CL
        for(auto itr = m_2h_mlcl_specific_structure_vector_copy.begin() ; itr != m_2h_mlcl_specific_structure_vector_copy.end();){
            //如果这个M-2H的CL只有PA
            if(!itr->m_fa_exist){
                if(m_2h_pa_merged_set.find(&*itr) == m_2h_pa_merged_set.end()){
                    itr = m_2h_mlcl_specific_structure_vector_copy.erase(itr);
                }
                else{
                    itr++;
                }
            }
            //如果这个M-2H的CL只有FA
            else if(!itr->m_pa_exist){
                if(m_2h_pa_merged_set.find(&*itr) == m_2h_pa_merged_set.end()){
                    itr = m_2h_mlcl_specific_structure_vector_copy.erase(itr);
                }
                else{
                    itr++;
                }
            }
            else{
                itr++;
            }
        }

        //把最终结果和M-2H中无法与M-H合并的部分，相连接
        merge_pair_mlcl_specific_structure_vector_ptr->splice(merge_pair_mlcl_specific_structure_vector_ptr->end() , m_2h_mlcl_specific_structure_vector_copy);

        //加入到哈希表中
        this->m_mlcl_merge_hash_table.insert({&mlcl_pair , merge_pair_mlcl_specific_structure_vector_ptr});
    }

    //如果是宽松模式
    else if(this->mode == "flexible"){
        set<MlclSpecificStructure*> m_2h_merged_set;//存储M-2H中可以被用来合并的PA
        //M-H和M-2H合并，M-2H被合并过的Cl加入set中
        for(auto first_itr = m_h_mlcl_specific_structure_vector_copy.begin() ; first_itr != m_h_mlcl_specific_structure_vector_copy.end(); first_itr++){
            for(auto second_itr = m_2h_mlcl_specific_structure_vector_copy.begin() ; second_itr != m_2h_mlcl_specific_structure_vector_copy.end() ; second_itr++){
                //合并
                if(first_itr->FlexibleMerge(*second_itr)){
                    m_2h_merged_set.insert(&*second_itr);//如果成功被合并，则加入set中
                }
                else{

                }
            }
            merge_pair_mlcl_specific_structure_vector_ptr->push_back(*first_itr);//加入到最终结果中
        }

        //M-2H去除被合并过的Cl
        for(auto itr = m_2h_mlcl_specific_structure_vector_copy.begin() ; itr != m_2h_mlcl_specific_structure_vector_copy.end();){
            //如果这个Cl被合并过
            if(m_2h_merged_set.find(&*itr) != m_2h_merged_set.end()){
                itr = m_2h_mlcl_specific_structure_vector_copy.erase(itr);//删除
            }
            else{
                itr++;
            }
        }

        //把M-2H不能被合并的部分加入到最终结果中
        merge_pair_mlcl_specific_structure_vector_ptr->splice(merge_pair_mlcl_specific_structure_vector_ptr->end() , m_2h_mlcl_specific_structure_vector_copy);

        //加入到哈希表中
        this->m_mlcl_merge_hash_table.insert({&mlcl_pair , merge_pair_mlcl_specific_structure_vector_ptr});
    }

}

void FragmentCombiner::MergeDlclPair(std::pair<Dlcl, Dlcl> &dlcl_pair)
{
    //如果M-H或者M-2H某个是空的，则把结果加入到m_cl_merge_hash_table中，退出函数
    if(dlcl_pair.first.CheckEmptyObject()){
        this->m_dlcl_merge_hash_table.insert({&dlcl_pair , &dlcl_pair.second.m_dlcl_specific_structure_vector});
        return;
    }
    else if(dlcl_pair.second.CheckEmptyObject())
    {
        this->m_dlcl_merge_hash_table.insert({&dlcl_pair , &dlcl_pair.first.m_dlcl_specific_structure_vector});
        return;
    }


    auto m_h_dlcl_specific_structure_vector_copy = dlcl_pair.first.m_dlcl_specific_structure_vector;
    auto m_2h_dlcl_specific_structure_vector_copy = dlcl_pair.second.m_dlcl_specific_structure_vector;

    //新建merge_pair_cl_specific_structure_vector_ptr指针，加入M-H和M-2H的拼接结果
    list<DlclSpecificStructure>* merge_pair_dlcl_specific_structure_vector_ptr = new list<DlclSpecificStructure>;


    //如果是严格模式，则跟CL中的正常合并的流程相同，但是打分不同，所以合并用StrictMerge
    if(this->mode == "strict"){
        set<DlclSpecificStructure*> m_2h_fa_merged_set;//存储M-2H中可以被用来合并的FA
        set<DlclSpecificStructure*> m_2h_pa_merged_set;//存储M-2H中可以被用来合并的PA
        for(auto first_itr = m_h_dlcl_specific_structure_vector_copy.begin() ; first_itr != m_h_dlcl_specific_structure_vector_copy.end(); first_itr++){
            //如果遍历的M-H是两者都有
            if(first_itr->m_fa_exist && first_itr->m_pa_exist){
                for(auto second_itr = m_2h_dlcl_specific_structure_vector_copy.begin() ; second_itr != m_2h_dlcl_specific_structure_vector_copy.end();){
                    //如果成功拼接
                    if(first_itr->StrictMerge(*second_itr)){
                       //如果这个M-2H的这个CL既有PA又有FA，则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                       if(second_itr->m_fa_exist && second_itr->m_pa_exist){
                           second_itr = m_2h_dlcl_specific_structure_vector_copy.erase(second_itr);
                       }
                       //如果M-2H的这个CL只有PA，那么M-H之后的FA仍然有可能与之拼接成Cl
                       else if(!second_itr->m_fa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                       }
                       //如果M-2H的这个CL只有FA，那么M-H之后的PA仍然有可能与之拼接成Cl
                       else if(!second_itr->m_pa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                       }
                       break;//退出循环
                   }
                    //如果无法成功拼接
                   else{
                       second_itr++;
                   }
                }
                merge_pair_dlcl_specific_structure_vector_ptr->push_back(*first_itr);
            }

            //如果遍历的M-H只有PA
            else if(!first_itr->m_fa_exist){
                bool splice_success = 0;
                for(auto second_itr = m_2h_dlcl_specific_structure_vector_copy.begin() ; second_itr != m_2h_dlcl_specific_structure_vector_copy.end();){
                    shared_ptr<DlclSpecificStructure> merge_result = first_itr->merge(&*second_itr , 1);
                    //如果成功拼接
                    if(merge_result != nullptr){
                        splice_success = 1;
                        merge_pair_dlcl_specific_structure_vector_ptr->push_back(*merge_result);
                        //如果M-2H的这个CL只有PA,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了；M-H后的FA不可能可以拼接成这个仅有PA的CL，因为如果有，那么这个FA肯定能和
                        //first_itr对应的PA早已经拼接了，所以没有
                        if(!second_itr->m_fa_exist){
                            second_itr = m_2h_dlcl_specific_structure_vector_copy.erase(second_itr);
                            break;
                        }
                        //如果M-2H的这个CL只有FA，M-H后面仍然存在只有PA的与之可以拼接，所以加入到集合
                        else if(!second_itr->m_pa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                        }
                        //如果M-2H的这个CL两者都有,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                        else if(second_itr->m_pa_exist && second_itr->m_fa_exist){
                            second_itr = m_2h_dlcl_specific_structure_vector_copy.erase(second_itr);
                        }
                    }
                    //如果无法成功拼接
                    else{
                        second_itr++;
                    }
                }
                //如果这个M-H的CL在M-2H总无法找到重复CL，则将这个CL加入到最终结果中
                if(!splice_success){
                    merge_pair_dlcl_specific_structure_vector_ptr->push_back(*first_itr);
                }
            }

            //如果遍历的M-H只有FA
            else if(!first_itr->m_pa_exist){
                bool splice_success = 0;
                for(auto second_itr = m_2h_dlcl_specific_structure_vector_copy.begin() ; second_itr != m_2h_dlcl_specific_structure_vector_copy.end();){
                    shared_ptr<DlclSpecificStructure> merge_result = first_itr->merge(&*second_itr , 1);
                    //如果成功合并
                    if(merge_result != nullptr){
                        splice_success = 1;
                        merge_pair_dlcl_specific_structure_vector_ptr->push_back(*merge_result);
                        //如果M-2H的这个CL只有FA
                        if(!second_itr->m_pa_exist){
                            second_itr = m_2h_dlcl_specific_structure_vector_copy.erase(second_itr);
                            break;
                        }
                        //如果M-2H的这个CL只有PA，M-H后面仍然存在只有FA的与之可以拼接，所以加入到集合
                        else if(!second_itr->m_fa_exist){
                            m_2h_pa_merged_set.insert(&*second_itr);
                            second_itr++;
                        }
                        //如果M-2H的这个CL两者都有,则去除M-2H的这个CL，因为M-2H的这个CL已经无法和M-H的其他CL重复了
                        else if(second_itr->m_pa_exist && second_itr->m_fa_exist){
                            second_itr = m_2h_dlcl_specific_structure_vector_copy.erase(second_itr);
                        }
                    }
                    //如果无法成功合并
                    else{
                        second_itr++;
                    }
                }
                //如果这个M-H的CL在M-2H总无法找到重复CL，则将这个CL加入到最终结果中
                if(!splice_success){
                    merge_pair_dlcl_specific_structure_vector_ptr->push_back(*first_itr);
                }
            }
        }

        //去除M-2H中被用于合并的仅有FA或仅有PA的CL
        for(auto itr = m_2h_dlcl_specific_structure_vector_copy.begin() ; itr != m_2h_dlcl_specific_structure_vector_copy.end();){
            //如果这个M-2H的CL只有PA
            if(!itr->m_fa_exist){
                if(m_2h_pa_merged_set.find(&*itr) == m_2h_pa_merged_set.end()){
                    itr = m_2h_dlcl_specific_structure_vector_copy.erase(itr);
                }
                else{
                    itr++;
                }
            }
            //如果这个M-2H的CL只有FA
            else if(!itr->m_pa_exist){
                if(m_2h_pa_merged_set.find(&*itr) == m_2h_pa_merged_set.end()){
                    itr = m_2h_dlcl_specific_structure_vector_copy.erase(itr);
                }
                else{
                    itr++;
                }
            }
            else{
                itr++;
            }
        }

        //把最终结果和M-2H中无法与M-H合并的部分，相连接
        merge_pair_dlcl_specific_structure_vector_ptr->splice(merge_pair_dlcl_specific_structure_vector_ptr->end() , m_2h_dlcl_specific_structure_vector_copy);

        //加入到哈希表中
        this->m_dlcl_merge_hash_table.insert({&dlcl_pair , merge_pair_dlcl_specific_structure_vector_ptr});
    }

    //如果是宽松模式
    else if(this->mode == "flexible"){
        set<DlclSpecificStructure*> m_2h_merged_set;//存储M-2H中可以被用来合并的PA
        //M-H和M-2H合并，M-2H被合并过的Cl加入set中
        for(auto first_itr = m_h_dlcl_specific_structure_vector_copy.begin() ; first_itr != m_h_dlcl_specific_structure_vector_copy.end(); first_itr++){
            for(auto second_itr = m_2h_dlcl_specific_structure_vector_copy.begin() ; second_itr != m_2h_dlcl_specific_structure_vector_copy.end() ; second_itr++){
                //合并
                if(first_itr->FlexibleMerge(*second_itr)){
                    m_2h_merged_set.insert(&*second_itr);//如果成功被合并，则加入set中
                }
                else{

                }
            }
            merge_pair_dlcl_specific_structure_vector_ptr->push_back(*first_itr);//加入到最终结果中
        }

        //M-2H去除被合并过的Cl
        for(auto itr = m_2h_dlcl_specific_structure_vector_copy.begin() ; itr != m_2h_dlcl_specific_structure_vector_copy.end();){
            //如果这个Cl被合并过
            if(m_2h_merged_set.find(&*itr) != m_2h_merged_set.end()){
                itr = m_2h_dlcl_specific_structure_vector_copy.erase(itr);//删除
            }
            else{
                itr++;
            }
        }

        //把M-2H不能被合并的部分加入到最终结果中
        merge_pair_dlcl_specific_structure_vector_ptr->splice(merge_pair_dlcl_specific_structure_vector_ptr->end() , m_2h_dlcl_specific_structure_vector_copy);

        //加入到哈希表中
        this->m_dlcl_merge_hash_table.insert({&dlcl_pair , merge_pair_dlcl_specific_structure_vector_ptr});
    }

}

void FragmentCombiner::OutputResultWithTxt(QString folder_path)
{
    this->Filter();//过滤结果

    QDateTime now_time = QDateTime::currentDateTime();
    QString fileName = folder_path + "/" + "计算机所有定性结果.txt";
    QFile file(fileName);
    file.open(QIODevice::WriteOnly | QIODevice::Text);//打开文件

    QTextStream stream(&file);//文本流
//    stream << "Hello, world!" << endl;
//    stream << "This is a sample text." << endl;

    //输出CL
    for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end();itr++){
        stream << "M-H Info->->" << endl;
        stream << itr->first->first.ShowInfo() << endl;//M-H

        stream << "M-2H Info->->" << endl;
        stream << itr->first->second.ShowInfo() << endl << endl;//M-2H

        stream << "splice Info" << endl;

        for(auto splice_itr = itr->second->begin() ; splice_itr != itr->second->end() ; splice_itr++){
            stream << splice_itr->ShowInfo() << endl;
        }

        stream << "----------------------------------------------------------" << endl;//分割线
    }

    //输出MLCL
    for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end();itr++){
        stream << "M-H Info->->" << endl;
        stream << itr->first->first.ShowInfo() << endl;//M-H

        stream << "M-2H Info->->" << endl;
        stream << itr->first->second.ShowInfo() << endl << endl;//M-2H

        stream << "splice Info" << endl;

        for(auto splice_itr = itr->second->begin() ; splice_itr != itr->second->end() ; splice_itr++){
            stream << splice_itr->ShowInfo() << endl;
        }

        stream << "----------------------------------------------------------" << endl;//分割线
    }

    //输出DLCL
    for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end();itr++){
        stream << "M-H Info->->" << endl;
        stream << itr->first->first.ShowInfo() << endl;//M-H

        stream << "M-2H Info->->" << endl;
        stream << itr->first->second.ShowInfo() << endl << endl;//M-2H

        stream << "splice Info" << endl;

        for(auto splice_itr = itr->second->begin() ; splice_itr != itr->second->end() ; splice_itr++){
            stream << splice_itr->ShowInfo() << endl;
        }

        stream << "----------------------------------------------------------" << endl;//分割线
    }
    file.close();//关闭文件
    qDebug() << "输出了:" << fileName;
}

void FragmentCombiner::OutPutWithCsv(QString folder_path)
{
    this->Filter();//过滤结果

    QString fileName = folder_path + "/" + "计算机所有定性结果.csv";
    //QString fileName = "D:/qt_project/R_analyze/计算机定性结果/dda/计算机所有定性结果.csv";
    QFile file(fileName);
    file.open(QIODevice::WriteOnly | QIODevice::Text);//打开文件

    QTextStream stream(&file);//文本流
    stream << "Compound" << "," << "M.H.m.z" << "," << "M.2H.m.z" << "," << "rt(sec)" << "," << "CL" << "," << "top-1" << "," << "intensity" << endl;
    //输出CL
    for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end();itr++){
        QString component_message = GetPairComponentWithQString<Cl>(*(itr->first));//获得pair的化合物信息
        pair<float,float> sample_mz = GetPairSampleMz<Cl>(*(itr->first));//获得pair的mz，有M-H和M-2H的mz
        float sample_rt = GetPairRt<Cl>(*(itr->first));//获得pair的rt，结果为M-H和M-2H的平均值
        QString best_splice = itr->second->begin()->ShowSimpleInfo();//获得最优拼接
        float sample_intensity = GetPairIntensity(*(itr->first));//获取pair的intensity，结果为M-H和M-2H的最大值
        stream << component_message << "," << sample_mz.first << "," << sample_mz.second << "," << sample_rt << "," << "CL" << "," << best_splice << "," << sample_intensity << endl;
    }

    //输出MLCL
    for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end();itr++){
        QString component_message = GetPairComponentWithQString<Mlcl>(*(itr->first));//获得pair的化合物信息
        pair<float,float> sample_mz = GetPairSampleMz<Mlcl>(*(itr->first));
        float sample_rt = GetPairRt<Mlcl>(*(itr->first));
        QString best_splice = itr->second->begin()->ShowSimpleInfo();
        float sample_intensity = GetPairIntensity(*(itr->first));//获取pair的intensity，结果为M-H和M-2H的最大值
        stream << component_message << "," << sample_mz.first << "," << sample_mz.second << "," << sample_rt << "," << "MLCL" << "," << best_splice << "," << sample_intensity << endl;
    }

    //输出DLCL
    for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end();itr++){
        QString component_message = GetPairComponentWithQString<Dlcl>(*(itr->first));//获得pair的化合物信息
        pair<float,float> sample_mz = GetPairSampleMz<Dlcl>(*(itr->first));
        float sample_rt = GetPairRt<Dlcl>(*(itr->first));
        QString best_splice = itr->second->begin()->ShowSimpleInfo();
        float sample_intensity = GetPairIntensity(*(itr->first));//获取pair的intensity，结果为M-H和M-2H的最大值
        stream << component_message << "," << sample_mz.first << "," << sample_mz.second << "," << sample_rt << "," << "DLCL" << "," << best_splice << "," << sample_intensity << endl;
    }

    file.close();//关闭文件
    qDebug() << "输出了:" << fileName;
}

void FragmentCombiner::Filter()
{
    ifstream fp("D:/R package/human_qualitative.csv"); //定义声明一个ifstream对象，指定文件路径，过滤掉人工定性中存在的部分
    string line;

    set<array<unsigned int,3>> human_qualitative_set;//人工定性的结果
    while (getline(fp,line)){ //循环读取每行数据
        int chain_length_pos = line.find_first_of(",");//找到的第一个','的位置
        unsigned int chain_length = stoul(line.substr(0, chain_length_pos));
        auto left_line = line.substr(chain_length_pos + 1);//取出剩余的字符串
        int unsaturation_pos = left_line.find_first_of(",");//找到的第二个','的位置
        unsigned int unsaturation = stoul(left_line.substr(0, unsaturation_pos));
        unsigned int oxygen = stoul(left_line.substr(unsaturation_pos + 1));
        human_qualitative_set.insert({chain_length , unsaturation , oxygen});//插入set中
    }

//    //过滤人工中不存在的
//    for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end();){
//        //如果M-H或者M-2H某一个找到了人工拼接成功的心磷脂，则删除，为了过滤
//        if(human_qualitative_set.find(itr->first->first.GetMs1DatabaseRecordComponent()) == human_qualitative_set.end() && human_qualitative_set.find(itr->first->second.GetMs1DatabaseRecordComponent()) == human_qualitative_set.end()){
//            itr = this->m_cl_merge_hash_table.erase(itr);//删除
//        }
//        else{
//            itr++;
//        }
//    }
//    for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end();){
//        //如果M-H或者M-2H某一个找到了人工拼接成功的心磷脂，则删除，为了过滤
//        if(human_qualitative_set.find(itr->first->first.GetMs1DatabaseRecordComponent()) == human_qualitative_set.end() && human_qualitative_set.find(itr->first->second.GetMs1DatabaseRecordComponent()) == human_qualitative_set.end()){
//            itr = this->m_mlcl_merge_hash_table.erase(itr);//删除
//        }
//        else{
//            itr++;
//        }
//    }
//    for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end();){
//        //如果M-H或者M-2H某一个找到了人工拼接成功的心磷脂，则删除，为了过滤
//        if(human_qualitative_set.find(itr->first->first.GetMs1DatabaseRecordComponent()) == human_qualitative_set.end() && human_qualitative_set.find(itr->first->second.GetMs1DatabaseRecordComponent()) == human_qualitative_set.end()){
//            itr = this->m_dlcl_merge_hash_table.erase(itr);//删除
//        }
//        else{
//            itr++;
//        }
//    }
//    //过滤人工中存在的
//    for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end();){
//        //如果M-H或者M-2H某一个找到了人工拼接成功的心磷脂，则删除，为了过滤
//        if(human_qualitative_set.find(itr->first->first.GetMs1DatabaseRecordComponent()) != human_qualitative_set.end() || human_qualitative_set.find(itr->first->second.GetMs1DatabaseRecordComponent()) != human_qualitative_set.end()){
//            itr = this->m_cl_merge_hash_table.erase(itr);//删除
//        }
//        else{
//            itr++;
//        }
//    }
//    for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end();){
//        //如果M-H或者M-2H某一个找到了人工拼接成功的心磷脂，则删除，为了过滤
//        if(human_qualitative_set.find(itr->first->first.GetMs1DatabaseRecordComponent()) != human_qualitative_set.end() || human_qualitative_set.find(itr->first->second.GetMs1DatabaseRecordComponent()) != human_qualitative_set.end()){
//            itr = this->m_mlcl_merge_hash_table.erase(itr);//删除
//        }
//        else{
//            itr++;
//        }
//    }
//    for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end();){
//        //如果M-H或者M-2H某一个找到了人工拼接成功的心磷脂，则删除，为了过滤
//        if(human_qualitative_set.find(itr->first->first.GetMs1DatabaseRecordComponent()) != human_qualitative_set.end() || human_qualitative_set.find(itr->first->second.GetMs1DatabaseRecordComponent()) != human_qualitative_set.end()){
//            itr = this->m_dlcl_merge_hash_table.erase(itr);//删除
//        }
//        else{
//            itr++;
//        }
//    }
//    //去除不含氧的心磷脂
//    for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end();){
//        if(itr->first->first.GetMs1DatabaseRecordComponent().at(2) == 0 && itr->first->second.GetMs1DatabaseRecordComponent().at(2) == 0){
//            itr = this->m_cl_merge_hash_table.erase(itr);
//        }
//        else{
//            itr++;
//        }
//    }
//    for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end();){
//        if(itr->first->first.GetMs1DatabaseRecordComponent().at(2) == 0 && itr->first->second.GetMs1DatabaseRecordComponent().at(2) == 0){
//            itr = this->m_mlcl_merge_hash_table.erase(itr);
//        }
//        else{
//            itr++;
//        }
//    }
//    for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end();){
//        if(itr->first->first.GetMs1DatabaseRecordComponent().at(2) == 0 && itr->first->second.GetMs1DatabaseRecordComponent().at(2) == 0){
//            itr = this->m_dlcl_merge_hash_table.erase(itr);
//        }
//        else{
//            itr++;
//        }
//    }

//    //去除含氧的心磷脂
//    for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end();){
//        if(itr->first->first.GetMs1DatabaseRecordComponent().at(2) > 0 || itr->first->second.GetMs1DatabaseRecordComponent().at(2) > 0){
//            itr = this->m_cl_merge_hash_table.erase(itr);
//        }
//        else{
//            itr++;
//        }
//    }
//    for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end();){
//        if(itr->first->first.GetMs1DatabaseRecordComponent().at(2) > 0 || itr->first->second.GetMs1DatabaseRecordComponent().at(2) > 0){
//            itr = this->m_mlcl_merge_hash_table.erase(itr);
//        }
//        else{
//            itr++;
//        }
//    }
//    for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end();){
//        if(itr->first->first.GetMs1DatabaseRecordComponent().at(2) > 0 || itr->first->second.GetMs1DatabaseRecordComponent().at(2) > 0){
//            itr = this->m_dlcl_merge_hash_table.erase(itr);
//        }
//        else{
//            itr++;
//        }
    //    }

    qDebug() << "拼接结果过滤成功";
}


void FragmentCombiner::ShowSpliceResult(QTableWidget *qtablewidget)
{
    //先清空信息
    qtablewidget->setRowCount(0);
    map<unsigned int , std::list<ClSpecificStructure>*>().swap(this->m_row_map_to_cl_splice_result);
    map<unsigned int , std::list<MlclSpecificStructure>*>().swap(this->m_row_map_to_mlcl_splice_result);
    map<unsigned int , std::list<DlclSpecificStructure>*>().swap(this->m_row_map_to_dlcl_splice_result);

    //在这里设置行，可以显著提升ui性能；而不是每个for中增加一行
    qtablewidget->setRowCount(this->m_cl_merge_hash_table.size() + this->m_mlcl_merge_hash_table.size() + this->m_dlcl_merge_hash_table.size());

    unsigned int row = 0;
    //显示CL
    for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end();itr++){
        this->m_row_map_to_cl_splice_result.insert({row , itr->second});//加入map
        if(itr->first->first.CheckEmptyObject()){
            qtablewidget->setItem(row , 0 , 0);
            qtablewidget->setItem(row , 1 , 0);
            qtablewidget->setItem(row , 2 , 0);
        }
        else{
            QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
            QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
            QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
            qtablewidget->setItem(row , 0 , M_H_mz);
            qtablewidget->setItem(row , 1 , M_H_rt);
            qtablewidget->setItem(row , 2 , M_H_struct);
        }

        if(itr->first->second.CheckEmptyObject()){
            qtablewidget->setItem(row , 3 , 0);
            qtablewidget->setItem(row , 4 , 0);
            qtablewidget->setItem(row , 5 , 0);
        }
        else{
            QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
            QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
            QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
            qtablewidget->setItem(row , 3 , M_2H_mz);
            qtablewidget->setItem(row , 4 , M_2H_rt);
            qtablewidget->setItem(row , 5 , M_2H_struct);
        }
        row++;
    }

    //显示MLCL
    for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end();itr++){
        this->m_row_map_to_mlcl_splice_result.insert({row , itr->second});//加入map
        if(itr->first->first.CheckEmptyObject()){
            qtablewidget->setItem(row , 0 , 0);
            qtablewidget->setItem(row , 1 , 0);
            qtablewidget->setItem(row , 2 , 0);
        }
        else{
            QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
            QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
            QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
            qtablewidget->setItem(row , 0 , M_H_mz);
            qtablewidget->setItem(row , 1 , M_H_rt);
            qtablewidget->setItem(row , 2 , M_H_struct);
        }

        if(itr->first->second.CheckEmptyObject()){
            qtablewidget->setItem(row , 3 , 0);
            qtablewidget->setItem(row , 4 , 0);
            qtablewidget->setItem(row , 5 , 0);
        }
        else{
            QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
            QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
            QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
            qtablewidget->setItem(row , 3 , M_2H_mz);
            qtablewidget->setItem(row , 4 , M_2H_rt);
            qtablewidget->setItem(row , 5 , M_2H_struct);
        }
        row++;
    }

    //显示DLCL
    for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end();itr++){
        this->m_row_map_to_dlcl_splice_result.insert({row , itr->second});//加入map
        if(itr->first->first.CheckEmptyObject()){
            qtablewidget->setItem(row , 0 , 0);
            qtablewidget->setItem(row , 1 , 0);
            qtablewidget->setItem(row , 2 , 0);
        }
        else{
            QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
            QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
            QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
            qtablewidget->setItem(row , 0 , M_H_mz);
            qtablewidget->setItem(row , 1 , M_H_rt);
            qtablewidget->setItem(row , 2 , M_H_struct);
        }

        if(itr->first->second.CheckEmptyObject()){
            qtablewidget->setItem(row , 3 , 0);
            qtablewidget->setItem(row , 4 , 0);
            qtablewidget->setItem(row , 5 , 0);
        }
        else{
            QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
            QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
            QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
            qtablewidget->setItem(row , 3 , M_2H_mz);
            qtablewidget->setItem(row , 4 , M_2H_rt);
            qtablewidget->setItem(row , 5 , M_2H_struct);
        }
        row++;
    }
}

void FragmentCombiner::ShowSpliceResult(QTableWidget *qtablewidget, float min_rt, float max_rt)
{
    //先清空信息
    qtablewidget->setRowCount(0);
    //在这里设置行，可以显著提升ui性能；而不是每个for中增加一行
    qtablewidget->setRowCount(this->m_cl_merge_hash_table.size() + this->m_mlcl_merge_hash_table.size() + this->m_dlcl_merge_hash_table.size());

    //清空哈希表
    map<unsigned int , std::list<ClSpecificStructure>*>().swap(this->m_row_map_to_cl_splice_result);
    map<unsigned int , std::list<MlclSpecificStructure>*>().swap(this->m_row_map_to_mlcl_splice_result);
    map<unsigned int , std::list<DlclSpecificStructure>*>().swap(this->m_row_map_to_dlcl_splice_result);

    unsigned int row = 0;

    //显示CL
    for(auto itr = this->m_cl_merge_hash_table.begin() ; itr != this->m_cl_merge_hash_table.end();itr++){
        bool success = 0;
        if(!itr->first->first.CheckEmptyObject() && !itr->first->second.CheckEmptyObject()){
            //只要M-H或M-2H一个满足时间范围即可
            if((itr->first->first.GetSampleRt() > min_rt && itr->first->first.GetSampleRt() < max_rt) || (itr->first->second.GetSampleRt() > min_rt && itr->first->second.GetSampleRt() < max_rt)){
                this->m_row_map_to_cl_splice_result.insert({row , itr->second});//把row和pair<Cl,Cl>*加入到哈希表
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }
        }
        //如果只存在M-H，不存在M-2H
        else if(!itr->first->first.CheckEmptyObject()){
            if(itr->first->first.GetSampleRt() > min_rt && itr->first->first.GetSampleRt() < max_rt){
                this->m_row_map_to_cl_splice_result.insert({row , itr->second});//把row和pair<Cl,Cl>*加入到哈希表
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                success = 1;
            }
        }
        //如果只存在M-2H，不存在M-H
        else if(!itr->first->second.CheckEmptyObject()){
            if(itr->first->second.GetSampleRt() > min_rt && itr->first->second.GetSampleRt() < max_rt){
                this->m_row_map_to_cl_splice_result.insert({row , itr->second});
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }

        }
        if(success){
            row++;
        }
    }//显示MLCL
    for(auto itr = this->m_mlcl_merge_hash_table.begin() ; itr != this->m_mlcl_merge_hash_table.end();itr++){
        bool success = 0;
        if(!itr->first->first.CheckEmptyObject() && !itr->first->second.CheckEmptyObject()){
            //只要M-H或M-2H一个满足时间范围即可
            if((itr->first->first.GetSampleRt() > min_rt && itr->first->first.GetSampleRt() < max_rt) || (itr->first->second.GetSampleRt() > min_rt && itr->first->second.GetSampleRt() < max_rt)){
                this->m_row_map_to_mlcl_splice_result.insert({row , itr->second});//把row和pair<Cl,Cl>*加入到哈希表
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }
        }
        else if(!itr->first->first.CheckEmptyObject()){
            if(itr->first->first.GetSampleRt() > min_rt && itr->first->first.GetSampleRt() < max_rt){
                this->m_row_map_to_mlcl_splice_result.insert({row , itr->second});
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                success = 1;
            }
        }
        else if(!itr->first->second.CheckEmptyObject()){
            if(itr->first->second.GetSampleRt() > min_rt && itr->first->second.GetSampleRt() < max_rt){
                this->m_row_map_to_mlcl_splice_result.insert({row , itr->second});
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }

        }
        if(success){
            row++;
        }
    }
    //显示DLCL
    for(auto itr = this->m_dlcl_merge_hash_table.begin() ; itr != this->m_dlcl_merge_hash_table.end();itr++){
        bool success = 0;
        if(!itr->first->first.CheckEmptyObject() && !itr->first->second.CheckEmptyObject()){
            //只要M-H或M-2H一个满足时间范围即可
            if((itr->first->first.GetSampleRt() > min_rt && itr->first->first.GetSampleRt() < max_rt) || (itr->first->second.GetSampleRt() > min_rt && itr->first->second.GetSampleRt() < max_rt)){
                this->m_row_map_to_dlcl_splice_result.insert({row , itr->second});//把row和pair<Cl,Cl>*加入到哈希表
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }
        }
        else if(!itr->first->first.CheckEmptyObject()){
            if(itr->first->first.GetSampleRt() > min_rt && itr->first->first.GetSampleRt() < max_rt){
                this->m_row_map_to_dlcl_splice_result.insert({row , itr->second});
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first->first.GetChainLength()) + ":" + QString::number(itr->first->first.GetUnsaturation()) + ":" + QString::number(itr->first->first.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                success = 1;
            }
        }
        else if(!itr->first->second.CheckEmptyObject()){
            if(itr->first->second.GetSampleRt() > min_rt && itr->first->second.GetSampleRt() < max_rt){
                this->m_row_map_to_dlcl_splice_result.insert({row , itr->second});//加入哈希表
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->first->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->first->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->first->second.GetChainLength()) + ":" + QString::number(itr->first->second.GetUnsaturation()) + ":" + QString::number(itr->first->second.GetOxygen()));
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }

        }
        if(success){
            row++;
        }
    }
    qtablewidget->setRowCount(row);//设置为真正有多少行
}

void FragmentCombiner::PrintSpliceResultByRow(unsigned int row)
{
    auto cl_match_itr = this->m_row_map_to_cl_splice_result.find(row);
    auto mlcl_match_itr = this->m_row_map_to_mlcl_splice_result.find(row);
    auto dlcl_match_itr = this->m_row_map_to_dlcl_splice_result.find(row);
    //如果这个行是CL的
    if(cl_match_itr != this->m_row_map_to_cl_splice_result.end()){
        for(auto itr = cl_match_itr->second->begin() ; itr != cl_match_itr->second->end() ; itr++){
            qDebug() << itr->ShowInfo();//打印信息
        }
        qDebug() << "-------------------------------";
        return;
    }

    if(mlcl_match_itr != this->m_row_map_to_mlcl_splice_result.end()){
        for(auto itr = mlcl_match_itr->second->begin() ; itr != mlcl_match_itr->second->end() ; itr++){
            qDebug() << itr->ShowInfo();//打印信息
        }
        qDebug() << "-------------------------------";
        return;
    }

    if(dlcl_match_itr != this->m_row_map_to_dlcl_splice_result.end()){
        for(auto itr = dlcl_match_itr->second->begin() ; itr != dlcl_match_itr->second->end() ; itr++){
            qDebug() << itr->ShowInfo();//打印信息
        }
        qDebug() << "-------------------------------";
        return;
    }
}





