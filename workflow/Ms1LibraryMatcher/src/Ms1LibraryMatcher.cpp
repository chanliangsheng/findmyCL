#include "Ms1LibraryMatcher.h"
#include <map>

float Ms1LibraryMatcher::m_ppm = 5;
float Ms1LibraryMatcher::m_tolerance_rt_scope = 6;
float Ms1LibraryMatcher::m_ppm_with_half_score = 5;
float Ms1LibraryMatcher::m_k = -((1.17741/pow((Ms1LibraryMatcher::m_ppm_with_half_score)/static_cast<double>(1000000), 2))/2);

using namespace std;
Ms1LibraryMatcher::Ms1LibraryMatcher()
{

}

Ms1LibraryMatcher::Ms1LibraryMatcher(float ppm, float tolerance_rt , float ppm_with_half_score)
{
    this->m_ppm = ppm;
    this->m_tolerance_rt_scope = tolerance_rt;
    this->m_k = -((1.17741/pow((ppm_with_half_score)/static_cast<double>(1000000), 2))/2);
}

std::shared_ptr<std::vector<std::pair<Cardiolipin , Cardiolipin>>> Ms1LibraryMatcher::MatchMs1With2Tables(std::vector<Ms1>& ms1_vector , std::pair<std::vector<DatabaseRecord>* , std::vector<DatabaseRecord>*> database_record_vector_pair)
{
    //std::vector<DatabaseRecord>中的mz已经按照从小到大排序了，所以不需要再排序
    //M-H和M-2H的同一个索引下的心磷脂是相同的，只是加和形式的区别
    vector<DatabaseRecord>* m_h_database_record = database_record_vector_pair.first;
    vector<DatabaseRecord>* m_2h_database_record = database_record_vector_pair.second;



    //结果ret，智能指针
    shared_ptr<vector<pair<Cardiolipin,Cardiolipin>>> ret = make_shared<vector<pair<Cardiolipin,Cardiolipin>>>();



    multimap<int , Ms1*> m_h_hash_map;//M-H配对上的结果，多对应哈希表，因为同一个int可能对应多个MS1*
    multimap<int , Ms1*> m_h_left_hash_map;//用来存储仅有M-H，没有M-2H的M-H
    //二分法，先搜索M-H，再搜索M-2H，然后去除同时再M-H的部分


    //对M-H进行搜索
    for(std::vector<Ms1>::iterator ms1_itr = ms1_vector.begin() ; ms1_itr != ms1_vector.end() ; ms1_itr++){
        //左索引和右索引
        int left = 0;
        int right = m_h_database_record->size() - 1;
//        if(ms1_itr->GetMz() > 1521.978 && ms1_itr->GetMz() < 1522 && ms1_itr->m_rt < 815.19 && ms1_itr->m_rt > 815.1){
//            qDebug() << "stop";
//        }
        while (left <= right) {
            int mid = (left + right) / 2;
            //基于ppm所给定的范围，需要样本中的mz位于m_h_min_mz和m_h_max_mz，才认为是符合的
            float m_h_min_mz = m_h_database_record->at(mid).GetMz() - (m_h_database_record->at(mid).GetMz() * (this->m_ppm)/1000000);
            float m_h_max_mz = m_h_database_record->at(mid).GetMz() + (m_h_database_record->at(mid).GetMz() * (this->m_ppm)/1000000);
            //如果符合要求
            if((ms1_itr->GetMz() >= m_h_min_mz) && (ms1_itr->GetMz() <= m_h_max_mz)){
                Ms1* ms1_pointer = &(*ms1_itr);//迭代器转换为指针
                m_h_hash_map.insert({mid , ms1_pointer});//加入到哈希表中
                //mid的左边一位和右边一位的各自的符合区间
                int left_t = mid - 1;
                int right_t = mid + 1;

                float m_h_left_t_min_mz;
                float m_h_left_t_max_mz;
                float m_h_right_t_min_mz;
                float m_h_right_t_max_mz;
                while(left_t >= 0){
                    //更新区间
                    m_h_left_t_min_mz = m_h_database_record->at(left_t).GetMz() - (m_h_database_record->at(left_t).GetMz() * (this->m_ppm)/1000000);
                    m_h_left_t_max_mz = m_h_database_record->at(left_t).GetMz() + (m_h_database_record->at(left_t).GetMz() * (this->m_ppm)/1000000);
                    if((ms1_itr->GetMz() >= m_h_left_t_min_mz) && (ms1_itr->GetMz() <= m_h_left_t_max_mz)){
                        m_h_hash_map.insert({left_t , ms1_pointer});//加入到哈希表中
                        left_t--;//向左继续找
                    }
                    else{
                        break;
                    }
                }
                while(right_t <= (m_h_database_record->size() - 1)){
                    //更新区间
                    m_h_right_t_min_mz = m_h_database_record->at(right_t).GetMz() - (m_h_database_record->at(right_t).GetMz() * (this->m_ppm)/1000000);
                    m_h_right_t_max_mz = m_h_database_record->at(right_t).GetMz() + (m_h_database_record->at(right_t).GetMz() * (this->m_ppm)/1000000);
                    if((ms1_itr->GetMz() >= m_h_right_t_min_mz) && (ms1_itr->GetMz() <= m_h_right_t_max_mz)){
                        m_h_hash_map.insert({right_t , ms1_pointer});
                        right_t++;//向右继续寻找
                    }
                    else{
                        break;
                    }
                }

                break;//找完左边和右边的，就说明ms1_itr的元素已经找完了，跳出这个while
            }
            //如果mz大于mid对应的mz，说明在右侧可能有配对的
            else if(ms1_itr->GetMz() >= m_h_max_mz){
                left = mid + 1;
            }
            //如果mz小于mid对应的mz，说明在左侧可能有配对的
            else if(ms1_itr->GetMz() <= m_h_min_mz){
                right = mid - 1;
            }
        }
    }

    m_h_left_hash_map = m_h_hash_map;//赋值，m_h_left_hash_map用来保留那些仅有M-H的


    //对M-2H进行搜索，查找既有M-H又有M-2H的一级，m_h_hash_map去除了既有M-H中有M-2H的部分
    for(vector<Ms1>::iterator ms1_itr = ms1_vector.begin() ; ms1_itr != ms1_vector.end() ; ms1_itr++){
        //左索引和右索引
        int left = 0;
        int right = m_2h_database_record->size() - 1;
        int mid;
        while (left <= right) {
            mid = (left + right)/2;//中值
            //基于ppm所给定的范围，需要样本中的mz位于m_2h_min_mz和m_2h_max_mz，才认为是符合的
            float m_2h_min_mz = m_2h_database_record->at(mid).GetMz() - (m_2h_database_record->at(mid).GetMz() * (this->m_ppm)/1000000);
            float m_2h_max_mz = m_2h_database_record->at(mid).GetMz() + (m_2h_database_record->at(mid).GetMz() * (this->m_ppm)/1000000);
            Ms1* ms1_pointer;//一级的指针


            if((ms1_itr->GetMz() >= m_2h_min_mz) && (ms1_itr->GetMz() <= m_2h_max_mz)){
                ms1_pointer = &*ms1_itr;
                //判断哪些M-H跟M-2H是一起的，加入到结果ret中，ret传的是引用，所以会直接在ret上修改
                this->MatchM_hWithM_2h(m_h_hash_map , m_h_left_hash_map , mid , ms1_pointer , m_h_database_record , m_2h_database_record , ret);
                int left_t = mid - 1;
                int right_t = mid + 1;
                //mid的左边一位和右边一位的各自的符合区间
                float m_2h_left_t_min_mz;
                float m_2h_left_t_max_mz;
                float m_2h_right_t_min_mz;
                float m_2h_right_t_max_mz;
                while(left_t >= 0){
                    m_2h_left_t_min_mz = m_2h_database_record->at(left_t).GetMz() - (m_2h_database_record->at(left_t).GetMz() * (this->m_ppm)/1000000);
                    m_2h_left_t_max_mz = m_2h_database_record->at(left_t).GetMz() + (m_2h_database_record->at(left_t).GetMz() * (this->m_ppm)/1000000);
                    if((ms1_itr->GetMz() >= m_2h_left_t_min_mz) && (ms1_itr->GetMz() <= m_2h_left_t_max_mz)){
                        this->MatchM_hWithM_2h(m_h_hash_map , m_h_left_hash_map , left_t , ms1_pointer , m_h_database_record , m_2h_database_record , ret);
                        left_t--;//向左继续找
                    }
                    else{
                        break;
                    }
                }

                while (right_t <= m_2h_database_record->size() - 1) {
                    m_2h_right_t_min_mz = m_2h_database_record->at(right_t).GetMz() - (m_2h_database_record->at(right_t).GetMz() * (this->m_ppm)/1000000);
                    m_2h_right_t_max_mz = m_2h_database_record->at(right_t).GetMz() + (m_2h_database_record->at(right_t).GetMz() * (this->m_ppm)/1000000);
                    if((ms1_itr->GetMz() >= m_2h_right_t_min_mz) && (ms1_itr->GetMz() <= m_2h_right_t_max_mz)){
                        this->MatchM_hWithM_2h(m_h_hash_map , m_h_left_hash_map , right_t , ms1_pointer , m_h_database_record , m_2h_database_record , ret);
                        right_t++;
                    }
                    else{
                        break;
                    }
                }

                break;
            }
            else if(ms1_itr->GetMz() >= m_2h_min_mz){
                left = mid + 1;
            }
            else if(ms1_itr->GetMz() <= m_2h_max_mz){
                right = mid - 1;
            }
        }
    }


//    for(auto itr = m_h_left_hash_map.begin() ; itr != m_h_left_hash_map.end() ; itr++){
//        if(itr->second->GetMz() > 1521.97 && itr->second->GetMz() < 1522 && itr->second->m_rt < 815.18 && itr->second->m_rt > 815.17){
//            qDebug() << "stop";
//        }
//    }

    //遍历哈希表，m_h_left_hash_map仅有M-H，其中既有M-H又有M-2H已经在Ms1LibraryMatcher::MatchM_hWithM_2h被排除了
    for(std::map<int , Ms1*>::iterator itr_map = m_h_left_hash_map.begin(); itr_map != m_h_left_hash_map.end() ; itr_map ++){
        Cardiolipin m_h(itr_map->second , &m_h_database_record->at(itr_map->first));//first指的是mz在m_h_database_record中匹配的索引，second指的是MS1*
        Cardiolipin m_2h;
        ret->emplace_back(pair<Cardiolipin,Cardiolipin>(m_h,m_2h));
    }


    //进行打分
    for(auto itr = ret->begin() ; itr!= ret->end() ; itr++){
        itr->first.ScoreMs1(this->m_k);
        itr->second.ScoreMs1(this->m_k);
    }

    return ret;
}

void Ms1LibraryMatcher::MatchMs1WithAllTables(Mzml& mzml, Database &database)
{
    vector<Ms1>& ms1_vector = mzml.GetLocalMs1Vector();


    //重新计算m_k
    this->m_k = -((1.17741/pow((this->m_ppm_with_half_score)/static_cast<double>(1000000), 2))/2);
    //清空上次配对的结果
    vector<pair<Cl , Cl>>().swap(this->m_cl_vector);
    vector<pair<Mlcl , Mlcl>>().swap(this->m_mlcl_vector);
    vector<pair<Dlcl , Dlcl>>().swap(this->m_dlcl_vector);

    //配对，返回值为std::shared_ptr<std::vector<std::pair<Cardiolipin , Cardiolipin>>>
    //MatchMs1With2Tables传入数据库的两个向量对的指针，而不是引用，引用离开函数，会结束生命周期；原本用引用指向的databaseRecord，但是结束后对应的地址会被销毁，所以需要用指针传入
    shared_ptr<vector<pair<Cardiolipin , Cardiolipin>>> cl_match_pointer = this->MatchMs1With2Tables(ms1_vector , database.GetLocalClPair());
    shared_ptr<vector<pair<Cardiolipin , Cardiolipin>>> mlcl_match_pointer = this->MatchMs1With2Tables(ms1_vector , database.GetLocalMlclPair());
    shared_ptr<vector<pair<Cardiolipin , Cardiolipin>>> dlcl_match_pointer = this->MatchMs1With2Tables(ms1_vector , database.GetLocalDlclPair());

    //重新设置大小
    this->m_cl_vector.resize(cl_match_pointer->size());
    this->m_mlcl_vector.resize(mlcl_match_pointer->size());
    this->m_dlcl_vector.resize(dlcl_match_pointer->size());

    //转换成子类
    transform(cl_match_pointer->begin() , cl_match_pointer->end() , this->m_cl_vector.begin() , [=](pair<Cardiolipin, Cardiolipin>(x)){return pair<Cl,Cl>(Cl(x.first),Cl(x.second));});
    transform(mlcl_match_pointer->begin() , mlcl_match_pointer->end() , this->m_mlcl_vector.begin() , [=](pair<Cardiolipin, Cardiolipin>(x)){return pair<Mlcl , Mlcl>(Mlcl(x.first),Mlcl(x.second));});
    transform(dlcl_match_pointer->begin() , dlcl_match_pointer->end() , this->m_dlcl_vector.begin() , [=](pair<Cardiolipin, Cardiolipin>(x)){return pair<Dlcl , Dlcl>(Dlcl(x.first),Dlcl(x.second));});
}



void Ms1LibraryMatcher::MatchM_hWithM_2h(std::multimap<int , Ms1*>& m_h_hash_map , std::multimap<int , Ms1*>& m_h_left_hash_map ,int key , Ms1* ms1_pointer ,std::vector<DatabaseRecord>* m_h_database_record , std::vector<DatabaseRecord>* m_2h_database_record , std::shared_ptr<std::vector<std::pair<Cardiolipin,Cardiolipin>>>& ret)
{
    auto m_h_hash_map_pair_itr =  m_h_hash_map.equal_range(key);//寻找M-H中是否匹配上同一个心磷脂，M-H和M-2H是同一个mid，多重哈希表返回第一个对应位置的迭代器
    //如果在M-H中没有找到同一个心磷脂，则只将m_2h加入到ret中
    if(m_h_hash_map_pair_itr.first == m_h_hash_map_pair_itr.second){
        Cardiolipin m_h;
        Cardiolipin m_2h(ms1_pointer , &m_2h_database_record->at(key));
        ret->emplace_back(pair<Cardiolipin,Cardiolipin>(m_h , m_2h));
    }
    //如果M-H找到了同一个心磷脂
    else{
        // 追加到结果中
        bool merge_success = 0;
        for(auto itr = m_h_hash_map_pair_itr.first ; itr != m_h_hash_map_pair_itr.second ; itr++){
            if(fabs(itr->second->GetRt() - ms1_pointer->GetRt()) <= this->m_tolerance_rt_scope){
                merge_success = 1;
                Cardiolipin m_h(itr->second , &m_h_database_record->at(key));
                Cardiolipin m_2h(ms1_pointer , &m_2h_database_record->at(key));
                ret->emplace_back(pair<Cardiolipin,Cardiolipin>(m_h , m_2h));
            }
        }
        //如果没有合并的结果，才把没有M-H的结果加入到ret中
        if(!merge_success){
            Cardiolipin m_h;
            Cardiolipin m_2h(ms1_pointer , &m_2h_database_record->at(key));
            ret->emplace_back(pair<Cardiolipin,Cardiolipin>(m_h , m_2h));
        }


        //剩下的部分寻找相同项
        auto m_h_left_hash_map_pair_itr = m_h_left_hash_map.equal_range(key);
        //如果没有找到
        if(m_h_left_hash_map_pair_itr.first == m_h_left_hash_map_pair_itr.second){

        }
        else{
            for(auto itr = m_h_left_hash_map_pair_itr.first ; itr != m_h_left_hash_map_pair_itr.second ;){
                if(fabs(itr->second->GetRt() - ms1_pointer->GetRt()) <= this->m_tolerance_rt_scope){
                    itr = m_h_left_hash_map.erase(itr);
                }
                else{
                    itr++;
                }
            }
        }
    }
}
