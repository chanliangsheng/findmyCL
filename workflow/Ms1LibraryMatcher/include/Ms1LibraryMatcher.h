#ifndef MS1LIBRARYMATCH_H
#define MS1LIBRARYMATCH_H

#include <utility>
#include <database.h>
#include <map>
#include <QTime>
#include <workflow.h>
#include <algorithm>
#include <QObject>
#include <QListWidget>
#include <QTableWidget>

class Ms1LibraryMatcher : public Workflow
{
    Q_OBJECT
public:
    Ms1LibraryMatcher();
    Ms1LibraryMatcher(float ppm , float tolerance_rt , float ppm_with_half_score);
public:
    static float m_ppm;//ppm，初始化为5
    static float m_tolerance_rt_scope;//可容忍的时间，初始化为6s

    static float m_ppm_with_half_score;//设定ppm为多少的时候分数为0.5，初始值为5
    static float m_k;//m_ppm_with_half_score对应的常数项
public:
    //对CL，MLCL，DLCL进行搜索
    void MatchMs1WithAllTables(Mzml& mzml , Database& database);
private:
    //一级与[M-H]-和[M-2H]2-的库进行配对，返回m_cl_match / m_mlcl_match / m_dlcl_match的智能指针
    std::shared_ptr<std::vector<std::pair<Cardiolipin , Cardiolipin>>> MatchMs1With2Tables(std::vector<Ms1>& ms1_vector , std::pair<std::vector<DatabaseRecord>* , std::vector<DatabaseRecord>*> database_record_vector_pair);
    //寻找相互匹配的M-H和M-2H
    void MatchM_hWithM_2h(std::multimap<int , Ms1*>& m_h_hash_map , std::multimap<int , Ms1*>& m_h_left_hash_map ,int key , Ms1* ms1_pointer ,std::vector<DatabaseRecord>* m_h_database_record , std::vector<DatabaseRecord>* m_2h_database_record , std::shared_ptr<std::vector<std::pair<Cardiolipin , Cardiolipin>>>& ret);
    //for QListWidget

};

#endif // MS1MATCH_H
