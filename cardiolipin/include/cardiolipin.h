#ifndef CARDIOLIPIN_H
#define CARDIOLIPIN_H

#include <database.h>
#include <mzml.h>
#include <memory>
#include <cmath>
#include <headgroup.h>
#include <set>
#include <pa.h>

//一个一级只和一个配对结果进行存储，因为要用于M-H和M-2H的配对，配对是利用配对结果进行计算的，一对一地配对
class Cardiolipin
{
public:
    Cardiolipin();//默认构造函数
    Cardiolipin(Ms1* Ms1_ptr , DatabaseRecord* DatabaseRecord_ptr);
public:
    void SetMs1Ptr(Ms1* ms1_ptr);
    void SetMs1DatabaseRecordPtr(DatabaseRecord* ms1_databaseRecord_ptr);
    void SetMs1MatchingScosre(float ms1_matching_score);
    void SetMs2VectorPtr(std::vector<Ms2*> ms2_vector_ptr);
    Ms1* GetMs1Ptr();
    DatabaseRecord* GetMs1DatabaseRecordPtr();
    float GetMs1MatchingScore();
    std::vector<Ms2*> GetMs2VectorPtr();
    unsigned int GetChainLength();
    unsigned int GetUnsaturation();
    unsigned int GetOxygen();
    void ClearMs2Info();

    QString ShowInfo();//展示信息
    std::array<unsigned int , 3> GetMs1DatabaseRecordComponent();//获得一级配对的chain length:unsaturation:oxygen
    QString GetMs1DatabaseRecordComponentWithQString();
    float GetSampleMz();
    float GetSampleRt();
    float GetSampleIntensity();
private:
    Ms1* m_Ms1_ptr;
    DatabaseRecord* m_Ms1_DatabaseRecord_ptr;//记录的一级配对结果
    float m_Ms1_Matching_Score;
    std::vector<Ms2*> m_Ms2_vector_ptr;//配对的二级结果
public:
    void ScoreMs1(float& k);//常数k

    //一级找对应的二级
public:
    void MatchCardiolipinWithMs2(std::vector<Ms2>& ms2_vector , float& dalton , float& tolerance_rt_scope);//与二级进行配对，std::vector<Ms2>需要已经根据precursor_ion_mz排序
    bool CheckMs2Exist();//检查一级配对完二级后，是否匹配上二级
    bool CheckEmptyObject();//检查这个对象是否为空的对象
    virtual void EmptyObject();//把对象清空

    //二级找头基
public:
    void CheckHeadgroup(float& ppm , float& mz_score_weight , float& k);//常数k，用于计算mz分数

    //二级找PA和FA
public:
    void FindPaAndFa(float& ppm , float& mz_score_weight , float& k , Database& database);//常数k，用于计算mz分数
public:
    virtual void splice();
};

#endif // SINGLEMS1MATCH_H
