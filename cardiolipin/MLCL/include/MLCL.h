#ifndef MLCL_H
#define MLCL_H

#include <cardiolipin.h>
#include "MlclSpecificStructure.h"

class Mlcl:public Cardiolipin
{
public:
    Mlcl();
    Mlcl(Cardiolipin& superclass);
public:
    std::list<MlclSpecificStructure> m_mlcl_specific_structure_vector;//存储所有二级的拼接结果
public:
    void splice() override;//覆写父类的splice函数
    void MergeSplice();//合并多个二级的拼接结果
    void ThreeFaSpliceMlcl(Ms2* ms2_ptr);//用3个FA来拼接成这个MLCl
    void OnePaThreeFaSpliceMlcl(Ms2* ms2_ptr);//1个PA和3个FA拼接这个Mlcl
public:
    void EmptyObject() override;
};

#endif // MLCL_H
