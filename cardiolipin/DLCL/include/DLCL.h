#ifndef DLCL_H
#define DLCL_H

#include <cardiolipin.h>
#include <map>
#include "DlclSpecificStructure.h"

class Dlcl:public Cardiolipin
{
public:
    Dlcl();
    Dlcl(Cardiolipin& superclass);
public:
    std::list<DlclSpecificStructure> m_dlcl_specific_structure_vector;//存储所有二级的拼接结果
public:
    void splice() override;//覆写父类的splice函数
    void MergeSplice();//合并多个二级的拼接结果
    void TwoFaSpliceDlcl(Ms2* ms2_ptr);//用2个FA来拼接成这个MLCl
    void OnePaTwoFaSpliceDlcl(Ms2* ms2_ptr);//1个PA和3个FA拼接这个Mlcl
public:
    void EmptyObject() override;
    void ClearSpliceResult();
    void DeleteRedundantSpliceResult();
};


#endif // DLCL_H
