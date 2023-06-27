#ifndef CL_H
#define CL_H
#include <array>
#include <list>
#include <unordered_map>
#include "ClSpecificStructure.h"
#include <cardiolipin.h>
class Cl: public Cardiolipin
{
public:
    Cl();
    Cl(Cardiolipin& superclass);
public:
    std::list<ClSpecificStructure> m_cl_specific_structure_vector;//存储所有二级的拼接结果
public:
    void splice() override;//覆写父类的splice函数
    void MergeSplice();//合并多个二级的拼接结果，只要四个FA的链长，不饱和度，氧个数相同，而不需要加和形式都相同，才认为是相同的。
    void FourFaSpliceCl(Ms2* ms2_ptr);//用4个FA来拼接成这个Cl
    void TwoPaSpliceCl( Ms2* ms2_ptr);//用2个PA来拼接这个Cl
    void TwoPaFourFaSpliceCl(Ms2* ms2_ptr);//两个PA和四个FA拼接这个Cl
    std::set<std::set<Fa*>>& FourFaSpliceTwoPa(Pa* pa_1_ptr , Pa* pa_2_ptr , std::vector<Fa>* fa_vector_ptr , std::set<std::set<Fa*>>& store , Ms2* ms2_ptr);//四个FA拼接两个PA，返回引用
public:
    void EmptyObject() override;
};




#endif // CL_H
