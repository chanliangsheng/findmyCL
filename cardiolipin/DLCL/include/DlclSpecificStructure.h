#ifndef DLCLSPECIFICSTRUCTURE_H
#define DLCLSPECIFICSTRUCTURE_H

#include <cardiolipin.h>
#include <map>
#include <PaNode.h>
#include <unordered_map>

//dlcl根节点
class DlclSpecificStructure
{
public:
    DlclSpecificStructure();
    DlclSpecificStructure(Fa* fa_1_ptr, Fa* fa_2_ptr , Ms2* ms2_ptr);
    DlclSpecificStructure(Pa *pa_1_ptr , Fa* fa_1_ptr, Fa* fa_2_ptr , Ms2* ms2_ptr);
    DlclSpecificStructure(Pa *pa_1_ptr , Ms2* ms2_ptr);
public:
    Ms2* m_ms2;
    PaNode* m_left_pa_ptr;
    float m_score;
public:
    bool m_pa_exist = 0;
    bool m_fa_exist = 0;
    float m_total_intensity = 0;
    std::set<std::array<unsigned int,3>> m_pa_info;
    std::set<std::array<unsigned int,3>> m_fa_info;
public:
    void score();//对这个拼接结果进行打分
    static float m_fragment_score_weight;//碎片分数的权重
    static float m_fa_consistency_score_weight;//FA一致性的权重
    static float m_pa_exist_score_weight;//PA存在的权重
    static float m_fa_intensity_variance_score_weight;//FA的强度的方差分数的权重
public:
    float GetTotalIntensity();
    void CalculateTotalIntensity();
    float GetMs2TotalIntensity();
public:
    bool operator==(const DlclSpecificStructure& other);//判定两个MlCl是否相同，如何相同，则改变this的Mlcl，如果不同，则不进行改变，并返回false；声明为非常量函数，意味着可以修改里面的值
    std::shared_ptr<DlclSpecificStructure> merge(DlclSpecificStructure* other , bool merge_m_h_m_2h = 0);//合并两个Cl的信息，但是不改变某一个Cl的信息，会返回一个新的Cl的共享指针
    void update();
public:
     //可用于M-H和M-2H的合并；也可用于Dlcl的多个二级合并中使用
    bool StrictMerge(DlclSpecificStructure& other);//用于M-H和M-2H的严格模式合并，与==不同的是，分数的最终计算方式不同
public:
    bool FlexibleMerge(DlclSpecificStructure& other);
    std::vector<DlclSpecificStructure*> m_flexible_mode_dif_pa_merge;//存储在flexible模式下，如果两都有PA和FA，但是PA不相同，则把另外的PA存储到m_flexible_mode_merge中

public:
    QString ShowInfo();
    QString ShowSimpleInfo();
public:
public:
    std::vector<int> GetTotalInfo();
    int GetTotalChainLength();
    int GetTotalUnsaturation();
    int GetTotalOxygen();
};


#endif // DLCLSPECIFICSTRUCTURE_H
