#ifndef MSLEVELMATCHER_H
#define MSLEVELMATCHER_H

#include <workflow.h>
#include <Ms1LibraryMatcher.h>
#include <mzml.h>
#include <algorithm>

//继承自Workflow
class MsLevelMatcher:public Workflow
{
public:
    MsLevelMatcher();
    MsLevelMatcher(Ms1LibraryMatcher & ms1_library_matcher , float dalton);
    MsLevelMatcher(Ms1LibraryMatcher & ms1_library_matcher);// 把ms1_match中的结果复制过来，不采用指针的方法是为了每一步在在执行之后，如果再执行，用的还是上一步的数据；所以workflow中没有使用指针
public:
    void CopyInfoFromMs1Match(Ms1LibraryMatcher & ms1_library_matcher);
    void CopyInfoFromMs1Match(Ms1LibraryMatcher & ms1_library_matcher , float dalton);
public:
    static float m_dalton;//初始化为0.5
    static float m_tolerance_rt_scope;//初始化为8s
//public:
//    std::vector<std::pair<Cl , Cl>> m_cl_vector;
//    std::vector<std::pair<Mlcl , Mlcl>> m_mlcl_vector;
//    std::vector<std::pair<Dlcl , Dlcl>> m_dlcl_vector;
private:
    //无二级的结果
    std::vector<std::pair<Cl , Cl>> m_cl_withoutMs2_vector;
    std::vector<std::pair<Mlcl , Mlcl>> m_mlcl_withoutMs2_vector;
    std::vector<std::pair<Dlcl , Dlcl>> m_dlcl_withoutMs2_vector;
public:
    void MatchCardiolipinWithMs2(Mzml& mzml);
    void OutPutWithoutMS2WithCsv(QString folder_path);//输出结果的csv文件
};

#endif // MS1MATCH_H
