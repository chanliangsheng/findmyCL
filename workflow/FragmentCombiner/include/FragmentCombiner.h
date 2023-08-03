#ifndef FRAGMENTCOMBINER_H
#define FRAGMENTCOMBINER_H

#include <workflow.h>
#include <FragmentFinder.h>
#include <QFile>


class FragmentCombiner:public Workflow
{
public:
    FragmentCombiner();
    FragmentCombiner(FragmentFinder& fragment_finder , std::string mode);
public:
    void CopyInfoFromMs2WithPaAndFa(FragmentFinder& fragment_finder);
public:
    void splice();
    void MergeClPair(std::pair<Cl , Cl>& cl_pair);
    void MergeMlclPair(std::pair<Mlcl , Mlcl>& mlcl_pair);
    void MergeDlclPair(std::pair<Dlcl , Dlcl>& dlcl_pair);
public:
    void OutputResultWithTxt(QString folder_path);//输出结果的txt文件
    void OutPutWithCsv(QString folder_path);//输出结果的csv文件
    void Filter();//过滤某些心磷脂
    void ClearSpliceResult();
public:
    static std::string mode;//拼接的模式，初始化为strict，还有另外一个模式flexible
public:
    //M-H和M-2H的合并结果
    std::map<std::pair<Cl , Cl>* , std::list<ClSpecificStructure>*> m_cl_merge_hash_table;
    std::map<std::pair<Mlcl , Mlcl>* , std::list<MlclSpecificStructure>*> m_mlcl_merge_hash_table;
    std::map<std::pair<Dlcl , Dlcl>* , std::list<DlclSpecificStructure>*> m_dlcl_merge_hash_table;

public slots:
    void ShowSpliceResult(QTableWidget* qtablewidget);
    void ShowSpliceResult(QTableWidget* qtablewidget , float min_rt , float max_rt);
    void PrintSpliceResultByRow(unsigned int row);
public:
    std::map<unsigned int , std::list<ClSpecificStructure>*> m_row_map_to_cl_splice_result;
    std::map<unsigned int , std::list<MlclSpecificStructure>*> m_row_map_to_mlcl_splice_result;
    std::map<unsigned int , std::list<DlclSpecificStructure>*> m_row_map_to_dlcl_splice_result;//存储row和数据的对应关系
};


#endif // MS1MATCH_H
