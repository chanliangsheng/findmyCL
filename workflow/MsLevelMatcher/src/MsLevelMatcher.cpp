#include <MsLevelMatcher.h>

float MsLevelMatcher::m_dalton = 0.5;
float MsLevelMatcher::m_tolerance_rt_scope = 8;

using namespace std;
MsLevelMatcher::MsLevelMatcher()
{
    this->m_cl_vector = std::vector<std::pair<Cl , Cl>>();
    this->m_mlcl_vector = std::vector<std::pair<Mlcl , Mlcl>>();
    this->m_dlcl_vector = std::vector<std::pair<Dlcl , Dlcl>>();

    this->m_cl_withoutMs2_vector = std::vector<std::pair<Cl , Cl>>();
    this->m_mlcl_withoutMs2_vector = std::vector<std::pair<Mlcl , Mlcl>>();
    this->m_dlcl_withoutMs2_vector = std::vector<std::pair<Dlcl , Dlcl>>();
}

MsLevelMatcher::MsLevelMatcher(Ms1LibraryMatcher & ms1_library_matcher, float dalton = 0.5)
{
    // 把ms1_library_matcher中的结果复制过来，不采用指针的方法是为了每一步在在执行之后，如果再执行，用的还是上一步的数据；所以workflow中没有使用指针
    this->m_cl_vector = ms1_library_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms1_library_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms1_library_matcher.GetCopyDlclPairVector();
    this->m_dalton = dalton;

    this->m_cl_withoutMs2_vector = std::vector<std::pair<Cl , Cl>>();
    this->m_mlcl_withoutMs2_vector = std::vector<std::pair<Mlcl , Mlcl>>();
    this->m_dlcl_withoutMs2_vector = std::vector<std::pair<Dlcl , Dlcl>>();
}

MsLevelMatcher::MsLevelMatcher(Ms1LibraryMatcher & ms1_library_matcher)
{
    // 把ms1_library_matcher中的结果复制过来，不采用指针的方法是为了每一步在在执行之后，如果再执行，用的还是上一步的数据；所以workflow中没有使用指针
    this->m_cl_vector = ms1_library_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms1_library_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms1_library_matcher.GetCopyDlclPairVector();

    this->m_cl_withoutMs2_vector = std::vector<std::pair<Cl , Cl>>();
    this->m_mlcl_withoutMs2_vector = std::vector<std::pair<Mlcl , Mlcl>>();
    this->m_dlcl_withoutMs2_vector = std::vector<std::pair<Dlcl , Dlcl>>();
}

void MsLevelMatcher::CopyInfoFromMs1Match(Ms1LibraryMatcher & ms1_library_matcher)
{
    // 把ms1_library_matcher中的结果复制过来，不采用指针的方法是为了每一步在在执行之后，如果再执行，用的还是上一步的数据；所以workflow中没有使用指针
    this->m_cl_vector = ms1_library_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms1_library_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms1_library_matcher.GetCopyDlclPairVector();

    this->m_cl_withoutMs2_vector = std::vector<std::pair<Cl , Cl>>();
    this->m_mlcl_withoutMs2_vector = std::vector<std::pair<Mlcl , Mlcl>>();
    this->m_dlcl_withoutMs2_vector = std::vector<std::pair<Dlcl , Dlcl>>();
}

void MsLevelMatcher::CopyInfoFromMs1Match(Ms1LibraryMatcher & ms1_library_matcher, float dalton)
{
    // 把ms1_library_matcher中的结果复制过来，不采用指针的方法是为了每一步在在执行之后，如果再执行，用的还是上一步的数据；所以workflow中没有使用指针
    this->m_cl_vector = ms1_library_matcher.GetCopyClPairVector();
    this->m_mlcl_vector = ms1_library_matcher.GetCopyMlclPairVector();
    this->m_dlcl_vector = ms1_library_matcher.GetCopyDlclPairVector();
    this->m_dalton = dalton;

    this->m_cl_withoutMs2_vector = std::vector<std::pair<Cl , Cl>>();
    this->m_mlcl_withoutMs2_vector = std::vector<std::pair<Mlcl , Mlcl>>();
    this->m_dlcl_withoutMs2_vector = std::vector<std::pair<Dlcl , Dlcl>>();
}

void MsLevelMatcher::MatchCardiolipinWithMs2(Mzml& mzml)
{
    //对二级vector进行排序
    mzml.SortMs2VectorByPrecursorIonMz();
    vector<Ms2>& ms2_vector = mzml.GetLocalMs2Vector();

    for(auto itr = this->m_cl_vector.begin() ; itr != this->m_cl_vector.end();){
        itr->first.MatchCardiolipinWithMs2(ms2_vector , this->m_dalton , this->m_tolerance_rt_scope);
        itr->second.MatchCardiolipinWithMs2(ms2_vector , this->m_dalton , this->m_tolerance_rt_scope);


        //如果M-H和M-2H都没有找到二级
        if(!itr->first.CheckMs2Exist() && !itr->second.CheckMs2Exist()){
            this->m_cl_withoutMs2_vector.push_back(*itr);//把结果加入到无二级的Cl向量中
            itr = this->m_cl_vector.erase(itr);
        }
        else{
             itr++;
        }
    }

    //    std::vector<std::pair<Mlcl , Mlcl>> m_mlcl_vector;
    for(auto itr = this->m_mlcl_vector.begin() ; itr != this->m_mlcl_vector.end();){
        itr->first.MatchCardiolipinWithMs2(ms2_vector , this->m_dalton , this->m_tolerance_rt_scope);
        itr->second.MatchCardiolipinWithMs2(ms2_vector , this->m_dalton , this->m_tolerance_rt_scope);
        //如果M-H和M-2H都没有找到二级
        if(!itr->first.CheckMs2Exist() && !itr->second.CheckMs2Exist()){
            this->m_mlcl_withoutMs2_vector.push_back(*itr);//把结果加入到无二级的Cl向量中
            itr = this->m_mlcl_vector.erase(itr);
        }
        else{
             itr++;
        }
    }
    //    std::vector<std::pair<Dlcl , Dlcl>> m_dlcl_vector;
    for(auto itr = this->m_dlcl_vector.begin() ; itr != this->m_dlcl_vector.end();){
        itr->first.MatchCardiolipinWithMs2(ms2_vector , this->m_dalton , this->m_tolerance_rt_scope);
        itr->second.MatchCardiolipinWithMs2(ms2_vector , this->m_dalton , this->m_tolerance_rt_scope);
        //如果M-H和M-2H都没有找到二级
        if(!itr->first.CheckMs2Exist() && !itr->second.CheckMs2Exist()){
            this->m_dlcl_withoutMs2_vector.push_back(*itr);//把结果加入到无二级的Cl向量中
            itr = this->m_dlcl_vector.erase(itr);
        }
        else{
             itr++;
        }
    }
    this->DeleteRedundantPair();//删除冗余心磷脂
}

void MsLevelMatcher::OutPutWithoutMS2WithCsv(QString folder_path)
{
    QString fileName = folder_path + "/" + "计算机定性中无二级部分.csv";
    //QString fileName = "D:/qt_project/R_analyze/计算机定性结果/dda/计算机所有定性结果.csv";
    QFile file(fileName);
    file.open(QIODevice::WriteOnly | QIODevice::Text);//打开文件

    QTextStream stream(&file);//文本流
//    stream << "Hello, world!" << endl;
//    stream << "This is a sample text." << endl;
    //输出CL
    for(auto itr = this->m_cl_withoutMs2_vector.begin() ; itr != this->m_cl_withoutMs2_vector.end();itr++){
        QString component_message = GetPairComponentWithQString<Cl>(*itr);//获得pair的化合物信息
        pair<float,float> sample_mz = GetPairSampleMz<Cl>(*itr);
        float sample_rt = GetPairRt<Cl>(*itr);
        stream << component_message << "," << sample_mz.first << "," << sample_mz.second << "," << sample_rt << "," << "CL" << endl;
    }
    //MLCL
    for(auto itr = this->m_mlcl_withoutMs2_vector.begin() ; itr != this->m_mlcl_withoutMs2_vector.end();itr++){
        QString component_message = GetPairComponentWithQString<Mlcl>(*itr);//获得pair的化合物信息
        pair<float,float> sample_mz = GetPairSampleMz<Mlcl>(*itr);
        float sample_rt = GetPairRt<Mlcl>(*itr);
        stream << component_message << "," << sample_mz.first << "," << sample_mz.second << "," << sample_rt << "," << "MLCL" << endl;
    }
    //DLCL
    for(auto itr = this->m_dlcl_withoutMs2_vector.begin() ; itr != this->m_dlcl_withoutMs2_vector.end();itr++){
        QString component_message = GetPairComponentWithQString<Dlcl>(*itr);//获得pair的化合物信息
        pair<float,float> sample_mz = GetPairSampleMz<Dlcl>(*itr);
        float sample_rt = GetPairRt<Dlcl>(*itr);
        stream << component_message << "," << sample_mz.first << "," << sample_mz.second << "," << sample_rt << "," << "DLCL" << endl;
    }
    file.close();//关闭文件
    qDebug() << "输出了:" << fileName;
}











