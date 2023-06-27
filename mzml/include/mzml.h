#ifndef MZML_H
#define MZML_H

#include <tinyxml2.h>
#include <ms1.h>
#include <ms2.h>
#include <QObject>
#include <QString>
#include <tinyxml2.h>
#include <QDebug>
#include <zlib.h>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <memory>
#include <fstream>
#include <list>
#include <QFile>
#include <base64.h>

class Mzml: public QObject
{
    Q_OBJECT
public:
    Mzml(QObject *parent = nullptr);
private:
    //所有的ms1和ms2存放的位置
    std::vector<Ms1> m_ms1_vector;
    std::vector<Ms2> m_ms2_vector;
    //数据的编码与压缩
    const char* m_bit_type_param = "32-bit float";
    const char* m_compression_param = "no compression";
//公共接口
public:
    void ReadMs2FromMzmls(QStringList file_names);//读取多个mzml文件中的二级
    void ReadMs1FromCsv(QString file_name);//从excel表格中获取一级
    void ReadMs2FromMzml(QString file_name);//读取mzml文件中的二级
    std::vector<Ms2>& GetLocalMs2Vector();//获取本地的MS2向量，不复制
    std::vector<Ms1>& GetLocalMs1Vector();//获取本地的MS1向量，不复制
    void ConvertMs1RtUnit(QString mode = "mintosec");//转换MS1保留时间的单位，默认是从分钟转化为秒
    void ConvertMs2RtUnit(QString mode = "mintosec");//转换MS2保留时间的单位，默认是从分钟转化为秒
    void DeleteMs2LowIntensityFragment(float radio);
private:
    void ParserMs1(tinyxml2::XMLElement *spectrum_node);//解析一级节点
    void ParserMs2(tinyxml2::XMLElement *spectrum_node);//解析一级节点
    char GetMsLevel(tinyxml2::XMLElement *spectrum_node);//获取该节点是属于一级节点还是二级节点
    void GetEncodeCompressionParam(tinyxml2::XMLElement *spectrum_node);//获取数据编码的压缩的方法
    std::shared_ptr<std::vector<float>> BytesToFloat(std::string &byte_array);//转化byte数组成float数组，智能指针
    std::shared_ptr<std::vector<double>> BytesToDouble(std::string &byte_array);//转化byte数组成double数组，智能指针
    std::string ZlibDecompress(const std::string& str);//解压zlib的字符串
public:
    void SortMs2VectorByPrecursorIonMz();
private:
    bool m_ms2_sort_by_precuisor_ion_mz = 0;

signals:
    void FinishReadAllMs2();
    void FinishReadAllMs1();
};

#endif // MZML_H
