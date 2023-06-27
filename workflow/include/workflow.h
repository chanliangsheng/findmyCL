#ifndef WORKFLOW_H
#define WORKFLOW_H
#include <utility>
#include <database.h>
#include <map>
#include <QTime>
#include <CL.h>
#include <MLCL.h>
#include <DLCL.h>
#include <MLCL.h>
#include <DLCL.h>
#include <QObject>
#include <QTableWidget>
#include <list>
#include <set>

class Workflow:public QObject
{
    Q_OBJECT
public:
    Workflow(QObject *parent = nullptr);
protected:
    std::vector<std::pair<Cl , Cl>> m_cl_vector;
    std::vector<std::pair<Mlcl , Mlcl>> m_mlcl_vector;
    std::vector<std::pair<Dlcl , Dlcl>> m_dlcl_vector;
public:
    std::vector<std::pair<Cl , Cl>> GetCopyClPairVector();
    std::vector<std::pair<Mlcl , Mlcl>> GetCopyMlclPairVector();
    std::vector<std::pair<Dlcl , Dlcl>> GetCopyDlclPairVector();
private:
    void DeleteRedundantClPair();
    void DeleteRedundantMlclPair();
    void DeleteRedundantDlclPair();
public:
    void DeleteRedundantPair();
public:
    //获得化合物形式
    template <typename T>
    QString GetPairComponentWithQString(std::pair<T , T> pair){
        QString first_result =  pair.first.GetMs1DatabaseRecordComponentWithQString();
        if(first_result != "0:0:0"){
            return first_result;
        }
        else{
            return pair.second.GetMs1DatabaseRecordComponentWithQString();
        }
    };//获得化合物的链长:不饱和度:氧个数
    //获得mz
    template <typename T>
    std::pair<float , float> GetPairSampleMz(std::pair<T , T> pair){
        return std::make_pair<float,float>(pair.first.GetSampleMz(),pair.second.GetSampleMz());
    };
    //获得rt
    template <typename T>
    float GetPairRt(std::pair<T , T> pair){
        float first_rt = pair.first.GetSampleRt();
        float second_rt = pair.second.GetSampleRt();
        //如果某个是空对象，则返回另一个；如果两者都不是空对象，则返回rt的平均值
        if(first_rt == 0){
            return second_rt;
        }
        else if(second_rt == 0){
            return first_rt;
        }
        else if(first_rt != 0 && second_rt != 0){
            return (first_rt + second_rt)/2 ;
        }
    };
    //获得intensity
    template <typename T>
    float GetPairIntensity(std::pair<T , T> pair){
        float first_intensity = pair.first.GetSampleIntensity();
        float second_intensity = pair.second.GetSampleIntensity();
        //如果某个是空对象，则返回另一个；如果两者都不是空对象，则返回intensity的最大值
        if(first_intensity == 0){
            return second_intensity;
        }
        else if(second_intensity == 0){
            return first_intensity;
        }
        else if(first_intensity != 0 && second_intensity != 0){
            return std::max(first_intensity,second_intensity) ;
        }
    };

public slots:
    virtual void ShowResult(QTableWidget* qtablewidget);
    virtual void ShowResult(QTableWidget* qtablewidget , float min_rt, float max_rt);
};

#endif // WORKFLOW_H
