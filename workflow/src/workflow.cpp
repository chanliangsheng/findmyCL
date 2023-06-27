#include "workflow.h"

using namespace std;
Workflow::Workflow(QObject *parent) : QObject(parent)
{

}

std::vector<std::pair<Cl, Cl> > Workflow::GetCopyClPairVector()
{
    return this->m_cl_vector;
}

std::vector<std::pair<Mlcl, Mlcl> > Workflow::GetCopyMlclPairVector()
{
    return this->m_mlcl_vector;
}

std::vector<std::pair<Dlcl, Dlcl> > Workflow::GetCopyDlclPairVector()
{
    return this->m_dlcl_vector;
}

void Workflow::DeleteRedundantClPair()
{
    //把仅有M-H或M-2H中重复的部分去除
    //方法：1.先将两者都有的pair的信息加入到cl_message中，2.然后仅有M-H或M-2H的一边在set中检查一边将信息加入到set中，防止仅有M-H或M-2H的自己也会重复
    //1.
    set<pair<Ms1* , DatabaseRecord*>> cl_message_set;
    for(auto itr = this->m_cl_vector.begin() ; itr != this->m_cl_vector.end();itr++){
        //如果M-H和M-2H都存在，把信息假如到cl_message中
        if(!itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject()){
            cl_message_set.insert({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()});
            cl_message_set.insert({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()});
        }
    }
    //2.
    for(auto itr = this->m_cl_vector.begin() ; itr != this->m_cl_vector.end();){
        //如果仅有M-H
        if(!itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            //如果找到重复值，删除该元素
            if(cl_message_set.find({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()}) != cl_message_set.end()){
                this->m_cl_vector.erase(itr);
            }
            //如果没有找到，将信息加入cl_message中
            else{
                cl_message_set.insert({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()});
                itr++;
            }
        }
        //如果仅有M-2H
        else if(itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject())
        {
            //如果找到重复值，删除该元素
            if(cl_message_set.find({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()}) != cl_message_set.end()){
                this->m_cl_vector.erase(itr);
            }
            //如果没有找到，将信息加入cl_message中
            else{
                cl_message_set.insert({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()});
                itr++;
            }
        }
        else{
            itr++;
        }
    }
}

void Workflow::DeleteRedundantMlclPair()
{
    //把仅有M-H或M-2H中重复的部分去除
    //方法：1.先将两者都有的pair的信息加入到cl_message中，2.然后仅有M-H或M-2H的一边在set中检查一边将信息加入到set中，防止仅有M-H或M-2H的自己也会重复
    //1.
    set<pair<Ms1* , DatabaseRecord*>> mlcl_message_set;
    for(auto itr = this->m_mlcl_vector.begin() ; itr != this->m_mlcl_vector.end();itr++){
        //如果M-H和M-2H都存在，把信息假如到cl_message中
        if(!itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject()){
            mlcl_message_set.insert({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()});
            mlcl_message_set.insert({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()});
        }
    }
    //2.
    for(auto itr = this->m_mlcl_vector.begin() ; itr != this->m_mlcl_vector.end();){
        //如果仅有M-H
        if(!itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            //如果找到重复值，删除该元素
            if(mlcl_message_set.find({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()}) != mlcl_message_set.end()){
                this->m_mlcl_vector.erase(itr);
            }
            //如果没有找到，将信息加入cl_message中
            else{
                mlcl_message_set.insert({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()});
                itr++;
            }
        }
        //如果仅有M-2H
        else if(itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject())
        {
            //如果找到重复值，删除该元素
            if(mlcl_message_set.find({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()}) != mlcl_message_set.end()){
                this->m_mlcl_vector.erase(itr);
            }
            //如果没有找到，将信息加入cl_message中
            else{
                mlcl_message_set.insert({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()});
                itr++;
            }
        }
        else{
            itr++;
        }
    }
}

void Workflow::DeleteRedundantDlclPair()
{
    //把仅有M-H或M-2H中重复的部分去除
    //方法：1.先将两者都有的pair的信息加入到cl_message中，2.然后仅有M-H或M-2H的一边在set中检查一边将信息加入到set中，防止仅有M-H或M-2H的自己也会重复
    //1.
    set<pair<Ms1* , DatabaseRecord*>> dlcl_message_set;
    for(auto itr = this->m_dlcl_vector.begin() ; itr != this->m_dlcl_vector.end();itr++){
        //如果M-H和M-2H都存在，把信息假如到dlcl_message中
        if(!itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject()){
            dlcl_message_set.insert({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()});
            dlcl_message_set.insert({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()});
        }
    }
    //2.
    for(auto itr = this->m_dlcl_vector.begin() ; itr != this->m_dlcl_vector.end();){
        //如果仅有M-H
        if(!itr->first.CheckEmptyObject() && itr->second.CheckEmptyObject()){
            //如果找到重复值，删除该元素
            if(dlcl_message_set.find({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()}) != dlcl_message_set.end()){
                this->m_dlcl_vector.erase(itr);
            }
            //如果没有找到，将信息加入dlcl_message中
            else{
                dlcl_message_set.insert({itr->first.GetMs1Ptr() , itr->first.GetMs1DatabaseRecordPtr()});
                itr++;
            }
        }
        //如果仅有M-2H
        else if(itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject())
        {
            //如果找到重复值，删除该元素
            if(dlcl_message_set.find({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()}) != dlcl_message_set.end()){
                this->m_dlcl_vector.erase(itr);
            }
            //如果没有找到，将信息加入dlcl_message中
            else{
                dlcl_message_set.insert({itr->second.GetMs1Ptr() , itr->second.GetMs1DatabaseRecordPtr()});
                itr++;
            }
        }
        else{
            itr++;
        }
    }
}

void Workflow::DeleteRedundantPair()
{
    this->DeleteRedundantClPair();
    this->DeleteRedundantMlclPair();
    this->DeleteRedundantDlclPair();
}

void Workflow::ShowResult(QTableWidget *qtablewidget)
{
    //先清空信息
    qtablewidget->setRowCount(0);

    //在这里设置行，可以显著提升ui性能；而不是每个for中增加一行
    qtablewidget->setRowCount(this->m_cl_vector.size() + this->m_mlcl_vector.size() + this->m_dlcl_vector.size());

    unsigned int row = 0;

    //显示CL
    while(row < this->m_cl_vector.size()){
        if(this->m_cl_vector[row].first.CheckEmptyObject()){
            qtablewidget->setItem(row , 0 , 0);
            qtablewidget->setItem(row , 1 , 0);
            qtablewidget->setItem(row , 2 , 0);
        }
        else{
            QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(this->m_cl_vector[row].first.GetSampleMz()));
            QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(this->m_cl_vector[row].first.GetSampleRt()));
            QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(this->m_cl_vector[row].first.GetChainLength()) + ":" + QString::number(this->m_cl_vector[row].first.GetUnsaturation()) + ":" + QString::number(this->m_cl_vector[row].first.GetOxygen()));
            qtablewidget->setItem(row , 0 , M_H_mz);
            qtablewidget->setItem(row , 1 , M_H_rt);
            qtablewidget->setItem(row , 2 , M_H_struct);

        }

        if(this->m_cl_vector[row].second.CheckEmptyObject()){
            qtablewidget->setItem(row , 3 , 0);
            qtablewidget->setItem(row , 4 , 0);
            qtablewidget->setItem(row , 5 , 0);
        }
        else{
            QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(this->m_cl_vector[row].second.GetSampleMz()));
            QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(this->m_cl_vector[row].second.GetSampleRt()));
            QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(this->m_cl_vector[row].second.GetChainLength()) + ":" + QString::number(this->m_cl_vector[row].second.GetUnsaturation()) + ":" + QString::number(this->m_cl_vector[row].second.GetOxygen()));
            qtablewidget->setItem(row , 3 , M_2H_mz);
            qtablewidget->setItem(row , 4 , M_2H_rt);
            qtablewidget->setItem(row , 5 , M_2H_struct);
        }
        row++;
    }

    //显示MLCL
    for(unsigned int mlcl_row = 0 ; mlcl_row < this->m_mlcl_vector.size() ; mlcl_row ++){
        if(this->m_mlcl_vector[mlcl_row].first.CheckEmptyObject()){
            qtablewidget->setItem(row , 0 , 0);
            qtablewidget->setItem(row , 1 , 0);
            qtablewidget->setItem(row , 2 , 0);
        }
        else{
            QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(this->m_mlcl_vector[mlcl_row].first.GetSampleMz()));
            QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(this->m_mlcl_vector[mlcl_row].first.GetSampleRt()));
            QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(this->m_mlcl_vector[mlcl_row].first.GetChainLength()) + ":" + QString::number(this->m_mlcl_vector[mlcl_row].first.GetUnsaturation()) + ":" + QString::number(this->m_mlcl_vector[mlcl_row].first.GetOxygen()));
            qtablewidget->setItem(row , 0 , M_H_mz);
            qtablewidget->setItem(row , 1 , M_H_rt);
            qtablewidget->setItem(row , 2 , M_H_struct);
        }

        if(this->m_mlcl_vector[mlcl_row].second.CheckEmptyObject()){
            qtablewidget->setItem(row , 3 , 0);
            qtablewidget->setItem(row , 4 , 0);
            qtablewidget->setItem(row , 5 , 0);
        }
        else{
            QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(this->m_mlcl_vector[mlcl_row].second.GetSampleMz()));
            QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(this->m_mlcl_vector[mlcl_row].second.GetSampleRt()));
            QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(this->m_mlcl_vector[mlcl_row].second.GetChainLength()) + ":" + QString::number(this->m_mlcl_vector[mlcl_row].second.GetUnsaturation()) + ":" + QString::number(this->m_mlcl_vector[mlcl_row].second.GetOxygen()));
            qtablewidget->setItem(row , 3 , M_2H_mz);
            qtablewidget->setItem(row , 4 , M_2H_rt);
            qtablewidget->setItem(row , 5 , M_2H_struct);
        }
        row++;
    }

    //显示DLCL
    for(unsigned int dlcl_row = 0 ; dlcl_row < this->m_dlcl_vector.size() ; dlcl_row++){
        if(this->m_dlcl_vector[dlcl_row].first.CheckEmptyObject()){
            qtablewidget->setItem(row , 0 , 0);
            qtablewidget->setItem(row , 1 , 0);
            qtablewidget->setItem(row , 2 , 0);
        }
        else{
            QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(this->m_dlcl_vector[dlcl_row].first.GetSampleMz()));
            QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(this->m_dlcl_vector[dlcl_row].first.GetSampleRt()));
            QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(this->m_dlcl_vector[dlcl_row].first.GetChainLength()) + ":" + QString::number(this->m_dlcl_vector[dlcl_row].first.GetUnsaturation()) + ":" + QString::number(this->m_dlcl_vector[dlcl_row].first.GetOxygen()));
            qtablewidget->setItem(row , 0 , M_H_mz);
            qtablewidget->setItem(row , 1 , M_H_rt);
            qtablewidget->setItem(row , 2 , M_H_struct);
        }

        if(this->m_dlcl_vector[dlcl_row].second.CheckEmptyObject()){
            qtablewidget->setItem(row , 3 , 0);
            qtablewidget->setItem(row , 4 , 0);
            qtablewidget->setItem(row , 5 , 0);
        }
        else{
            QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(this->m_dlcl_vector[dlcl_row].second.GetSampleMz()));
            QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(this->m_dlcl_vector[dlcl_row].second.GetSampleRt()));
            QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(this->m_dlcl_vector[dlcl_row].second.GetChainLength()) + ":" + QString::number(this->m_dlcl_vector[dlcl_row].second.GetUnsaturation()) + ":" + QString::number(this->m_dlcl_vector[dlcl_row].second.GetOxygen()));
            qtablewidget->setItem(row , 3 , M_2H_mz);
            qtablewidget->setItem(row , 4 , M_2H_rt);
            qtablewidget->setItem(row , 5 , M_2H_struct);
        }
        row++;
    }
}

void Workflow::ShowResult(QTableWidget *qtablewidget, float min_rt, float max_rt)
{
    //先清空信息
    qtablewidget->setRowCount(0);
    //在这里设置行，可以显著提升ui性能；而不是每个for中增加一行
    qtablewidget->setRowCount(this->m_cl_vector.size() + this->m_mlcl_vector.size() + this->m_dlcl_vector.size());


    unsigned int row = 0;

    //显示CL
    for(auto itr = this->m_cl_vector.begin() ; itr != this->m_cl_vector.end();itr++){
        bool success = 0;
        if(!itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject()){
            //只要M-H或M-2H一个满足时间范围即可
            if((itr->first.GetSampleRt() > min_rt && itr->first.GetSampleRt() < max_rt) || (itr->second.GetSampleRt() > min_rt && itr->second.GetSampleRt() < max_rt)){
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first.GetChainLength()) + ":" + QString::number(itr->first.GetUnsaturation()) + ":" + QString::number(itr->first.GetOxygen()));
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->second.GetChainLength()) + ":" + QString::number(itr->second.GetUnsaturation()) + ":" + QString::number(itr->second.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }
        }
        //如果只存在M-H，不存在M-2H
        else if(!itr->first.CheckEmptyObject()){
            if(itr->first.GetSampleRt() > min_rt && itr->first.GetSampleRt() < max_rt){
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first.GetChainLength()) + ":" + QString::number(itr->first.GetUnsaturation()) + ":" + QString::number(itr->first.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                success = 1;
            }
        }
        //如果只存在M-2H，不存在M-H
        else if(!itr->second.CheckEmptyObject()){
            if(itr->second.GetSampleRt() > min_rt && itr->second.GetSampleRt() < max_rt){
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->second.GetChainLength()) + ":" + QString::number(itr->second.GetUnsaturation()) + ":" + QString::number(itr->second.GetOxygen()));
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }

        }
        if(success){
            row++;
        }
    }
    //显示MLCL
    for(auto itr = this->m_mlcl_vector.begin() ; itr != this->m_mlcl_vector.end();itr++){
        bool success = 0;
        if(!itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject()){
            if((itr->first.GetSampleRt() > min_rt && itr->first.GetSampleRt() < max_rt) || (itr->second.GetSampleRt() > min_rt && itr->second.GetSampleRt() < max_rt)){
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first.GetChainLength()) + ":" + QString::number(itr->first.GetUnsaturation()) + ":" + QString::number(itr->first.GetOxygen()));
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->second.GetChainLength()) + ":" + QString::number(itr->second.GetUnsaturation()) + ":" + QString::number(itr->second.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }

        }
        else if(!itr->first.CheckEmptyObject()){
            if(itr->first.GetSampleRt() > min_rt && itr->first.GetSampleRt() < max_rt){
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first.GetChainLength()) + ":" + QString::number(itr->first.GetUnsaturation()) + ":" + QString::number(itr->first.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                success = 1;
            }
        }
        else if(!itr->second.CheckEmptyObject()){
            if(itr->second.GetSampleRt() > min_rt && itr->second.GetSampleRt() < max_rt){
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->second.GetChainLength()) + ":" + QString::number(itr->second.GetUnsaturation()) + ":" + QString::number(itr->second.GetOxygen()));
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }

        }
        if(success){
            row++;
        }
    }
    //显示DLCL
    for(auto itr = this->m_dlcl_vector.begin() ; itr != this->m_dlcl_vector.end();itr++){
        bool success = 0;
        if(!itr->first.CheckEmptyObject() && !itr->second.CheckEmptyObject()){
            if((itr->first.GetSampleRt() > min_rt && itr->first.GetSampleRt() < max_rt) || (itr->second.GetSampleRt() > min_rt && itr->second.GetSampleRt() < max_rt)){
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first.GetChainLength()) + ":" + QString::number(itr->first.GetUnsaturation()) + ":" + QString::number(itr->first.GetOxygen()));
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->second.GetChainLength()) + ":" + QString::number(itr->second.GetUnsaturation()) + ":" + QString::number(itr->second.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }

        }
        else if(!itr->first.CheckEmptyObject()){
            if(itr->first.GetSampleRt() > min_rt && itr->first.GetSampleRt() < max_rt){
                QTableWidgetItem* M_H_mz = new QTableWidgetItem(QString::number(itr->first.GetSampleMz()));
                QTableWidgetItem* M_H_rt = new QTableWidgetItem(QString::number(itr->first.GetSampleRt()));
                QTableWidgetItem* M_H_struct =  new QTableWidgetItem(QString::number(itr->first.GetChainLength()) + ":" + QString::number(itr->first.GetUnsaturation()) + ":" + QString::number(itr->first.GetOxygen()));
                qtablewidget->setItem(row , 0 , M_H_mz);
                qtablewidget->setItem(row , 1 , M_H_rt);
                qtablewidget->setItem(row , 2 , M_H_struct);
                success = 1;
            }
        }
        else if(!itr->second.CheckEmptyObject()){
            if(itr->second.GetSampleRt() > min_rt && itr->second.GetSampleRt() < max_rt){
                QTableWidgetItem* M_2H_mz = new QTableWidgetItem(QString::number(itr->second.GetSampleMz()));
                QTableWidgetItem* M_2H_rt = new QTableWidgetItem(QString::number(itr->second.GetSampleRt()));
                QTableWidgetItem* M_2H_struct =  new QTableWidgetItem(QString::number(itr->second.GetChainLength()) + ":" + QString::number(itr->second.GetUnsaturation()) + ":" + QString::number(itr->second.GetOxygen()));
                qtablewidget->setItem(row , 3 , M_2H_mz);
                qtablewidget->setItem(row , 4 , M_2H_rt);
                qtablewidget->setItem(row , 5 , M_2H_struct);
                success = 1;
            }

        }
        if(success){
            row++;
        }
    }
    qtablewidget->setRowCount(row);//设置为真正有多少行
}


