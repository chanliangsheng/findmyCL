#ifndef DATABASE_H
#define DATABASE_H

#include <vector>
#include <string>
#include <QString>
#include <QSqlDatabase>
#include <QSqlError>
#include <QSqlQuery>
#include <QDebug>
#include <databaserecord.h>
#include <memory>

class Database:public QObject
{
    Q_OBJECT
public:
    Database(QObject *parent = nullptr);
public:
    void LoadAllTable();//加载数据库中的所有数据表
    std::pair<std::vector<DatabaseRecord>* , std::vector<DatabaseRecord>*> GetLocalClPair();//返回CL的M-H和M-2H的向量对
    std::pair<std::vector<DatabaseRecord>* , std::vector<DatabaseRecord>*> GetLocalMlclPair();//返回MLCL的M-H和M-2H的向量对
    std::pair<std::vector<DatabaseRecord>* , std::vector<DatabaseRecord>*> GetLocalDlclPair();//返回DLCL的M-H和M-2H的向量对
    std::vector<DatabaseRecord>* GetLocalPA();
    std::vector<DatabaseRecord>* GetLocalFA();
private:
    //数据库
    QSqlDatabase m_database;
    std::vector<DatabaseRecord> m_CL_M_H;
    std::vector<DatabaseRecord> m_CL_M_2H;
    std::vector<DatabaseRecord> m_MLCL_M_H;
    std::vector<DatabaseRecord> m_MLCL_M_2H;
    std::vector<DatabaseRecord> m_DLCL_M_H;
    std::vector<DatabaseRecord> m_DLCL_M_2H;
    std::vector<DatabaseRecord> m_PA;
    std::vector<DatabaseRecord> m_FA;
private:
    void LoadDatebase(QString file_name = "database.db");//加载数据库
    std::shared_ptr<std::vector<DatabaseRecord>> LoadSingelTable(QString table);//加载数据库中的单个数据表
};







#endif // DATABASE_H
