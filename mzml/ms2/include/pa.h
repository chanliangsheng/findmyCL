#ifndef PA_H
#define PA_H
#include <databaserecord.h>

class Pa
{
public:
    Pa();
    Pa(float mz ,float intensity ,float score , DatabaseRecord* database_record);
public:
    void SetMz(float mz);
    void SetIntensity(float intensity);
    void SetScore(float score);
    void SetDatabaseRecordPtr(DatabaseRecord* database_record);
    float GetMz();
    float GetIntensity();
    float GetScore();
    unsigned int GetChainLength();
    unsigned int GetUnsaturation();
    unsigned int GetOxygen();
    std::string GetAdditiveForm();
    std::string GetFormula();
private:
    float m_mz;
    float m_intensity;
    float m_score;
    DatabaseRecord* m_database_record;
};


class Fa:public Pa
{
public:
    Fa();
    Fa(float mz ,float intensity ,float score , DatabaseRecord* database_record);
};

#endif // PA_H
