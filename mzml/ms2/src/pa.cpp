#include <pa.h>

Pa::Pa()
{
    this->m_mz = 0;
    this->m_intensity = 0;
    this->m_score = 0;
    this->m_database_record = nullptr;
}

Pa::Pa(float mz, float intensity, float score, DatabaseRecord *database_record)
{
    this->m_mz = mz;
    this->m_intensity = intensity;
    this->m_score = score;
    this->m_database_record = database_record;
}



void Pa::SetMz(float mz)
{
    this->m_mz = mz;
}

void Pa::SetIntensity(float intensity)
{
    this->m_intensity = intensity;
}

void Pa::SetScore(float score)
{
    this->m_score = score;
}

void Pa::SetDatabaseRecordPtr(DatabaseRecord *database_record)
{
    this->m_database_record = database_record;
}

float Pa::GetMz()
{
    return this->m_mz;
}

float Pa::GetIntensity()
{
    return this->m_intensity;
}

float Pa::GetScore()
{
    return this->m_score;
}

unsigned int Pa::GetChainLength()
{
    return this->m_database_record->GetChainLength();
}

unsigned int Pa::GetUnsaturation()
{
    return this->m_database_record->GetUnsaturation();
}

unsigned int Pa::GetOxygen()
{
    return this->m_database_record->GetOxygen();
}

std::string Pa::GetAdditiveForm()
{
    return this->m_database_record->GetAdditiveForm();
}

std::string Pa::GetFormula()
{
    return  this->m_database_record->GetFormula();
}


Fa::Fa()
{
    this->SetMz(0);
    this->SetIntensity(0);
    this->SetScore(0);
    this->SetDatabaseRecordPtr(nullptr);
}

Fa::Fa(float mz, float intensity, float score, DatabaseRecord *database_record)
{
    this->SetMz(mz);
    this->SetIntensity(intensity);
    this->SetScore(score);
    this->SetDatabaseRecordPtr(database_record);
}
