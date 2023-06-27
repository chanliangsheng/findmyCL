#ifndef PANODE_H
#define PANODE_H
#include <database.h>
#include <mzml.h>
#include <memory>
#include <cmath>
#include <headgroup.h>
#include <set>
#include <pa.h>

//心磷脂精细结构的PA节点信息
class PaNode
{
//心磷脂精细结构的PA节点信息
public:
    Pa* m_pa_ptr;
    Fa* m_left_fa_ptr;
    Fa* m_right_fa_ptr;
};

#endif // PANODE_H
