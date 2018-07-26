#include <math.h>

#include "Behavior.h"

NEAT::Behavior::Behavior() : PhenotypeBehavior()
{
  m_Data = {std::vector<double>()};
}

double NEAT::Behavior::Distance_To(NEAT::PhenotypeBehavior* other)
{  
  double sum = 0;
  
  

  //if (m_Data[0].size() != other->m_Data[0].size())
    std::cout << "size 1: " << m_Data[0].size() << ", size 2: " << other->m_Data[0].size() << std::endl;
  
  for (int i = 0; i < m_Data[0].size(); i++) {
    sum += (m_Data.at(0).at(i) - other->m_Data.at(0).at(i)) * (m_Data.at(0).at(i) - other->m_Data.at(0).at(i));
  }
  
  return sqrt(sum);
}

