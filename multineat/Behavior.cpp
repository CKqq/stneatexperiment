#include <math.h>

#include "Behavior.h"

NEAT::Behavior::Behavior() : PhenotypeBehavior()
{
   m_Data = {std::vector<double>(6)};
}

double NEAT::Behavior::Distance_To(NEAT::PhenotypeBehavior* other)
{  
  double sum = 0.0;
  
  Behavior* b = (Behavior*) other;
  
  
  if (m_Data[0].size() != other->m_Data[0].size()) {
    std::cout << "Different sizes, " << m_Data[0].size() << " and " << b->m_Data[0].size() << std::endl;
  }
  
  for (int i = 0; i < m_Data[0].size() && i < b->m_Data[0].size(); i++) {
    double x1 = m_Data[0][i];
    double x2 = b->m_Data[0][i];
    double difference = x1 - x2;
    difference *= difference;
    sum += difference;
  }
  
  return sqrt(sum);
}

