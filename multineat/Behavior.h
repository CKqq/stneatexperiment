#ifndef BEHAVIOR_H
#define BEHAVIOR_H

#include "PhenotypeBehavior.h"

namespace NEAT 
{
  
class Behavior : public PhenotypeBehavior
{
public:
  Behavior();
  
  double Distance_To(PhenotypeBehavior* other) override;
};

};

#endif