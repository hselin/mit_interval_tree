#include "interval_tree.h"
#include <stdio.h>


class SimpleInterval : public Interval {
public:
  SimpleInterval(const uint64_t low,const uint64_t high)
    :_low(low), _high(high)
    { }
  
  uint64_t GetLowPoint() const { return _low;}
  uint64_t GetHighPoint() const { return _high;}
  IntervalTreeNode * GetNode() { return _node;}
  void SetNode(IntervalTreeNode * node) {_node = node;}
protected:
  uint64_t _low;
  uint64_t _high;
  IntervalTreeNode * _node;
};


Interval *addInterval(IntervalTree *tree, uint64_t low, uint64_t high)
{
  SimpleInterval * newSimpleInterval = new SimpleInterval(low,high);
  newSimpleInterval->SetNode(tree->Insert(newSimpleInterval));

  return newSimpleInterval;
}

std::queue<void *> *queryInterval(IntervalTree *tree, uint64_t low, uint64_t high)
{
  return tree->Enumerate(low, high);
}

void printIntervals(std::queue<Interval *> *intervals)
{
  while(!intervals->empty())
  {
    Interval *i = intervals->front();

    printf("[%lu - %lu]\n", i->GetLowPoint(), i->GetHighPoint());
    intervals->pop();
  }
}

void deleteInterval(IntervalTree *tree, uint64_t low, uint64_t high)
{
  std::queue<Interval *> *qr = (std::queue<Interval *> *)queryInterval(tree, low, high);

  while(!qr->empty())
  {
    SimpleInterval *i = (SimpleInterval *)qr->front();

    if((i->GetLowPoint() == low) && (i->GetHighPoint() == high))
    {
      //found our interval
      tree->DeleteNode(i->GetNode());
    }

    qr->pop();
  }
}


int main(int argc, char ** argv)
{
  IntervalTree *tree = new IntervalTree();

  addInterval(tree, 99, 100);
  addInterval(tree, 100, 101);
  addInterval(tree, 100, 100);
  addInterval(tree, 100, 100);

  deleteInterval(tree, 0, 0);

  std::queue<Interval *> *qr = (std::queue<Interval *> *) queryInterval(tree, 0, 0);
  printIntervals(qr);

  tree->Print();

  return 0;
}