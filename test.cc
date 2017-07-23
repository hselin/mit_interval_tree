#include "interval_tree.h"
#include <stdio.h>


class SimpleInterval : public Interval {
public:
  SimpleInterval(const int low,const int high)
    :_low(low), _high(high)
    { }
  
  int GetLowPoint() const { return _low;}
  int GetHighPoint() const { return _high;}
  IntervalTreeNode * GetNode() { return _node;}
  void SetNode(IntervalTreeNode * node) {_node = node;}
protected:
  int _low;
  int _high;
  IntervalTreeNode * _node;
};


Interval *addInterval(IntervalTree *tree, int low, int high)
{
  SimpleInterval * newSimpleInterval = new SimpleInterval(low,high);
  newSimpleInterval->SetNode(tree->Insert(newSimpleInterval));

  return newSimpleInterval;
}

std::queue<void *> *queryInterval(IntervalTree *tree, int low, int high)
{
  return tree->Enumerate(low, high);
}

void printIntervals(std::queue<Interval *> *intervals)
{
  while(!intervals->empty())
  {
    Interval *i = intervals->front();

    printf("[%d - %d]\n", i->GetLowPoint(), i->GetHighPoint());
    intervals->pop();
  }
}

void deleteInterval(IntervalTree *tree, int low, int high)
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

  deleteInterval(tree, 100, 100);

  std::queue<Interval *> *qr = (std::queue<Interval *> *) queryInterval(tree, 100, 100);
  printIntervals(qr);

  tree->Print();

  return 0;
}