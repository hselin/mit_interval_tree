#include "interval_tree.h"
#include <stdio.h>
#include <assert.h>


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

int queryInterval(IntervalTree *tree, uint64_t low, uint64_t high, std::queue<void *> *qr)
{
  return tree->Enumerate(low, high, qr);
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
  std::queue<Interval *> qr;
  int count = queryInterval(tree, low, high, (std::queue<void *> *)&qr);

  //printf("count: %d\n", count);
  //printf("qr.size(): %lu\n", qr.size());

  assert(count == qr.size());

  while(!qr.empty())
  {
    SimpleInterval *i = (SimpleInterval *)qr.front();

    if((i->GetLowPoint() == low) && (i->GetHighPoint() == high))
    {
      //found our interval
      tree->DeleteNode(i->GetNode());
    }

    qr.pop();
  }
}


int main(int argc, char ** argv)
{
  IntervalTree *tree = new IntervalTree();

  addInterval(tree, 99, 100);
  addInterval(tree, 100, 101);
  addInterval(tree, 100, 100);
  addInterval(tree, 100, 100);

  deleteInterval(tree, 100, 101);

  std::queue<Interval *> qr;
  int count = queryInterval(tree, 100, 100, (std::queue<void *> *)&qr);
  assert(count == qr.size());
  printIntervals(&qr);

  tree->Print();

  return 0;
}