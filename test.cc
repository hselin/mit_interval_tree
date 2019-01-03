#include "interval_tree.h"
#include <stdio.h>
#include <assert.h>


class SimpleInterval : public Interval {
public:
  SimpleInterval(const int64_t low,const int64_t high)
    :_low(low), _high(high)
    { }
  
  int64_t GetLowPoint() const { return _low;}
  int64_t GetHighPoint() const { return _high;}
  IntervalTreeNode * GetNode() { return _node;}
  void SetNode(IntervalTreeNode * node) {_node = node;}
protected:
  int64_t _low;
  int64_t _high;
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
  return tree->Enumerate((int64_t)low, (int64_t)high, qr);
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

  /*
  SimpleInterval *i = (SimpleInterval *)addInterval(tree, 0, 9);
  printf("i: %p %lu %lu\n", (void *)i, i->GetLowPoint(), i->GetHighPoint());
  Interval *deletedInterval = tree->DeleteNode(i->GetNode());
  printf("deletedInterval: %p %lu %lu\n", (void *)deletedInterval, deletedInterval->GetLowPoint(), deletedInterval->GetHighPoint());

  //delete twice will be bad
  //deletedInterval = tree->DeleteNode(i->GetNode());
  //printf("deletedInterval: %p %lu %lu\n", (void *)deletedInterval, deletedInterval->GetLowPoint(), deletedInterval->GetHighPoint());


  return 0;
  */


  /*
  addInterval(tree, 99, 100);
  addInterval(tree, 100, 101);
  addInterval(tree, 100, 100);
  addInterval(tree, 100, 100);

  deleteInterval(tree, 100, 101);
  */

  /*
  printf("INT64_MIN: %ld\n", INT64_MIN);
  printf("INT64_MAX: %ld\n", INT64_MAX);

  printf("LLONG_MIN: %lld\n", LLONG_MIN);
  printf("LLONG_MAX: %lld\n", LLONG_MAX);
  */

  addInterval(tree, 0, 0);
  addInterval(tree, 1, 1);
  addInterval(tree, 2, 2);
  addInterval(tree, 3, 3);
  addInterval(tree, 4, 4);
  addInterval(tree, 5, 5);
  addInterval(tree, 6, 6);
  addInterval(tree, 7, 7);
  addInterval(tree, 8, 8);
  addInterval(tree, 9, 9);
  
  //std::queue<Interval *> qr;
  //int count = queryInterval(tree, 0, 7047, (std::queue<void *> *)&qr);
  //assert(count == qr.size());
  //printIntervals(&qr);

  tree->Print();

  return 0;
}