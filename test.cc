#include "interval_tree.h"
#include <stdio.h>
#include <string>
#include <cstdint>

template <class T>
class SimpleInterval : public Interval <T> {
public:
  SimpleInterval(const T low,const T high)
    :_low(low), _high(high)
    { }
  
  T GetLowPoint() const { return _low;}
  T GetHighPoint() const { return _high;}
  IntervalTreeNode<T> * GetNode() { return _node;}
  void SetNode(IntervalTreeNode<T> *node) {_node = node;}
protected:
  T _low;
  T _high;
  IntervalTreeNode<T> * _node;
};

/*
Interval<uint64_t> *addInterval(IntervalTree<uint64_t> *tree, uint64_t low, uint64_t high)
{
  SimpleInterval<uint64_t> *newSimpleInterval = new SimpleInterval<uint64_t>(low,high);
  newSimpleInterval->SetNode(tree->Insert(newSimpleInterval));

  return newSimpleInterval;
}

std::queue<void *> *queryInterval(IntervalTree<uint64_t> *tree, uint64_t low, uint64_t high)
{
  return tree->Enumerate(low, high);
}

void printIntervals(std::queue<Interval<uint64_t> *> *intervals)
{
  while(!intervals->empty())
  {
    Interval<uint64_t> *i = intervals->front();

    printf("[%lu - %lu]\n", i->GetLowPoint(), i->GetHighPoint());
    intervals->pop();
  }
}

void deleteInterval(IntervalTree<uint64_t> *tree, uint64_t low, uint64_t high)
{
  std::queue<Interval<uint64_t> *> *qr = (std::queue<Interval<uint64_t> *> *)queryInterval(tree, low, high);

  while(!qr->empty())
  {
    SimpleInterval<uint64_t> *i = (SimpleInterval<uint64_t> *)qr->front();

    if((i->GetLowPoint() == low) && (i->GetHighPoint() == high))
    {
      //found our interval
      tree->DeleteNode(i->GetNode());
    }

    qr->pop();
  }
}
*/


int main(int argc, char ** argv)
{
  IntervalTree<uint64_t> *tree = new IntervalTree<uint64_t>();

/*
  addInterval(tree, 99, 100);
  addInterval(tree, 100, 101);
  addInterval(tree, 100, 100);
  addInterval(tree, 100, 100);

  deleteInterval(tree, 100, 100);

  std::queue<Interval<uint64_t> *> *qr = (std::queue<Interval<uint64_t> *> *) queryInterval(tree, 100, 100);
  printIntervals(qr);

  tree->Print();
*/

  return 0;
}
