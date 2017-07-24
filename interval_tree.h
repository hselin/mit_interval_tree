#ifndef E_INTERVAL_TREE
#define E_INTERVAL_TREE


//#include "misc.h"
//#include "TemplateStack.H"
#include <math.h>
#include <limits.h>
#include <iostream>
#include <stack>
#include <queue>




//  The interval_tree.h and interval_tree.cc files contain code for 
//  interval trees implemented using red-black-trees as described in 
//  the book _Introduction_To_Algorithms_ by Cormen, Leisserson, 
//  and Rivest.  

//  CONVENTIONS:  
//                Function names: Each word in a function name begins with 
//                a capital letter.  An example funcntion name is  
//                CreateRedTree(a,b,c). Furthermore, each function name 
//                should begin with a capital letter to easily distinguish 
//                them from variables. 
//                                                                     
//                Variable names: Each word in a variable name begins with 
//                a capital letter EXCEPT the first letter of the variable 
//                name.  For example, int newLongInt.  Global variables have 
//                names beginning with "g".  An example of a global 
//                variable name is gNewtonsConstant. 


#ifndef MAX_INT
#define MAX_INT INT_MAX // some architechturs define INT_MAX not MAX_INT
#endif

// The Interval class is an Abstract Base Class.  This means that no
// instance of the Interval class can exist.  Only classes which
// inherit from the Interval class can exist.  Furthermore any class
// which inherits from the Interval class must define the member
// functions GetLowPoint and GetHighPoint.
//
// The GetLowPoint should return the lowest point of the interval and
// the GetHighPoint should return the highest point of the interval.  

class IntervalTree;
class Interval;
class IntervalTreeNode;

template <class T>
class Interval {
public:
  Interval();
  virtual ~Interval();
  virtual T GetLowPoint() const = 0;
  virtual T GetHighPoint() const = 0;
  virtual void Print() const;
};

template <class T>
class IntervalTreeNode {
  friend class IntervalTree;
public:
  void Print(IntervalTreeNode *,
	     IntervalTreeNode*) const;
  IntervalTreeNode();
  IntervalTreeNode(Interval<T> *);
  ~IntervalTreeNode();
protected:
  Interval<T> * storedInterval;
  T key;
  T high;
  T maxHigh;
  int red; /* if red=0 then the node is black */
  IntervalTreeNode * left;
  IntervalTreeNode * right;
  IntervalTreeNode * parent;
};

template <class T>
struct it_recursion_node {
public:
  /*  this structure stores the information needed when we take the */
  /*  right branch in searching for intervals but possibly come back */
  /*  and check the left branch as well. */

  IntervalTreeNode<T> * start_node;
  unsigned int parentIndex;
  int tryRightBranch;
} ;

template <class T>
class IntervalTree {
public:
  IntervalTree();
  ~IntervalTree();
  void Print() const;
  Interval<T> * DeleteNode(IntervalTreeNode<T> *);
  IntervalTreeNode<T> * Insert(Interval<T> *);
  IntervalTreeNode<T> * GetPredecessorOf(IntervalTreeNode<T> *) const;
  IntervalTreeNode<T> * GetSuccessorOf(IntervalTreeNode<T> *) const;
  //TemplateStack<void *> * Enumerate(int low, int high) ;
  std::queue<void *> * Enumerate(T low, T high);
  void CheckAssumptions() const;
protected:
  /*  A sentinel is used for root and for nil.  These sentinels are */
  /*  created when ITTreeCreate is caled.  root->left should always */
  /*  point to the node which is the root of the tree.  nil points to a */
  /*  node which should always be black but has aribtrary children and */
  /*  parent and no key or info.  The point of using these sentinels is so */
  /*  that the root and nil nodes do not require special cases in the code */
  IntervalTreeNode<T> * root;
  IntervalTreeNode<T> * nil;
  void LeftRotate(IntervalTreeNode<T> *);
  void RightRotate(IntervalTreeNode<T> *);
  void TreeInsertHelp(IntervalTreeNode<T> *);
  void TreePrintHelper(IntervalTreeNode<T> *) const;
  void FixUpMaxHigh(IntervalTreeNode<T> *);
  void DeleteFixUp(IntervalTreeNode<T> *);
  void CheckMaxHighFields(IntervalTreeNode<T> *) const;
  int CheckMaxHighFieldsHelper(IntervalTreeNode<T> * y, 
			const T currentHigh,
			int match) const;
private:
  unsigned int recursionNodeStackSize;
  it_recursion_node<T> * recursionNodeStack;
  unsigned int currentParent;
  unsigned int recursionNodeStackTop;
};









#include "interval_tree.h"
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// If the symbol CHECK_INTERVAL_TREE_ASSUMPTIONS is defined then the
// code does a lot of extra checking to make sure certain assumptions
// are satisfied.  This only needs to be done if you suspect bugs are
// present or if you make significant changes and want to make sure
// your changes didn't mess anything up.

#define CHECK_INTERVAL_TREE_ASSUMPTIONS 1


#define VERIFY(x) assert(x)

const int MIN_INT=-MAX_INT;

// define a function to find the maximum of two objects.
#define ITMax(a, b) ( (a > b) ? a : b )



#define EXIT(reason) { \
printf("Error: "); printf("%s ", reason); \
printf("Exiting from line %i in file %s\n",__LINE__,__FILE__); \
exit(-1); \
}







////////////////////////////






template <class T>
IntervalTreeNode<T>::IntervalTreeNode(){}

template <class T>
IntervalTreeNode<T>::IntervalTreeNode(Interval<T> * newInterval) 
  : storedInterval (newInterval) ,
    key(newInterval->GetLowPoint()), 
    high(newInterval->GetHighPoint()) , 
    maxHigh(high) {
}

template <class T>
IntervalTreeNode<T>::~IntervalTreeNode(){}

template <class T>
Interval<T>::Interval(){}

template <class T>
Interval<T>::~Interval(){}

template <class T>
void Interval<T>::Print() const {
  std::cout << "No Print Method defined for instance of Interval" << std::endl;
}

template <class T>
IntervalTree<T>::IntervalTree()
{
  nil = new IntervalTreeNode<T>;
  nil->left = nil->right = nil->parent = nil;
  nil->red = 0;
  nil->key = nil->high = nil->maxHigh = MIN_INT;
  nil->storedInterval = NULL;
  
  root = new IntervalTreeNode<T>;
  root->parent = root->left = root->right = nil;
  root->key = root->high = root->maxHigh = MAX_INT;
  root->red=0;
  root->storedInterval = NULL;

  /* the following are used for the Enumerate function */
  recursionNodeStackSize = 128;
  recursionNodeStack = (it_recursion_node<T> *) 
    malloc(recursionNodeStackSize * sizeof(it_recursion_node<T>));

  assert(recursionNodeStack);

  recursionNodeStackTop = 1;
  recursionNodeStack[0].start_node = NULL;
  
}

/***********************************************************************/
/*  FUNCTION:  LeftRotate */
/**/
/*  INPUTS:  the node to rotate on */
/**/
/*  OUTPUT:  None */
/**/
/*  Modifies Input: this, x */
/**/
/*  EFFECTS:  Rotates as described in _Introduction_To_Algorithms by */
/*            Cormen, Leiserson, Rivest (Chapter 14).  Basically this */
/*            makes the parent of x be to the left of x, x the parent of */
/*            its parent before the rotation and fixes other pointers */
/*            accordingly. Also updates the maxHigh fields of x and y */
/*            after rotation. */
/***********************************************************************/

template <class T>
void IntervalTree<T>::LeftRotate(IntervalTreeNode<T>* x) {
  IntervalTreeNode<T>* y;
 
  /*  I originally wrote this function to use the sentinel for */
  /*  nil to avoid checking for nil.  However this introduces a */
  /*  very subtle bug because sometimes this function modifies */
  /*  the parent pointer of nil.  This can be a problem if a */
  /*  function which calls LeftRotate also uses the nil sentinel */
  /*  and expects the nil sentinel's parent pointer to be unchanged */
  /*  after calling this function.  For example, when DeleteFixUP */
  /*  calls LeftRotate it expects the parent pointer of nil to be */
  /*  unchanged. */

  y=x->right;
  x->right=y->left;

  if (y->left != nil) y->left->parent=x; /* used to use sentinel here */
  /* and do an unconditional assignment instead of testing for nil */
  
  y->parent=x->parent;   

  /* instead of checking if x->parent is the root as in the book, we */
  /* count on the root sentinel to implicitly take care of this case */
  if( x == x->parent->left) {
    x->parent->left=y;
  } else {
    x->parent->right=y;
  }
  y->left=x;
  x->parent=y;

  x->maxHigh=ITMax(x->left->maxHigh,ITMax(x->right->maxHigh,x->high));
  y->maxHigh=ITMax(x->maxHigh,ITMax(y->right->maxHigh,y->high));
#ifdef CHECK_INTERVAL_TREE_ASSUMPTIONS
  CheckAssumptions();
#elif defined(DEBUG_ASSERT)
  Assert(!nil->red,"nil not red in ITLeftRotate");
  Assert((nil->maxHigh=MIN_INT),
   "nil->maxHigh != MIN_INT in ITLeftRotate");
#endif
}


/***********************************************************************/
/*  FUNCTION:  RighttRotate */
/**/
/*  INPUTS:  node to rotate on */
/**/
/*  OUTPUT:  None */
/**/
/*  Modifies Input?: this, y */
/**/
/*  EFFECTS:  Rotates as described in _Introduction_To_Algorithms by */
/*            Cormen, Leiserson, Rivest (Chapter 14).  Basically this */
/*            makes the parent of x be to the left of x, x the parent of */
/*            its parent before the rotation and fixes other pointers */
/*            accordingly. Also updates the maxHigh fields of x and y */
/*            after rotation. */
/***********************************************************************/

template <class T>
void IntervalTree<T>::RightRotate(IntervalTreeNode<T> *y) {
  IntervalTreeNode<T> *x;

  /*  I originally wrote this function to use the sentinel for */
  /*  nil to avoid checking for nil.  However this introduces a */
  /*  very subtle bug because sometimes this function modifies */
  /*  the parent pointer of nil.  This can be a problem if a */
  /*  function which calls LeftRotate also uses the nil sentinel */
  /*  and expects the nil sentinel's parent pointer to be unchanged */
  /*  after calling this function.  For example, when DeleteFixUP */
  /*  calls LeftRotate it expects the parent pointer of nil to be */
  /*  unchanged. */

  x=y->left;
  y->left=x->right;

  if (nil != x->right)  x->right->parent=y; /*used to use sentinel here */
  /* and do an unconditional assignment instead of testing for nil */

  /* instead of checking if x->parent is the root as in the book, we */
  /* count on the root sentinel to implicitly take care of this case */
  x->parent=y->parent;
  if( y == y->parent->left) {
    y->parent->left=x;
  } else {
    y->parent->right=x;
  }
  x->right=y;
  y->parent=x;

  y->maxHigh=ITMax(y->left->maxHigh,ITMax(y->right->maxHigh,y->high));
  x->maxHigh=ITMax(x->left->maxHigh,ITMax(y->maxHigh,x->high));
#ifdef CHECK_INTERVAL_TREE_ASSUMPTIONS
  CheckAssumptions();
#elif defined(DEBUG_ASSERT)
  Assert(!nil->red,"nil not red in ITRightRotate");
  Assert((nil->maxHigh=MIN_INT),
   "nil->maxHigh != MIN_INT in ITRightRotate");
#endif
}

/***********************************************************************/
/*  FUNCTION:  TreeInsertHelp  */
/**/
/*  INPUTS:  z is the node to insert */
/**/
/*  OUTPUT:  none */
/**/
/*  Modifies Input:  this, z */
/**/
/*  EFFECTS:  Inserts z into the tree as if it were a regular binary tree */
/*            using the algorithm described in _Introduction_To_Algorithms_ */
/*            by Cormen et al.  This funciton is only intended to be called */
/*            by the InsertTree function and not by the user */
/***********************************************************************/

template <class T>
void IntervalTree<T>::TreeInsertHelp(IntervalTreeNode<T> *z) {
  /*  This function should only be called by InsertITTree (see above) */
  IntervalTreeNode<T> *x;
  IntervalTreeNode<T> *y;
    
  z->left=z->right=nil;
  y=root;
  x=root->left;
  while( x != nil) {
    y=x;
    if ( x->key > z->key) { 
      x=x->left;
    } else { /* x->key <= z->key */
      x=x->right;
    }
  }
  z->parent=y;
  if ( (y == root) ||
       (y->key > z->key) ) { 
    y->left=z;
  } else {
    y->right=z;
  }


#if defined(DEBUG_ASSERT)
  Assert(!nil->red,"nil not red in ITTreeInsertHelp");
  Assert((nil->maxHigh=MIN_INT),
   "nil->maxHigh != MIN_INT in ITTreeInsertHelp");
#endif
}


/***********************************************************************/
/*  FUNCTION:  FixUpMaxHigh  */
/**/
/*  INPUTS:  x is the node to start from*/
/**/
/*  OUTPUT:  none */
/**/
/*  Modifies Input:  this */
/**/
/*  EFFECTS:  Travels up to the root fixing the maxHigh fields after */
/*            an insertion or deletion */
/***********************************************************************/

template <class T>
void IntervalTree<T>::FixUpMaxHigh(IntervalTreeNode<T> * x) {
  while(x != root) {
    x->maxHigh=ITMax(x->high,ITMax(x->left->maxHigh,x->right->maxHigh));
    x=x->parent;
  }
#ifdef CHECK_INTERVAL_TREE_ASSUMPTIONS
  CheckAssumptions();
#endif
}

/*  Before calling InsertNode  the node x should have its key set */

/***********************************************************************/
/*  FUNCTION:  InsertNode */
/**/
/*  INPUTS:  newInterval is the interval to insert*/
/**/
/*  OUTPUT:  This function returns a pointer to the newly inserted node */
/*           which is guarunteed to be valid until this node is deleted. */
/*           What this means is if another data structure stores this */
/*           pointer then the tree does not need to be searched when this */
/*           is to be deleted. */
/**/
/*  Modifies Input: tree */
/**/
/*  EFFECTS:  Creates a node node which contains the appropriate key and */
/*            info pointers and inserts it into the tree. */
/***********************************************************************/

template <class T>
IntervalTreeNode<T> * IntervalTree<T>::Insert(Interval<T> * newInterval)
{
  IntervalTreeNode<T> * y;
  IntervalTreeNode<T> * x;
  IntervalTreeNode<T> * newNode;

  x = new IntervalTreeNode<T>(newInterval);
  TreeInsertHelp(x);
  FixUpMaxHigh(x->parent);
  newNode = x;
  x->red=1;
  while(x->parent->red) { /* use sentinel instead of checking for root */
    if (x->parent == x->parent->parent->left) {
      y=x->parent->parent->right;
      if (y->red) {
  x->parent->red=0;
  y->red=0;
  x->parent->parent->red=1;
  x=x->parent->parent;
      } else {
  if (x == x->parent->right) {
    x=x->parent;
    LeftRotate(x);
  }
  x->parent->red=0;
  x->parent->parent->red=1;
  RightRotate(x->parent->parent);
      } 
    } else { /* case for x->parent == x->parent->parent->right */
             /* this part is just like the section above with */
             /* left and right interchanged */
      y=x->parent->parent->left;
      if (y->red) {
  x->parent->red=0;
  y->red=0;
  x->parent->parent->red=1;
  x=x->parent->parent;
      } else {
  if (x == x->parent->left) {
    x=x->parent;
    RightRotate(x);
  }
  x->parent->red=0;
  x->parent->parent->red=1;
  LeftRotate(x->parent->parent);
      } 
    }
  }
  root->left->red=0;
  return(newNode);

#ifdef CHECK_INTERVAL_TREE_ASSUMPTIONS
  CheckAssumptions();
#elif defined(DEBUG_ASSERT)
  Assert(!nil->red,"nil not red in ITTreeInsert");
  Assert(!root->red,"root not red in ITTreeInsert");
  Assert((nil->maxHigh=MIN_INT),
   "nil->maxHigh != MIN_INT in ITTreeInsert");
#endif
}

/***********************************************************************/
/*  FUNCTION:  GetSuccessorOf  */
/**/
/*    INPUTS:  x is the node we want the succesor of */
/**/
/*    OUTPUT:  This function returns the successor of x or NULL if no */
/*             successor exists. */
/**/
/*    Modifies Input: none */
/**/
/*    Note:  uses the algorithm in _Introduction_To_Algorithms_ */
/***********************************************************************/
  
template <class T>
IntervalTreeNode<T> * IntervalTree<T>::GetSuccessorOf(IntervalTreeNode<T> * x) const
{ 
  IntervalTreeNode<T> * y;

  if (nil != (y = x->right)) { /* assignment to y is intentional */
    while(y->left != nil) { /* returns the minium of the right subtree of x */
      y=y->left;
    }
    return(y);
  } else {
    y=x->parent;
    while(x == y->right) { /* sentinel used instead of checking for nil */
      x=y;
      y=y->parent;
    }
    if (y == root) return(nil);
    return(y);
  }
}

/***********************************************************************/
/*  FUNCTION:  GetPredecessorOf  */
/**/
/*    INPUTS:  x is the node to get predecessor of */
/**/
/*    OUTPUT:  This function returns the predecessor of x or NULL if no */
/*             predecessor exists. */
/**/
/*    Modifies Input: none */
/**/
/*    Note:  uses the algorithm in _Introduction_To_Algorithms_ */
/***********************************************************************/

template <class T>
IntervalTreeNode<T> *IntervalTree<T>::GetPredecessorOf(IntervalTreeNode<T> * x) const {
  IntervalTreeNode<T> *y;

  if (nil != (y = x->left)) { /* assignment to y is intentional */
    while(y->right != nil) { /* returns the maximum of the left subtree of x */
      y=y->right;
    }
    return(y);
  } else {
    y=x->parent;
    while(x == y->left) { 
      if (y == root) return(nil); 
      x=y;
      y=y->parent;
    }
    return(y);
  }
}

/***********************************************************************/
/*  FUNCTION:  Print */
/**/
/*    INPUTS:  none */
/**/
/*    OUTPUT:  none  */
/**/
/*    EFFECTS:  This function recursively prints the nodes of the tree */
/*              inorder. */
/**/
/*    Modifies Input: none */
/**/
/*    Note:    This function should only be called from ITTreePrint */
/***********************************************************************/

template <class T>
void IntervalTreeNode<T>::Print(IntervalTreeNode<T> * nil,
           IntervalTreeNode<T> * root) const {
  storedInterval->Print();
  printf(", k=%i, h=%i, mH=%i",key,high,maxHigh);
  printf("  l->key=");
  if( left == nil) printf("NULL"); else printf("%i",left->key);
  printf("  r->key=");
  if( right == nil) printf("NULL"); else printf("%i",right->key);
  printf("  p->key=");
  if( parent == root) printf("NULL"); else printf("%i",parent->key);
  printf("  red=%i\n",red);
}

template <class T>
void IntervalTree<T>::TreePrintHelper(IntervalTreeNode<T> *x) const {
  
  if (x != nil) {
    TreePrintHelper(x->left);
    x->Print(nil,root);
    TreePrintHelper(x->right);
  }
}

template <class T>
IntervalTree<T>::~IntervalTree() {
  IntervalTreeNode<T> * x = root->left;
   std::stack<IntervalTreeNode<T> *> stuffToFree;

  if (x != nil) {
    if (x->left != nil) {
      stuffToFree.push(x->left);
    }
    if (x->right != nil) {
      stuffToFree.push(x->right);
    }
    // delete x->storedInterval;
    delete x;
    while( !stuffToFree.empty() ) {
      x = stuffToFree.top();
      stuffToFree.pop();
      if (x->left != nil) {
  stuffToFree.push(x->left);
      }
      if (x->right != nil) {
  stuffToFree.push(x->right);
      }
      // delete x->storedInterval;
      delete x;
    }
  }
  delete nil;
  delete root;
  free(recursionNodeStack);
}


/***********************************************************************/
/*  FUNCTION:  Print */
/**/
/*    INPUTS:  none */
/**/
/*    OUTPUT:  none */
/**/
/*    EFFECT:  This function recursively prints the nodes of the tree */
/*             inorder. */
/**/
/*    Modifies Input: none */
/**/
/***********************************************************************/

template <class T>
void IntervalTree<T>::Print() const {
  TreePrintHelper(root->left);
}

/***********************************************************************/
/*  FUNCTION:  DeleteFixUp */
/**/
/*    INPUTS:  x is the child of the spliced */
/*             out node in DeleteNode. */
/**/
/*    OUTPUT:  none */
/**/
/*    EFFECT:  Performs rotations and changes colors to restore red-black */
/*             properties after a node is deleted */
/**/
/*    Modifies Input: this, x */
/**/
/*    The algorithm from this function is from _Introduction_To_Algorithms_ */
/***********************************************************************/

template <class T>
void IntervalTree<T>::DeleteFixUp(IntervalTreeNode<T> *x) {
  IntervalTreeNode<T> * w;
  IntervalTreeNode<T> * rootLeft = root->left;

  while( (!x->red) && (rootLeft != x)) {
    if (x == x->parent->left) {
      w=x->parent->right;
      if (w->red) {
  w->red=0;
  x->parent->red=1;
  LeftRotate(x->parent);
  w=x->parent->right;
      }
      if ( (!w->right->red) && (!w->left->red) ) { 
  w->red=1;
  x=x->parent;
      } else {
  if (!w->right->red) {
    w->left->red=0;
    w->red=1;
    RightRotate(w);
    w=x->parent->right;
  }
  w->red=x->parent->red;
  x->parent->red=0;
  w->right->red=0;
  LeftRotate(x->parent);
  x=rootLeft; /* this is to exit while loop */
      }
    } else { /* the code below is has left and right switched from above */
      w=x->parent->left;
      if (w->red) {
  w->red=0;
  x->parent->red=1;
  RightRotate(x->parent);
  w=x->parent->left;
      }
      if ( (!w->right->red) && (!w->left->red) ) { 
  w->red=1;
  x=x->parent;
      } else {
  if (!w->left->red) {
    w->right->red=0;
    w->red=1;
    LeftRotate(w);
    w=x->parent->left;
  }
  w->red=x->parent->red;
  x->parent->red=0;
  w->left->red=0;
  RightRotate(x->parent);
  x=rootLeft; /* this is to exit while loop */
      }
    }
  }
  x->red=0;

#ifdef CHECK_INTERVAL_TREE_ASSUMPTIONS
  CheckAssumptions();
#elif defined(DEBUG_ASSERT)
  Assert(!nil->red,"nil not black in ITDeleteFixUp");
  Assert((nil->maxHigh=MIN_INT),
   "nil->maxHigh != MIN_INT in ITDeleteFixUp");
#endif
}


/***********************************************************************/
/*  FUNCTION:  DeleteNode */
/**/
/*    INPUTS:  tree is the tree to delete node z from */
/**/
/*    OUTPUT:  returns the Interval stored at deleted node */
/**/
/*    EFFECT:  Deletes z from tree and but don't call destructor */
/*             Then calls FixUpMaxHigh to fix maxHigh fields then calls */
/*             ITDeleteFixUp to restore red-black properties */
/**/
/*    Modifies Input:  z */
/**/
/*    The algorithm from this function is from _Introduction_To_Algorithms_ */
/***********************************************************************/

template <class T>
Interval<T> * IntervalTree<T>::DeleteNode(IntervalTreeNode<T> * z){
  IntervalTreeNode<T> * y;
  IntervalTreeNode<T> * x;
  Interval<T> * returnValue = z->storedInterval;

  y= ((z->left == nil) || (z->right == nil)) ? z : GetSuccessorOf(z);
  x= (y->left == nil) ? y->right : y->left;
  if (root == (x->parent = y->parent)) { /* assignment of y->p to x->p is intentional */
    root->left=x;
  } else {
    if (y == y->parent->left) {
      y->parent->left=x;
    } else {
      y->parent->right=x;
    }
  }
  if (y != z) { /* y should not be nil in this case */

#ifdef DEBUG_ASSERT
    Assert( (y!=nil),"y is nil in DeleteNode \n");
#endif
    /* y is the node to splice out and x is its child */
  
    y->maxHigh = MIN_INT;
    y->left=z->left;
    y->right=z->right;
    y->parent=z->parent;
    z->left->parent=z->right->parent=y;
    if (z == z->parent->left) {
      z->parent->left=y; 
    } else {
      z->parent->right=y;
    }
    FixUpMaxHigh(x->parent); 
    if (!(y->red)) {
      y->red = z->red;
      DeleteFixUp(x);
    } else
      y->red = z->red; 
    delete z;
#ifdef CHECK_INTERVAL_TREE_ASSUMPTIONS
    CheckAssumptions();
#elif defined(DEBUG_ASSERT)
    Assert(!nil->red,"nil not black in ITDelete");
    Assert((nil->maxHigh=MIN_INT),"nil->maxHigh != MIN_INT in ITDelete");
#endif
  } else {
    FixUpMaxHigh(x->parent);
    if (!(y->red)) DeleteFixUp(x);
    delete y;
#ifdef CHECK_INTERVAL_TREE_ASSUMPTIONS
    CheckAssumptions();
#elif defined(DEBUG_ASSERT)
    Assert(!nil->red,"nil not black in ITDelete");
    Assert((nil->maxHigh=MIN_INT),"nil->maxHigh != MIN_INT in ITDelete");
#endif
  }
  return returnValue;
}


/***********************************************************************/
/*  FUNCTION:  Overlap */
/**/
/*    INPUTS:  [a1,a2] and [b1,b2] are the low and high endpoints of two */
/*             closed intervals.  */
/**/
/*    OUTPUT:  stack containing pointers to the nodes between [low,high] */
/**/
/*    Modifies Input: none */
/**/
/*    EFFECT:  returns 1 if the intervals overlap, and 0 otherwise */
/***********************************************************************/

int Overlap(int a1, int a2, int b1, int b2) {
  if (a1 <= b1) {
    return( (b1 <= a2) );
  } else {
    return( (a1 <= b2) );
  }
}


/***********************************************************************/
/*  FUNCTION:  Enumerate */
/**/
/*    INPUTS:  tree is the tree to look for intervals overlapping the */
/*             closed interval [low,high]  */
/**/
/*    OUTPUT:  stack containing pointers to the nodes overlapping */
/*             [low,high] */
/**/
/*    Modifies Input: none */
/**/
/*    EFFECT:  Returns a stack containing pointers to nodes containing */
/*             intervals which overlap [low,high] in O(max(N,k*log(N))) */
/*             where N is the number of intervals in the tree and k is  */
/*             the number of overlapping intervals                      */
/**/
/*    Note:    This basic idea for this function comes from the  */
/*              _Introduction_To_Algorithms_ book by Cormen et al, but */
/*             modifications were made to return all overlapping intervals */
/*             instead of just the first overlapping interval as in the */
/*             book.  The natural way to do this would require recursive */
/*             calls of a basic search function.  I translated the */
/*             recursive version into an interative version with a stack */
/*             as described below. */
/***********************************************************************/



/*  The basic idea for the function below is to take the IntervalSearch */
/*  function from the book and modify to find all overlapping intervals */
/*  instead of just one.  This means that any time we take the left */
/*  branch down the tree we must also check the right branch if and only if */
/*  we find an overlapping interval in that left branch.  Note this is a */
/*  recursive condition because if we go left at the root then go left */
/*  again at the first left child and find an overlap in the left subtree */
/*  of the left child of root we must recursively check the right subtree */
/*  of the left child of root as well as the right child of root. */

template <class T>
std::queue<void *> * IntervalTree<T>::Enumerate(T low, 
              T high)  {
  std::queue<void *> * enumResultStack;
  IntervalTreeNode<T> * x=root->left;
  int stuffToDo = (x != nil);
  
  // Possible speed up: add min field to prune right searches //

#ifdef DEBUG_ASSERT
  Assert((recursionNodeStackTop == 1),
   "recursionStack not empty when entering IntervalTree::Enumerate");
#endif
  currentParent = 0;
  enumResultStack = new std::queue<void *>;

  while(stuffToDo) {
    if (Overlap(low,high,x->key,x->high) ) {
      enumResultStack->push(x->storedInterval);
      recursionNodeStack[currentParent].tryRightBranch=1;
    }
    if(x->left->maxHigh >= low) { // implies x != nil 
      if ( recursionNodeStackTop == recursionNodeStackSize ) {
  recursionNodeStackSize *= 2;
  recursionNodeStack = (it_recursion_node<T> *) 
    realloc(recursionNodeStack,
      recursionNodeStackSize * sizeof(it_recursion_node<T>));
  if (recursionNodeStack == NULL) 
    EXIT("realloc failed in IntervalTree::Enumerate\n");
      }
      recursionNodeStack[recursionNodeStackTop].start_node = x;
      recursionNodeStack[recursionNodeStackTop].tryRightBranch = 0;
      recursionNodeStack[recursionNodeStackTop].parentIndex = currentParent;
      currentParent = recursionNodeStackTop++;
      x = x->left;
    } else {
      x = x->right;
    }
    stuffToDo = (x != nil);
    while( (!stuffToDo) && (recursionNodeStackTop > 1) ) {
  if(recursionNodeStack[--recursionNodeStackTop].tryRightBranch) {
    x=recursionNodeStack[recursionNodeStackTop].start_node->right;
    currentParent=recursionNodeStack[recursionNodeStackTop].parentIndex;
    recursionNodeStack[currentParent].tryRightBranch=1;
    stuffToDo = ( x != nil);
  }
    }
  }
#ifdef DEBUG_ASSERT
  Assert((recursionNodeStackTop == 1),
   "recursionStack not empty when exiting IntervalTree::Enumerate");
#endif
  return(enumResultStack);   
}
  

template <class T>
int IntervalTree<T>::CheckMaxHighFieldsHelper(IntervalTreeNode<T> * y, 
            const T currentHigh,
            int match) const
{
  if (y != nil) {
    match = CheckMaxHighFieldsHelper(y->left,currentHigh,match) ?
      1 : match;
    VERIFY(y->high <= currentHigh);
    if (y->high == currentHigh)
      match = 1;
    match = CheckMaxHighFieldsHelper(y->right,currentHigh,match) ?
      1 : match;
  }
  return match;
}

    

/* Make sure the maxHigh fields for everything makes sense. *
 * If something is wrong, print a warning and exit */
template <class T>
void IntervalTree<T>::CheckMaxHighFields(IntervalTreeNode<T> * x) const {
  if (x != nil) {
    CheckMaxHighFields(x->left);
    if(!(CheckMaxHighFieldsHelper(x,x->maxHigh,0) > 0)) {
      EXIT("error found in CheckMaxHighFields.\n");
    }
    CheckMaxHighFields(x->right);
  }
}

template <class T>
void IntervalTree<T>::CheckAssumptions() const {
 VERIFY(nil->key == MIN_INT);
 VERIFY(nil->high == MIN_INT);
 VERIFY(nil->maxHigh == MIN_INT);
 VERIFY(root->key == MAX_INT);
 VERIFY(root->high == MAX_INT);
 VERIFY(root->maxHigh == MAX_INT);
 VERIFY(nil->storedInterval == NULL);
 VERIFY(root->storedInterval == NULL);
 VERIFY(nil->red == 0);
 VERIFY(root->red == 0);
 CheckMaxHighFields(root->left);
}


 




#endif



