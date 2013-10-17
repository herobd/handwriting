#include "dinstancecounter.h"


/*! \class DInstanceCounter
    \brief Provides debug functions for tracking object instance counts

    When compiling with DEBUG defined, the static functions of this class
    are used to count how many instances of a given object have been
    constructed and destructed, and allows an application to print out a
    report.  This is useful for debugging to make sure the application is
    cleaning up after itself.  When DEBUG is not defined, the functions for
    adding and removing instances (upon construction/destruction) become
    empty functions that should be optimized away by the compiler so that
    release code is not slow.  Calls to addInstance() should be made from
    within constructors, and calls to removeInstance should be made from
    within the destructor.  For example:

    \code
      #include "dinstancecounter.h"

      class Foo{
        public:
        Foo();
        Foo(const Foo &src);
        ~Foo();
      };
      inline Foo::Foo(){DInstanceCounter::addInstance("Foo");...}
      inline Foo::Foo(const Foo &src){DInstanceCounter::addInstance("Foo");...}
      inline Foo::~Foo(){DInstanceCounter::removeInstance("Foo");...}
    \endcode

    \warning Default constructors and copy constructors should be over-ridden so that
    all object instantiations are counted.
**/



// static member variables:

std::map<const char*, int/*, DInstanceCounter::lstr*/> 
  DInstanceCounter::mapNumConstructed;
std::map<const char*, int/*, DInstanceCounter::lstr*/> 
  DInstanceCounter::mapNumDestructed;


/// Clean up the DInstanceCounter memory.
/** DInstanceCounter stores the class name and count of constructed/destructed
    objects.  This function is used by the application to release that memory
    when it is no longer in use.
**/
void DInstanceCounter::cleanup(){
  mapNumConstructed.clear();
  mapNumDestructed.clear();
}


