//
// File generated by rootcint at Fri Apr 25 16:41:23 2014

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME DictOutput
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "DictOutput.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
      #if !(defined(R__ACCESS_IN_SYMBOL) || defined(R__USE_SHADOW_CLASS))
      typedef ::SimpleClass SimpleClass;
      #else
      class SimpleClass  {
         public:
         //friend XX;
         int an_obvious_variable; //
         int a_secret_variable; //
      };
      #endif

      #if !(defined(R__ACCESS_IN_SYMBOL) || defined(R__USE_SHADOW_CLASS))
      typedef ::runevtobj runevtobj;
      #else
      class runevtobj  {
         public:
         //friend XX;
         int run_num; //
         int evt_num; //
         bool is_lowest; //
      };
      #endif

   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void SimpleClass_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void SimpleClass_Dictionary();
   static void *new_SimpleClass(void *p = 0);
   static void *newArray_SimpleClass(Long_t size, void *p);
   static void delete_SimpleClass(void *p);
   static void deleteArray_SimpleClass(void *p);
   static void destruct_SimpleClass(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimpleClass*)
   {
      // Make sure the shadow class has the right sizeof
      R__ASSERT(sizeof(::SimpleClass) == sizeof(::ROOT::Shadow::SimpleClass));
      ::SimpleClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SimpleClass),0);
      static ::ROOT::TGenericClassInfo 
         instance("SimpleClass", "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/SimpleClass.h", 4,
                  typeid(::SimpleClass), DefineBehavior(ptr, ptr),
                  &SimpleClass_ShowMembers, &SimpleClass_Dictionary, isa_proxy, 4,
                  sizeof(::SimpleClass) );
      instance.SetNew(&new_SimpleClass);
      instance.SetNewArray(&newArray_SimpleClass);
      instance.SetDelete(&delete_SimpleClass);
      instance.SetDeleteArray(&deleteArray_SimpleClass);
      instance.SetDestructor(&destruct_SimpleClass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimpleClass*)
   {
      return GenerateInitInstanceLocal((::SimpleClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::SimpleClass*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void SimpleClass_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::SimpleClass*)0x0)->GetClass();
   }

} // end of namespace ROOT

namespace ROOT {
   void runevtobj_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void runevtobj_Dictionary();
   static void *new_runevtobj(void *p = 0);
   static void *newArray_runevtobj(Long_t size, void *p);
   static void delete_runevtobj(void *p);
   static void deleteArray_runevtobj(void *p);
   static void destruct_runevtobj(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::runevtobj*)
   {
      // Make sure the shadow class has the right sizeof
      R__ASSERT(sizeof(::runevtobj) == sizeof(::ROOT::Shadow::runevtobj));
      ::runevtobj *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::runevtobj),0);
      static ::ROOT::TGenericClassInfo 
         instance("runevtobj", "/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/runevtobj.h", 5,
                  typeid(::runevtobj), DefineBehavior(ptr, ptr),
                  &runevtobj_ShowMembers, &runevtobj_Dictionary, isa_proxy, 4,
                  sizeof(::runevtobj) );
      instance.SetNew(&new_runevtobj);
      instance.SetNewArray(&newArray_runevtobj);
      instance.SetDelete(&delete_runevtobj);
      instance.SetDeleteArray(&deleteArray_runevtobj);
      instance.SetDestructor(&destruct_runevtobj);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::runevtobj*)
   {
      return GenerateInitInstanceLocal((::runevtobj*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::runevtobj*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void runevtobj_Dictionary() {
      ::ROOT::GenerateInitInstanceLocal((const ::runevtobj*)0x0)->GetClass();
   }

} // end of namespace ROOT

//______________________________________________________________________________
namespace ROOT {
   void SimpleClass_ShowMembers(void *obj, TMemberInspector &R__insp)
   {
      // Inspect the data members of an object of class SimpleClass.
      typedef ::ROOT::Shadow::SimpleClass ShadowClass;
      ShadowClass *sobj = (ShadowClass*)obj;
      if (sobj) { } // Dummy usage just in case there is no datamember.

      TClass *R__cl  = ::ROOT::GenerateInitInstanceLocal((const ::SimpleClass*)0x0)->GetClass();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "an_obvious_variable", &sobj->an_obvious_variable);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "a_secret_variable", &sobj->a_secret_variable);
   }

}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimpleClass(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::SimpleClass : new ::SimpleClass;
   }
   static void *newArray_SimpleClass(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::SimpleClass[nElements] : new ::SimpleClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimpleClass(void *p) {
      delete ((::SimpleClass*)p);
   }
   static void deleteArray_SimpleClass(void *p) {
      delete [] ((::SimpleClass*)p);
   }
   static void destruct_SimpleClass(void *p) {
      typedef ::SimpleClass current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimpleClass

//______________________________________________________________________________
namespace ROOT {
   void runevtobj_ShowMembers(void *obj, TMemberInspector &R__insp)
   {
      // Inspect the data members of an object of class runevtobj.
      typedef ::ROOT::Shadow::runevtobj ShadowClass;
      ShadowClass *sobj = (ShadowClass*)obj;
      if (sobj) { } // Dummy usage just in case there is no datamember.

      TClass *R__cl  = ::ROOT::GenerateInitInstanceLocal((const ::runevtobj*)0x0)->GetClass();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "run_num", &sobj->run_num);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "evt_num", &sobj->evt_num);
      R__insp.Inspect(R__cl, R__insp.GetParent(), "is_lowest", &sobj->is_lowest);
   }

}

namespace ROOT {
   // Wrappers around operator new
   static void *new_runevtobj(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::runevtobj : new ::runevtobj;
   }
   static void *newArray_runevtobj(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::runevtobj[nElements] : new ::runevtobj[nElements];
   }
   // Wrapper around operator delete
   static void delete_runevtobj(void *p) {
      delete ((::runevtobj*)p);
   }
   static void deleteArray_runevtobj(void *p) {
      delete [] ((::runevtobj*)p);
   }
   static void destruct_runevtobj(void *p) {
      typedef ::runevtobj current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::runevtobj

/********************************************************
* DictOutput.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableDictOutput();

extern "C" void G__set_cpp_environmentDictOutput() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/SimpleClass.h");
  G__add_compiledheader("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/events_comparison.h");
  G__add_compiledheader("/direct/phenix+u/workarea/danielj/cvs_code/offline/analysis/danielj/run13_analysis/utility/cross_check/events_comparison/runevtobj.h");
  G__cpp_reset_tagtableDictOutput();
}
#include <new>
extern "C" int G__cpp_dllrevDictOutput() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* SimpleClass */
static int G__DictOutput_162_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   SimpleClass* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new SimpleClass[n];
     } else {
       p = new((void*) gvp) SimpleClass[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new SimpleClass;
     } else {
       p = new((void*) gvp) SimpleClass;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DictOutputLN_SimpleClass));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DictOutput_162_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) ((SimpleClass*) G__getstructoffset())->SetSecret((int) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__DictOutput_162_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   SimpleClass* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new SimpleClass(*(SimpleClass*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DictOutputLN_SimpleClass));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef SimpleClass G__TSimpleClass;
static int G__DictOutput_162_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 0
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (SimpleClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((SimpleClass*) (soff+(sizeof(SimpleClass)*i)))->~G__TSimpleClass();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (SimpleClass*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((SimpleClass*) (soff))->~G__TSimpleClass();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__DictOutput_162_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   SimpleClass* dest = (SimpleClass*) G__getstructoffset();
   *dest = *(SimpleClass*) libp->para[0].ref;
   const SimpleClass& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* runevtobj */
static int G__DictOutput_163_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   runevtobj* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new runevtobj[n];
     } else {
       p = new((void*) gvp) runevtobj[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new runevtobj;
     } else {
       p = new((void*) gvp) runevtobj;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DictOutputLN_runevtobj));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__DictOutput_163_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   runevtobj* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 2
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new runevtobj((int) G__int(libp->para[0]), (int) G__int(libp->para[1]));
   } else {
     p = new((void*) gvp) runevtobj((int) G__int(libp->para[0]), (int) G__int(libp->para[1]));
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DictOutputLN_runevtobj));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__DictOutput_163_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   runevtobj* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new runevtobj(*(runevtobj*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__DictOutputLN_runevtobj));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef runevtobj G__Trunevtobj;
static int G__DictOutput_163_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 0
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (runevtobj*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((runevtobj*) (soff+(sizeof(runevtobj)*i)))->~G__Trunevtobj();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (runevtobj*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((runevtobj*) (soff))->~G__Trunevtobj();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__DictOutput_163_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   runevtobj* dest = (runevtobj*) G__getstructoffset();
   *dest = *(runevtobj*) libp->para[0].ref;
   const runevtobj& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */
static int G__DictOutput__0_1290(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      events_comparison();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/*********************************************************
* Member function Stub
*********************************************************/

/* SimpleClass */

/* runevtobj */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncDictOutput {
 public:
  G__Sizep2memfuncDictOutput(): p(&G__Sizep2memfuncDictOutput::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncDictOutput::*p)();
};

size_t G__get_sizep2memfuncDictOutput()
{
  G__Sizep2memfuncDictOutput a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceDictOutput() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableDictOutput() {

   /* Setting up typedef entry */
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__DictOutputLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DictOutputLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DictOutputLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DictOutputLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DictOutputLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__DictOutputLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__DictOutputLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DictOutputLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__DictOutputLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__DictOutputLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* SimpleClass */
static void G__setup_memvarSimpleClass(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__DictOutputLN_SimpleClass));
   { SimpleClass *p; p=(SimpleClass*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->an_obvious_variable)-(long)(p)),105,0,0,-1,-1,-1,1,"an_obvious_variable=",0,(char*)NULL);
   G__memvar_setup((void*)0,105,0,0,-1,-1,-1,4,"a_secret_variable=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}


   /* runevtobj */
static void G__setup_memvarrunevtobj(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__DictOutputLN_runevtobj));
   { runevtobj *p; p=(runevtobj*)0x1000; if (p) { }
   G__memvar_setup((void*)((long)(&p->run_num)-(long)(p)),105,0,0,-1,-1,-1,1,"run_num=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->evt_num)-(long)(p)),105,0,0,-1,-1,-1,1,"evt_num=",0,(char*)NULL);
   G__memvar_setup((void*)((long)(&p->is_lowest)-(long)(p)),103,0,0,-1,-1,-1,1,"is_lowest=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarDictOutput() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncSimpleClass(void) {
   /* SimpleClass */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__DictOutputLN_SimpleClass));
   G__memfunc_setup("SimpleClass",1120,G__DictOutput_162_0_1, 105, G__get_linked_tagnum(&G__DictOutputLN_SimpleClass), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("SetSecret",914,G__DictOutput_162_0_2, 105, -1, -1, 0, 1, 1, 1, 0, "i - - 0 - -", (char*)NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("SimpleClass", 1120, G__DictOutput_162_0_3, (int) ('i'), G__get_linked_tagnum(&G__DictOutputLN_SimpleClass), -1, 0, 1, 1, 1, 0, "u 'SimpleClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~SimpleClass", 1246, G__DictOutput_162_0_4, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__DictOutput_162_0_5, (int) ('u'), G__get_linked_tagnum(&G__DictOutputLN_SimpleClass), -1, 1, 1, 1, 1, 0, "u 'SimpleClass' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}

static void G__setup_memfuncrunevtobj(void) {
   /* runevtobj */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__DictOutputLN_runevtobj));
   G__memfunc_setup("runevtobj",991,G__DictOutput_163_0_1, 105, G__get_linked_tagnum(&G__DictOutputLN_runevtobj), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("runevtobj",991,G__DictOutput_163_0_2, 105, G__get_linked_tagnum(&G__DictOutputLN_runevtobj), -1, 0, 2, 1, 1, 0, 
"i - - 0 - run i - - 0 - evt", (char*)NULL, (void*) NULL, 0);
   // automatic copy constructor
   G__memfunc_setup("runevtobj", 991, G__DictOutput_163_0_3, (int) ('i'), G__get_linked_tagnum(&G__DictOutputLN_runevtobj), -1, 0, 1, 1, 1, 0, "u 'runevtobj' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~runevtobj", 1117, G__DictOutput_163_0_4, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__DictOutput_163_0_5, (int) ('u'), G__get_linked_tagnum(&G__DictOutputLN_runevtobj), -1, 1, 1, 1, 1, 0, "u 'runevtobj' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncDictOutput() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalDictOutput() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {
   G__memfunc_setup("events_comparison", 1839, G__DictOutput__0_1290, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL
, (void*) NULL, 0);

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcDictOutput() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__DictOutputLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__DictOutputLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DictOutputLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__DictOutputLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__DictOutputLN_SimpleClass = { "SimpleClass" , 99 , -1 };
G__linked_taginfo G__DictOutputLN_runevtobj = { "runevtobj" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableDictOutput() {
  G__DictOutputLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__DictOutputLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DictOutputLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__DictOutputLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__DictOutputLN_SimpleClass.tagnum = -1 ;
  G__DictOutputLN_runevtobj.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableDictOutput() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__DictOutputLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__DictOutputLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__DictOutputLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__DictOutputLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__DictOutputLN_SimpleClass),sizeof(SimpleClass),-1,263424,(char*)NULL,G__setup_memvarSimpleClass,G__setup_memfuncSimpleClass);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__DictOutputLN_runevtobj),sizeof(runevtobj),-1,296192,(char*)NULL,G__setup_memvarrunevtobj,G__setup_memfuncrunevtobj);
}
extern "C" void G__cpp_setupDictOutput(void) {
  G__check_setup_version(30051515,"G__cpp_setupDictOutput()");
  G__set_cpp_environmentDictOutput();
  G__cpp_setup_tagtableDictOutput();

  G__cpp_setup_inheritanceDictOutput();

  G__cpp_setup_typetableDictOutput();

  G__cpp_setup_memvarDictOutput();

  G__cpp_setup_memfuncDictOutput();
  G__cpp_setup_globalDictOutput();
  G__cpp_setup_funcDictOutput();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncDictOutput();
  return;
}
class G__cpp_setup_initDictOutput {
  public:
    G__cpp_setup_initDictOutput() { G__add_setup_func("DictOutput",(G__incsetup)(&G__cpp_setupDictOutput)); G__call_setup_funcs(); }
   ~G__cpp_setup_initDictOutput() { G__remove_setup_func("DictOutput"); }
};
G__cpp_setup_initDictOutput G__cpp_setup_initializerDictOutput;

