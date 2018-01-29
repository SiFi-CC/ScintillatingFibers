// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__ScintillatingFibers

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/kasia/ScintillatingFibers/sources/include/SFData.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_SFData(void *p = 0);
   static void *newArray_SFData(Long_t size, void *p);
   static void delete_SFData(void *p);
   static void deleteArray_SFData(void *p);
   static void destruct_SFData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SFData*)
   {
      ::SFData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SFData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SFData", ::SFData::Class_Version(), "SFData.hh", 18,
                  typeid(::SFData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SFData::Dictionary, isa_proxy, 4,
                  sizeof(::SFData) );
      instance.SetNew(&new_SFData);
      instance.SetNewArray(&newArray_SFData);
      instance.SetDelete(&delete_SFData);
      instance.SetDeleteArray(&deleteArray_SFData);
      instance.SetDestructor(&destruct_SFData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SFData*)
   {
      return GenerateInitInstanceLocal((::SFData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SFData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr SFData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SFData::Class_Name()
{
   return "SFData";
}

//______________________________________________________________________________
const char *SFData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SFData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SFData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SFData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SFData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SFData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SFData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SFData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void SFData::Streamer(TBuffer &R__b)
{
   // Stream an object of class SFData.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SFData::Class(),this);
   } else {
      R__b.WriteClassBuffer(SFData::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SFData(void *p) {
      return  p ? new(p) ::SFData : new ::SFData;
   }
   static void *newArray_SFData(Long_t nElements, void *p) {
      return p ? new(p) ::SFData[nElements] : new ::SFData[nElements];
   }
   // Wrapper around operator delete
   static void delete_SFData(void *p) {
      delete ((::SFData*)p);
   }
   static void deleteArray_SFData(void *p) {
      delete [] ((::SFData*)p);
   }
   static void destruct_SFData(void *p) {
      typedef ::SFData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SFData

namespace {
  void TriggerDictionaryInitialization_libScintillatingFibers_Impl() {
    static const char* headers[] = {
"/home/kasia/ScintillatingFibers/sources/include/SFData.hh",
0
    };
    static const char* includePaths[] = {
"/home/kasia/Install/root-git/root/build/include",
"/home/kasia/ScintillatingFibers/sources/include",
"/home/kasia/Install/root-git/root/build/include",
"/home/kasia/ScintillatingFibers/build/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libScintillatingFibers dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$/home/kasia/ScintillatingFibers/sources/include/SFData.hh")))  SFData;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libScintillatingFibers dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "/home/kasia/ScintillatingFibers/sources/include/SFData.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"SFData", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libScintillatingFibers",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libScintillatingFibers_Impl, {}, classesHeaders, /*has C++ module?*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libScintillatingFibers_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libScintillatingFibers() {
  TriggerDictionaryInitialization_libScintillatingFibers_Impl();
}
