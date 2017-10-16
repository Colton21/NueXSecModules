// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME xsecAna_dict

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
#include "TpcObjectContainer.h"
#include "ParticleContainer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
static TClass *xsecAnacLcLParticleContainer_Dictionary();
static void xsecAnacLcLParticleContainer_TClassManip(TClass*);
static void *new_xsecAnacLcLParticleContainer(void *p = 0);
static void *newArray_xsecAnacLcLParticleContainer(Long_t size, void *p);
static void delete_xsecAnacLcLParticleContainer(void *p);
static void deleteArray_xsecAnacLcLParticleContainer(void *p);
static void destruct_xsecAnacLcLParticleContainer(void *p);

// Function generating the singleton type initializer
static TGenericClassInfo *GenerateInitInstanceLocal(const ::xsecAna::ParticleContainer*)
{
	::xsecAna::ParticleContainer *ptr = 0;
	static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::xsecAna::ParticleContainer));
	static ::ROOT::TGenericClassInfo
	        instance("xsecAna::ParticleContainer", "ParticleContainer.h", 10,
	                 typeid(::xsecAna::ParticleContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
	                 &xsecAnacLcLParticleContainer_Dictionary, isa_proxy, 0,
	                 sizeof(::xsecAna::ParticleContainer) );
	instance.SetNew(&new_xsecAnacLcLParticleContainer);
	instance.SetNewArray(&newArray_xsecAnacLcLParticleContainer);
	instance.SetDelete(&delete_xsecAnacLcLParticleContainer);
	instance.SetDeleteArray(&deleteArray_xsecAnacLcLParticleContainer);
	instance.SetDestructor(&destruct_xsecAnacLcLParticleContainer);
	return &instance;
}
TGenericClassInfo *GenerateInitInstance(const ::xsecAna::ParticleContainer*)
{
	return GenerateInitInstanceLocal((::xsecAna::ParticleContainer*)0);
}
// Static variable to force the class initialization
static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::xsecAna::ParticleContainer*)0x0); R__UseDummy(_R__UNIQUE_(Init));

// Dictionary for non-ClassDef classes
static TClass *xsecAnacLcLParticleContainer_Dictionary() {
	TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::xsecAna::ParticleContainer*)0x0)->GetClass();
	xsecAnacLcLParticleContainer_TClassManip(theClass);
	return theClass;
}

static void xsecAnacLcLParticleContainer_TClassManip(TClass* ){
}

} // end of namespace ROOT

namespace ROOT {
static TClass *xsecAnacLcLTPCObjectContainer_Dictionary();
static void xsecAnacLcLTPCObjectContainer_TClassManip(TClass*);
static void *new_xsecAnacLcLTPCObjectContainer(void *p = 0);
static void *newArray_xsecAnacLcLTPCObjectContainer(Long_t size, void *p);
static void delete_xsecAnacLcLTPCObjectContainer(void *p);
static void deleteArray_xsecAnacLcLTPCObjectContainer(void *p);
static void destruct_xsecAnacLcLTPCObjectContainer(void *p);

// Function generating the singleton type initializer
static TGenericClassInfo *GenerateInitInstanceLocal(const ::xsecAna::TPCObjectContainer*)
{
	::xsecAna::TPCObjectContainer *ptr = 0;
	static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::xsecAna::TPCObjectContainer));
	static ::ROOT::TGenericClassInfo
	        instance("xsecAna::TPCObjectContainer", "TpcObjectContainer.h", 12,
	                 typeid(::xsecAna::TPCObjectContainer), ::ROOT::Internal::DefineBehavior(ptr, ptr),
	                 &xsecAnacLcLTPCObjectContainer_Dictionary, isa_proxy, 0,
	                 sizeof(::xsecAna::TPCObjectContainer) );
	instance.SetNew(&new_xsecAnacLcLTPCObjectContainer);
	instance.SetNewArray(&newArray_xsecAnacLcLTPCObjectContainer);
	instance.SetDelete(&delete_xsecAnacLcLTPCObjectContainer);
	instance.SetDeleteArray(&deleteArray_xsecAnacLcLTPCObjectContainer);
	instance.SetDestructor(&destruct_xsecAnacLcLTPCObjectContainer);
	return &instance;
}
TGenericClassInfo *GenerateInitInstance(const ::xsecAna::TPCObjectContainer*)
{
	return GenerateInitInstanceLocal((::xsecAna::TPCObjectContainer*)0);
}
// Static variable to force the class initialization
static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::xsecAna::TPCObjectContainer*)0x0); R__UseDummy(_R__UNIQUE_(Init));

// Dictionary for non-ClassDef classes
static TClass *xsecAnacLcLTPCObjectContainer_Dictionary() {
	TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::xsecAna::TPCObjectContainer*)0x0)->GetClass();
	xsecAnacLcLTPCObjectContainer_TClassManip(theClass);
	return theClass;
}

static void xsecAnacLcLTPCObjectContainer_TClassManip(TClass* ){
}

} // end of namespace ROOT

namespace ROOT {
// Wrappers around operator new
static void *new_xsecAnacLcLParticleContainer(void *p) {
	return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::xsecAna::ParticleContainer : new ::xsecAna::ParticleContainer;
}
static void *newArray_xsecAnacLcLParticleContainer(Long_t nElements, void *p) {
	return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::xsecAna::ParticleContainer[nElements] : new ::xsecAna::ParticleContainer[nElements];
}
// Wrapper around operator delete
static void delete_xsecAnacLcLParticleContainer(void *p) {
	delete ((::xsecAna::ParticleContainer*)p);
}
static void deleteArray_xsecAnacLcLParticleContainer(void *p) {
	delete [] ((::xsecAna::ParticleContainer*)p);
}
static void destruct_xsecAnacLcLParticleContainer(void *p) {
	typedef ::xsecAna::ParticleContainer current_t;
	((current_t*)p)->~current_t();
}
} // end of namespace ROOT for class ::xsecAna::ParticleContainer

namespace ROOT {
// Wrappers around operator new
static void *new_xsecAnacLcLTPCObjectContainer(void *p) {
	return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::xsecAna::TPCObjectContainer : new ::xsecAna::TPCObjectContainer;
}
static void *newArray_xsecAnacLcLTPCObjectContainer(Long_t nElements, void *p) {
	return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::xsecAna::TPCObjectContainer[nElements] : new ::xsecAna::TPCObjectContainer[nElements];
}
// Wrapper around operator delete
static void delete_xsecAnacLcLTPCObjectContainer(void *p) {
	delete ((::xsecAna::TPCObjectContainer*)p);
}
static void deleteArray_xsecAnacLcLTPCObjectContainer(void *p) {
	delete [] ((::xsecAna::TPCObjectContainer*)p);
}
static void destruct_xsecAnacLcLTPCObjectContainer(void *p) {
	typedef ::xsecAna::TPCObjectContainer current_t;
	((current_t*)p)->~current_t();
}
} // end of namespace ROOT for class ::xsecAna::TPCObjectContainer

namespace ROOT {
static TClass *vectorlExsecAnacLcLTPCObjectContainergR_Dictionary();
static void vectorlExsecAnacLcLTPCObjectContainergR_TClassManip(TClass*);
static void *new_vectorlExsecAnacLcLTPCObjectContainergR(void *p = 0);
static void *newArray_vectorlExsecAnacLcLTPCObjectContainergR(Long_t size, void *p);
static void delete_vectorlExsecAnacLcLTPCObjectContainergR(void *p);
static void deleteArray_vectorlExsecAnacLcLTPCObjectContainergR(void *p);
static void destruct_vectorlExsecAnacLcLTPCObjectContainergR(void *p);

// Function generating the singleton type initializer
static TGenericClassInfo *GenerateInitInstanceLocal(const vector<xsecAna::TPCObjectContainer>*)
{
	vector<xsecAna::TPCObjectContainer> *ptr = 0;
	static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<xsecAna::TPCObjectContainer>));
	static ::ROOT::TGenericClassInfo
	        instance("vector<xsecAna::TPCObjectContainer>", -2, "vector", 457,
	                 typeid(vector<xsecAna::TPCObjectContainer>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
	                 &vectorlExsecAnacLcLTPCObjectContainergR_Dictionary, isa_proxy, 4,
	                 sizeof(vector<xsecAna::TPCObjectContainer>) );
	instance.SetNew(&new_vectorlExsecAnacLcLTPCObjectContainergR);
	instance.SetNewArray(&newArray_vectorlExsecAnacLcLTPCObjectContainergR);
	instance.SetDelete(&delete_vectorlExsecAnacLcLTPCObjectContainergR);
	instance.SetDeleteArray(&deleteArray_vectorlExsecAnacLcLTPCObjectContainergR);
	instance.SetDestructor(&destruct_vectorlExsecAnacLcLTPCObjectContainergR);
	instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<xsecAna::TPCObjectContainer> >()));
	return &instance;
}
// Static variable to force the class initialization
static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<xsecAna::TPCObjectContainer>*) 0x0); R__UseDummy(_R__UNIQUE_(Init));

// Dictionary for non-ClassDef classes
static TClass *vectorlExsecAnacLcLTPCObjectContainergR_Dictionary() {
	TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<xsecAna::TPCObjectContainer>*) 0x0)->GetClass();
	vectorlExsecAnacLcLTPCObjectContainergR_TClassManip(theClass);
	return theClass;
}

static void vectorlExsecAnacLcLTPCObjectContainergR_TClassManip(TClass* ){
}

} // end of namespace ROOT

namespace ROOT {
// Wrappers around operator new
static void *new_vectorlExsecAnacLcLTPCObjectContainergR(void *p) {
	return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p)vector<xsecAna::TPCObjectContainer> : new vector<xsecAna::TPCObjectContainer>;
}
static void *newArray_vectorlExsecAnacLcLTPCObjectContainergR(Long_t nElements, void *p) {
	return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p)vector<xsecAna::TPCObjectContainer>[nElements] : new vector<xsecAna::TPCObjectContainer>[nElements];
}
// Wrapper around operator delete
static void delete_vectorlExsecAnacLcLTPCObjectContainergR(void *p) {
	delete ((vector<xsecAna::TPCObjectContainer>*)p);
}
static void deleteArray_vectorlExsecAnacLcLTPCObjectContainergR(void *p) {
	delete [] ((vector<xsecAna::TPCObjectContainer>*)p);
}
static void destruct_vectorlExsecAnacLcLTPCObjectContainergR(void *p) {
	typedef vector<xsecAna::TPCObjectContainer> current_t;
	((current_t*)p)->~current_t();
}
} // end of namespace ROOT for class vector<xsecAna::TPCObjectContainer>

namespace ROOT {
static TClass *vectorlExsecAnacLcLParticleContainergR_Dictionary();
static void vectorlExsecAnacLcLParticleContainergR_TClassManip(TClass*);
static void *new_vectorlExsecAnacLcLParticleContainergR(void *p = 0);
static void *newArray_vectorlExsecAnacLcLParticleContainergR(Long_t size, void *p);
static void delete_vectorlExsecAnacLcLParticleContainergR(void *p);
static void deleteArray_vectorlExsecAnacLcLParticleContainergR(void *p);
static void destruct_vectorlExsecAnacLcLParticleContainergR(void *p);

// Function generating the singleton type initializer
static TGenericClassInfo *GenerateInitInstanceLocal(const vector<xsecAna::ParticleContainer>*)
{
	vector<xsecAna::ParticleContainer> *ptr = 0;
	static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<xsecAna::ParticleContainer>));
	static ::ROOT::TGenericClassInfo
	        instance("vector<xsecAna::ParticleContainer>", -2, "vector", 457,
	                 typeid(vector<xsecAna::ParticleContainer>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
	                 &vectorlExsecAnacLcLParticleContainergR_Dictionary, isa_proxy, 4,
	                 sizeof(vector<xsecAna::ParticleContainer>) );
	instance.SetNew(&new_vectorlExsecAnacLcLParticleContainergR);
	instance.SetNewArray(&newArray_vectorlExsecAnacLcLParticleContainergR);
	instance.SetDelete(&delete_vectorlExsecAnacLcLParticleContainergR);
	instance.SetDeleteArray(&deleteArray_vectorlExsecAnacLcLParticleContainergR);
	instance.SetDestructor(&destruct_vectorlExsecAnacLcLParticleContainergR);
	instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<xsecAna::ParticleContainer> >()));
	return &instance;
}
// Static variable to force the class initialization
static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<xsecAna::ParticleContainer>*) 0x0); R__UseDummy(_R__UNIQUE_(Init));

// Dictionary for non-ClassDef classes
static TClass *vectorlExsecAnacLcLParticleContainergR_Dictionary() {
	TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<xsecAna::ParticleContainer>*) 0x0)->GetClass();
	vectorlExsecAnacLcLParticleContainergR_TClassManip(theClass);
	return theClass;
}

static void vectorlExsecAnacLcLParticleContainergR_TClassManip(TClass* ){
}

} // end of namespace ROOT

namespace ROOT {
// Wrappers around operator new
static void *new_vectorlExsecAnacLcLParticleContainergR(void *p) {
	return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p)vector<xsecAna::ParticleContainer> : new vector<xsecAna::ParticleContainer>;
}
static void *newArray_vectorlExsecAnacLcLParticleContainergR(Long_t nElements, void *p) {
	return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p)vector<xsecAna::ParticleContainer>[nElements] : new vector<xsecAna::ParticleContainer>[nElements];
}
// Wrapper around operator delete
static void delete_vectorlExsecAnacLcLParticleContainergR(void *p) {
	delete ((vector<xsecAna::ParticleContainer>*)p);
}
static void deleteArray_vectorlExsecAnacLcLParticleContainergR(void *p) {
	delete [] ((vector<xsecAna::ParticleContainer>*)p);
}
static void destruct_vectorlExsecAnacLcLParticleContainergR(void *p) {
	typedef vector<xsecAna::ParticleContainer> current_t;
	((current_t*)p)->~current_t();
}
} // end of namespace ROOT for class vector<xsecAna::ParticleContainer>

namespace {
void TriggerDictionaryInitialization_xsecAna_dict_Impl() {
	static const char* headers[] = {
		"TpcObjectContainer.h",
		"ParticleContainer.h",
		0
	};
	static const char* includePaths[] = {
		"/usr/local/Cellar/root6/6.08.02/include/root",
		"/Users/chill/PhD/nue_analysis_modules/xsecAna/",
		0
	};
	static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "xsecAna_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace xsecAna{class __attribute__((annotate("$clingAutoload$ParticleContainer.h")))  __attribute__((annotate("$clingAutoload$TpcObjectContainer.h")))  ParticleContainer;}
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$memory")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
namespace xsecAna{class __attribute__((annotate("$clingAutoload$TpcObjectContainer.h")))  TPCObjectContainer;}
)DICTFWDDCLS";
	static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "xsecAna_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TpcObjectContainer.h"
#include "ParticleContainer.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
	static const char* classesHeaders[]={
		"xsecAna::ParticleContainer", payloadCode, "@",
		"xsecAna::TPCObjectContainer", payloadCode, "@",
		nullptr
	};

	static bool isInitialized = false;
	if (!isInitialized) {
		TROOT::RegisterModule("xsecAna_dict",
		                      headers, includePaths, payloadCode, fwdDeclCode,
		                      TriggerDictionaryInitialization_xsecAna_dict_Impl, {}, classesHeaders);
		isInitialized = true;
	}
}
static struct DictInit {
	DictInit() {
		TriggerDictionaryInitialization_xsecAna_dict_Impl();
	}
} __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_xsecAna_dict() {
	TriggerDictionaryInitialization_xsecAna_dict_Impl();
}
