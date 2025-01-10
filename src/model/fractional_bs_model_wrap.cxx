#ifndef SWIGPYTHON
#define SWIGPYTHON
#endif

#define SWIG_PYTHON_DIRECTOR_NO_VTABLE


#ifdef __cplusplus
/* SwigValueWrapper is described in swig.swg */
template<typename T> class SwigValueWrapper {
  struct SwigMovePointer {
    T *ptr;
    SwigMovePointer(T *p) : ptr(p) { }
    ~SwigMovePointer() { delete ptr; }
    SwigMovePointer& operator=(SwigMovePointer& rhs) { T* oldptr = ptr; ptr = 0; delete oldptr; ptr = rhs.ptr; rhs.ptr = 0; return *this; }
  } pointer;
  SwigValueWrapper& operator=(const SwigValueWrapper<T>& rhs);
  SwigValueWrapper(const SwigValueWrapper<T>& rhs);
public:
  SwigValueWrapper() : pointer(0) { }
  SwigValueWrapper& operator=(const T& t) { SwigMovePointer tmp(new T(t)); pointer = tmp; return *this; }
  operator T&() const { return *pointer.ptr; }
  T *operator&() { return pointer.ptr; }
};

template <typename T> T SwigValueInit() {
  return T();
}
#endif

/* -----------------------------------------------------------------------------
 *  This section contains generic SWIG labels for method/variable
 *  declarations/attributes, and other compiler dependent labels.
 * ----------------------------------------------------------------------------- */

/* template workaround for compilers that cannot correctly implement the C++ standard */
#ifndef SWIGTEMPLATEDISAMBIGUATOR
# if defined(__SUNPRO_CC) && (__SUNPRO_CC <= 0x560)
#  define SWIGTEMPLATEDISAMBIGUATOR template
# elif defined(__HP_aCC)
/* Needed even with `aCC -AA' when `aCC -V' reports HP ANSI C++ B3910B A.03.55 */
/* If we find a maximum version that requires this, the test would be __HP_aCC <= 35500 for A.03.55 */
#  define SWIGTEMPLATEDISAMBIGUATOR template
# else
#  define SWIGTEMPLATEDISAMBIGUATOR
# endif
#endif

/* inline attribute */
#ifndef SWIGINLINE
# if defined(__cplusplus) || (defined(__GNUC__) && !defined(__STRICT_ANSI__))
#   define SWIGINLINE inline
# else
#   define SWIGINLINE
# endif
#endif

/* attribute recognised by some compilers to avoid 'unused' warnings */
#ifndef SWIGUNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define SWIGUNUSED __attribute__ ((__unused__))
#   else
#     define SWIGUNUSED
#   endif
# elif defined(__ICC)
#   define SWIGUNUSED __attribute__ ((__unused__))
# else
#   define SWIGUNUSED
# endif
#endif

#ifndef SWIG_MSC_UNSUPPRESS_4505
# if defined(_MSC_VER)
#   pragma warning(disable : 4505) /* unreferenced local function has been removed */
# endif
#endif

#ifndef SWIGUNUSEDPARM
# ifdef __cplusplus
#   define SWIGUNUSEDPARM(p)
# else
#   define SWIGUNUSEDPARM(p) p SWIGUNUSED
# endif
#endif

/* internal SWIG method */
#ifndef SWIGINTERN
# define SWIGINTERN static SWIGUNUSED
#endif

/* internal inline SWIG method */
#ifndef SWIGINTERNINLINE
# define SWIGINTERNINLINE SWIGINTERN SWIGINLINE
#endif

/* exporting methods */
#if defined(__GNUC__)
#  if (__GNUC__ >= 4) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#    ifndef GCC_HASCLASSVISIBILITY
#      define GCC_HASCLASSVISIBILITY
#    endif
#  endif
#endif

#ifndef SWIGEXPORT
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   if defined(STATIC_LINKED)
#     define SWIGEXPORT
#   else
#     define SWIGEXPORT __declspec(dllexport)
#   endif
# else
#   if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#     define SWIGEXPORT __attribute__ ((visibility("default")))
#   else
#     define SWIGEXPORT
#   endif
# endif
#endif

/* calling conventions for Windows */
#ifndef SWIGSTDCALL
# if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#   define SWIGSTDCALL __stdcall
# else
#   define SWIGSTDCALL
# endif
#endif

/* Deal with Microsoft's attempt at deprecating C standard runtime functions */
#if !defined(SWIG_NO_CRT_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_CRT_SECURE_NO_DEPRECATE)
# define _CRT_SECURE_NO_DEPRECATE
#endif

/* Deal with Microsoft's attempt at deprecating methods in the standard C++ library */
#if !defined(SWIG_NO_SCL_SECURE_NO_DEPRECATE) && defined(_MSC_VER) && !defined(_SCL_SECURE_NO_DEPRECATE)
# define _SCL_SECURE_NO_DEPRECATE
#endif

/* Deal with Apple's deprecated 'AssertMacros.h' from Carbon-framework */
#if defined(__APPLE__) && !defined(__ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES)
# define __ASSERT_MACROS_DEFINE_VERSIONS_WITHOUT_UNDERSCORES 0
#endif

/* Intel's compiler complains if a variable which was never initialised is
 * cast to void, which is a common idiom which we use to indicate that we
 * are aware a variable isn't used.  So we just silence that warning.
 * See: https://github.com/swig/swig/issues/192 for more discussion.
 */
#ifdef __INTEL_COMPILER
# pragma warning disable 592
#endif


#if defined(__GNUC__) && defined(_WIN32) && !defined(SWIG_PYTHON_NO_HYPOT_WORKAROUND)
/* Workaround for '::hypot' has not been declared', see https://bugs.python.org/issue11566 */
# include <math.h>
#endif

#if defined(_DEBUG) && defined(SWIG_PYTHON_INTERPRETER_NO_DEBUG)
/* Use debug wrappers with the Python release dll */
# undef _DEBUG
# include <Python.h>
# define _DEBUG 1
#else
# include <Python.h>
#endif

/* -----------------------------------------------------------------------------
 * swigrun.swg
 *
 * This file contains generic C API SWIG runtime support for pointer
 * type checking.
 * ----------------------------------------------------------------------------- */

/* This should only be incremented when either the layout of swig_type_info changes,
   or for whatever reason, the runtime changes incompatibly */
#define SWIG_RUNTIME_VERSION "4"

/* define SWIG_TYPE_TABLE_NAME as "SWIG_TYPE_TABLE" */
#ifdef SWIG_TYPE_TABLE
# define SWIG_QUOTE_STRING(x) #x
# define SWIG_EXPAND_AND_QUOTE_STRING(x) SWIG_QUOTE_STRING(x)
# define SWIG_TYPE_TABLE_NAME SWIG_EXPAND_AND_QUOTE_STRING(SWIG_TYPE_TABLE)
#else
# define SWIG_TYPE_TABLE_NAME
#endif

/*
  You can use the SWIGRUNTIME and SWIGRUNTIMEINLINE macros for
  creating a static or dynamic library from the SWIG runtime code.
  In 99.9% of the cases, SWIG just needs to declare them as 'static'.

  But only do this if strictly necessary, ie, if you have problems
  with your compiler or suchlike.
*/

#ifndef SWIGRUNTIME
# define SWIGRUNTIME SWIGINTERN
#endif

#ifndef SWIGRUNTIMEINLINE
# define SWIGRUNTIMEINLINE SWIGRUNTIME SWIGINLINE
#endif

/*  Generic buffer size */
#ifndef SWIG_BUFFER_SIZE
# define SWIG_BUFFER_SIZE 1024
#endif

/* Flags for pointer conversions */
#define SWIG_POINTER_DISOWN        0x1
#define SWIG_CAST_NEW_MEMORY       0x2
#define SWIG_POINTER_NO_NULL       0x4

/* Flags for new pointer objects */
#define SWIG_POINTER_OWN           0x1


/*
   Flags/methods for returning states.

   The SWIG conversion methods, as ConvertPtr, return an integer
   that tells if the conversion was successful or not. And if not,
   an error code can be returned (see swigerrors.swg for the codes).

   Use the following macros/flags to set or process the returning
   states.

   In old versions of SWIG, code such as the following was usually written:

     if (SWIG_ConvertPtr(obj,vptr,ty.flags) != -1) {
       // success code
     } else {
       //fail code
     }

   Now you can be more explicit:

    int res = SWIG_ConvertPtr(obj,vptr,ty.flags);
    if (SWIG_IsOK(res)) {
      // success code
    } else {
      // fail code
    }

   which is the same really, but now you can also do

    Type *ptr;
    int res = SWIG_ConvertPtr(obj,(void **)(&ptr),ty.flags);
    if (SWIG_IsOK(res)) {
      // success code
      if (SWIG_IsNewObj(res) {
        ...
	delete *ptr;
      } else {
        ...
      }
    } else {
      // fail code
    }

   I.e., now SWIG_ConvertPtr can return new objects and you can
   identify the case and take care of the deallocation. Of course that
   also requires SWIG_ConvertPtr to return new result values, such as

      int SWIG_ConvertPtr(obj, ptr,...) {
        if (<obj is ok>) {
          if (<need new object>) {
            *ptr = <ptr to new allocated object>;
            return SWIG_NEWOBJ;
          } else {
            *ptr = <ptr to old object>;
            return SWIG_OLDOBJ;
          }
        } else {
          return SWIG_BADOBJ;
        }
      }

   Of course, returning the plain '0(success)/-1(fail)' still works, but you can be
   more explicit by returning SWIG_BADOBJ, SWIG_ERROR or any of the
   SWIG errors code.

   Finally, if the SWIG_CASTRANK_MODE is enabled, the result code
   allows to return the 'cast rank', for example, if you have this

       int food(double)
       int fooi(int);

   and you call

      food(1)   // cast rank '1'  (1 -> 1.0)
      fooi(1)   // cast rank '0'

   just use the SWIG_AddCast()/SWIG_CheckState()
*/

#define SWIG_OK                    (0)
#define SWIG_ERROR                 (-1)
#define SWIG_IsOK(r)               (r >= 0)
#define SWIG_ArgError(r)           ((r != SWIG_ERROR) ? r : SWIG_TypeError)

/* The CastRankLimit says how many bits are used for the cast rank */
#define SWIG_CASTRANKLIMIT         (1 << 8)
/* The NewMask denotes the object was created (using new/malloc) */
#define SWIG_NEWOBJMASK            (SWIG_CASTRANKLIMIT  << 1)
/* The TmpMask is for in/out typemaps that use temporal objects */
#define SWIG_TMPOBJMASK            (SWIG_NEWOBJMASK << 1)
/* Simple returning values */
#define SWIG_BADOBJ                (SWIG_ERROR)
#define SWIG_OLDOBJ                (SWIG_OK)
#define SWIG_NEWOBJ                (SWIG_OK | SWIG_NEWOBJMASK)
#define SWIG_TMPOBJ                (SWIG_OK | SWIG_TMPOBJMASK)
/* Check, add and del mask methods */
#define SWIG_AddNewMask(r)         (SWIG_IsOK(r) ? (r | SWIG_NEWOBJMASK) : r)
#define SWIG_DelNewMask(r)         (SWIG_IsOK(r) ? (r & ~SWIG_NEWOBJMASK) : r)
#define SWIG_IsNewObj(r)           (SWIG_IsOK(r) && (r & SWIG_NEWOBJMASK))
#define SWIG_AddTmpMask(r)         (SWIG_IsOK(r) ? (r | SWIG_TMPOBJMASK) : r)
#define SWIG_DelTmpMask(r)         (SWIG_IsOK(r) ? (r & ~SWIG_TMPOBJMASK) : r)
#define SWIG_IsTmpObj(r)           (SWIG_IsOK(r) && (r & SWIG_TMPOBJMASK))

/* Cast-Rank Mode */
#if defined(SWIG_CASTRANK_MODE)
#  ifndef SWIG_TypeRank
#    define SWIG_TypeRank             unsigned long
#  endif
#  ifndef SWIG_MAXCASTRANK            /* Default cast allowed */
#    define SWIG_MAXCASTRANK          (2)
#  endif
#  define SWIG_CASTRANKMASK          ((SWIG_CASTRANKLIMIT) -1)
#  define SWIG_CastRank(r)           (r & SWIG_CASTRANKMASK)
SWIGINTERNINLINE int SWIG_AddCast(int r) {
  return SWIG_IsOK(r) ? ((SWIG_CastRank(r) < SWIG_MAXCASTRANK) ? (r + 1) : SWIG_ERROR) : r;
}
SWIGINTERNINLINE int SWIG_CheckState(int r) {
  return SWIG_IsOK(r) ? SWIG_CastRank(r) + 1 : 0;
}
#else /* no cast-rank mode */
#  define SWIG_AddCast(r) (r)
#  define SWIG_CheckState(r) (SWIG_IsOK(r) ? 1 : 0)
#endif


#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void *(*swig_converter_func)(void *, int *);
typedef struct swig_type_info *(*swig_dycast_func)(void **);

/* Structure to store information on one type */
typedef struct swig_type_info {
  const char             *name;			/* mangled name of this type */
  const char             *str;			/* human readable name of this type */
  swig_dycast_func        dcast;		/* dynamic cast function down a hierarchy */
  struct swig_cast_info  *cast;			/* linked list of types that can cast into this type */
  void                   *clientdata;		/* language specific type data */
  int                    owndata;		/* flag if the structure owns the clientdata */
} swig_type_info;

/* Structure to store a type and conversion function used for casting */
typedef struct swig_cast_info {
  swig_type_info         *type;			/* pointer to type that is equivalent to this type */
  swig_converter_func     converter;		/* function to cast the void pointers */
  struct swig_cast_info  *next;			/* pointer to next cast in linked list */
  struct swig_cast_info  *prev;			/* pointer to the previous cast */
} swig_cast_info;

/* Structure used to store module information
 * Each module generates one structure like this, and the runtime collects
 * all of these structures and stores them in a circularly linked list.*/
typedef struct swig_module_info {
  swig_type_info         **types;		/* Array of pointers to swig_type_info structures that are in this module */
  size_t                 size;		        /* Number of types in this module */
  struct swig_module_info *next;		/* Pointer to next element in circularly linked list */
  swig_type_info         **type_initial;	/* Array of initially generated type structures */
  swig_cast_info         **cast_initial;	/* Array of initially generated casting structures */
  void                    *clientdata;		/* Language specific module data */
} swig_module_info;

/*
  Compare two type names skipping the space characters, therefore
  "char*" == "char *" and "Class<int>" == "Class<int >", etc.

  Return 0 when the two name types are equivalent, as in
  strncmp, but skipping ' '.
*/
SWIGRUNTIME int
SWIG_TypeNameComp(const char *f1, const char *l1,
		  const char *f2, const char *l2) {
  for (;(f1 != l1) && (f2 != l2); ++f1, ++f2) {
    while ((*f1 == ' ') && (f1 != l1)) ++f1;
    while ((*f2 == ' ') && (f2 != l2)) ++f2;
    if (*f1 != *f2) return (*f1 > *f2) ? 1 : -1;
  }
  return (int)((l1 - f1) - (l2 - f2));
}

/*
  Check type equivalence in a name list like <name1>|<name2>|...
  Return 0 if equal, -1 if nb < tb, 1 if nb > tb
*/
SWIGRUNTIME int
SWIG_TypeCmp(const char *nb, const char *tb) {
  int equiv = 1;
  const char* te = tb + strlen(tb);
  const char* ne = nb;
  while (equiv != 0 && *ne) {
    for (nb = ne; *ne; ++ne) {
      if (*ne == '|') break;
    }
    equiv = SWIG_TypeNameComp(nb, ne, tb, te);
    if (*ne) ++ne;
  }
  return equiv;
}

/*
  Check type equivalence in a name list like <name1>|<name2>|...
  Return 0 if not equal, 1 if equal
*/
SWIGRUNTIME int
SWIG_TypeEquiv(const char *nb, const char *tb) {
  return SWIG_TypeCmp(nb, tb) == 0 ? 1 : 0;
}

/*
  Check the typename
*/
SWIGRUNTIME swig_cast_info *
SWIG_TypeCheck(const char *c, swig_type_info *ty) {
  if (ty) {
    swig_cast_info *iter = ty->cast;
    while (iter) {
      if (strcmp(iter->type->name, c) == 0) {
        if (iter == ty->cast)
          return iter;
        /* Move iter to the top of the linked list */
        iter->prev->next = iter->next;
        if (iter->next)
          iter->next->prev = iter->prev;
        iter->next = ty->cast;
        iter->prev = 0;
        if (ty->cast) ty->cast->prev = iter;
        ty->cast = iter;
        return iter;
      }
      iter = iter->next;
    }
  }
  return 0;
}

/*
  Identical to SWIG_TypeCheck, except strcmp is replaced with a pointer comparison
*/
SWIGRUNTIME swig_cast_info *
SWIG_TypeCheckStruct(swig_type_info *from, swig_type_info *ty) {
  if (ty) {
    swig_cast_info *iter = ty->cast;
    while (iter) {
      if (iter->type == from) {
        if (iter == ty->cast)
          return iter;
        /* Move iter to the top of the linked list */
        iter->prev->next = iter->next;
        if (iter->next)
          iter->next->prev = iter->prev;
        iter->next = ty->cast;
        iter->prev = 0;
        if (ty->cast) ty->cast->prev = iter;
        ty->cast = iter;
        return iter;
      }
      iter = iter->next;
    }
  }
  return 0;
}

/*
  Cast a pointer up an inheritance hierarchy
*/
SWIGRUNTIMEINLINE void *
SWIG_TypeCast(swig_cast_info *ty, void *ptr, int *newmemory) {
  return ((!ty) || (!ty->converter)) ? ptr : (*ty->converter)(ptr, newmemory);
}

/*
   Dynamic pointer casting. Down an inheritance hierarchy
*/
SWIGRUNTIME swig_type_info *
SWIG_TypeDynamicCast(swig_type_info *ty, void **ptr) {
  swig_type_info *lastty = ty;
  if (!ty || !ty->dcast) return ty;
  while (ty && (ty->dcast)) {
    ty = (*ty->dcast)(ptr);
    if (ty) lastty = ty;
  }
  return lastty;
}

/*
  Return the name associated with this type
*/
SWIGRUNTIMEINLINE const char *
SWIG_TypeName(const swig_type_info *ty) {
  return ty->name;
}

/*
  Return the pretty name associated with this type,
  that is an unmangled type name in a form presentable to the user.
*/
SWIGRUNTIME const char *
SWIG_TypePrettyName(const swig_type_info *type) {
  /* The "str" field contains the equivalent pretty names of the
     type, separated by vertical-bar characters.  We choose
     to print the last name, as it is often (?) the most
     specific. */
  if (!type) return NULL;
  if (type->str != NULL) {
    const char *last_name = type->str;
    const char *s;
    for (s = type->str; *s; s++)
      if (*s == '|') last_name = s+1;
    return last_name;
  }
  else
    return type->name;
}

/*
   Set the clientdata field for a type
*/
SWIGRUNTIME void
SWIG_TypeClientData(swig_type_info *ti, void *clientdata) {
  swig_cast_info *cast = ti->cast;
  /* if (ti->clientdata == clientdata) return; */
  ti->clientdata = clientdata;

  while (cast) {
    if (!cast->converter) {
      swig_type_info *tc = cast->type;
      if (!tc->clientdata) {
	SWIG_TypeClientData(tc, clientdata);
      }
    }
    cast = cast->next;
  }
}
SWIGRUNTIME void
SWIG_TypeNewClientData(swig_type_info *ti, void *clientdata) {
  SWIG_TypeClientData(ti, clientdata);
  ti->owndata = 1;
}

/*
  Search for a swig_type_info structure only by mangled name
  Search is a O(log #types)

  We start searching at module start, and finish searching when start == end.
  Note: if start == end at the beginning of the function, we go all the way around
  the circular list.
*/
SWIGRUNTIME swig_type_info *
SWIG_MangledTypeQueryModule(swig_module_info *start,
                            swig_module_info *end,
		            const char *name) {
  swig_module_info *iter = start;
  do {
    if (iter->size) {
      size_t l = 0;
      size_t r = iter->size - 1;
      do {
	/* since l+r >= 0, we can (>> 1) instead (/ 2) */
	size_t i = (l + r) >> 1;
	const char *iname = iter->types[i]->name;
	if (iname) {
	  int compare = strcmp(name, iname);
	  if (compare == 0) {
	    return iter->types[i];
	  } else if (compare < 0) {
	    if (i) {
	      r = i - 1;
	    } else {
	      break;
	    }
	  } else if (compare > 0) {
	    l = i + 1;
	  }
	} else {
	  break; /* should never happen */
	}
      } while (l <= r);
    }
    iter = iter->next;
  } while (iter != end);
  return 0;
}

/*
  Search for a swig_type_info structure for either a mangled name or a human readable name.
  It first searches the mangled names of the types, which is a O(log #types)
  If a type is not found it then searches the human readable names, which is O(#types).

  We start searching at module start, and finish searching when start == end.
  Note: if start == end at the beginning of the function, we go all the way around
  the circular list.
*/
SWIGRUNTIME swig_type_info *
SWIG_TypeQueryModule(swig_module_info *start,
                     swig_module_info *end,
		     const char *name) {
  /* STEP 1: Search the name field using binary search */
  swig_type_info *ret = SWIG_MangledTypeQueryModule(start, end, name);
  if (ret) {
    return ret;
  } else {
    /* STEP 2: If the type hasn't been found, do a complete search
       of the str field (the human readable name) */
    swig_module_info *iter = start;
    do {
      size_t i = 0;
      for (; i < iter->size; ++i) {
	if (iter->types[i]->str && (SWIG_TypeEquiv(iter->types[i]->str, name)))
	  return iter->types[i];
      }
      iter = iter->next;
    } while (iter != end);
  }

  /* neither found a match */
  return 0;
}

/*
   Pack binary data into a string
*/
SWIGRUNTIME char *
SWIG_PackData(char *c, void *ptr, size_t sz) {
  static const char hex[17] = "0123456789abcdef";
  const unsigned char *u = (unsigned char *) ptr;
  const unsigned char *eu =  u + sz;
  for (; u != eu; ++u) {
    unsigned char uu = *u;
    *(c++) = hex[(uu & 0xf0) >> 4];
    *(c++) = hex[uu & 0xf];
  }
  return c;
}

/*
   Unpack binary data from a string
*/
SWIGRUNTIME const char *
SWIG_UnpackData(const char *c, void *ptr, size_t sz) {
  unsigned char *u = (unsigned char *) ptr;
  const unsigned char *eu = u + sz;
  for (; u != eu; ++u) {
    char d = *(c++);
    unsigned char uu;
    if ((d >= '0') && (d <= '9'))
      uu = (unsigned char)((d - '0') << 4);
    else if ((d >= 'a') && (d <= 'f'))
      uu = (unsigned char)((d - ('a'-10)) << 4);
    else
      return (char *) 0;
    d = *(c++);
    if ((d >= '0') && (d <= '9'))
      uu |= (unsigned char)(d - '0');
    else if ((d >= 'a') && (d <= 'f'))
      uu |= (unsigned char)(d - ('a'-10));
    else
      return (char *) 0;
    *u = uu;
  }
  return c;
}

/*
   Pack 'void *' into a string buffer.
*/
SWIGRUNTIME char *
SWIG_PackVoidPtr(char *buff, void *ptr, const char *name, size_t bsz) {
  char *r = buff;
  if ((2*sizeof(void *) + 2) > bsz) return 0;
  *(r++) = '_';
  r = SWIG_PackData(r,&ptr,sizeof(void *));
  if (strlen(name) + 1 > (bsz - (r - buff))) return 0;
  strcpy(r,name);
  return buff;
}

SWIGRUNTIME const char *
SWIG_UnpackVoidPtr(const char *c, void **ptr, const char *name) {
  if (*c != '_') {
    if (strcmp(c,"NULL") == 0) {
      *ptr = (void *) 0;
      return name;
    } else {
      return 0;
    }
  }
  return SWIG_UnpackData(++c,ptr,sizeof(void *));
}

SWIGRUNTIME char *
SWIG_PackDataName(char *buff, void *ptr, size_t sz, const char *name, size_t bsz) {
  char *r = buff;
  size_t lname = (name ? strlen(name) : 0);
  if ((2*sz + 2 + lname) > bsz) return 0;
  *(r++) = '_';
  r = SWIG_PackData(r,ptr,sz);
  if (lname) {
    strncpy(r,name,lname+1);
  } else {
    *r = 0;
  }
  return buff;
}

SWIGRUNTIME const char *
SWIG_UnpackDataName(const char *c, void *ptr, size_t sz, const char *name) {
  if (*c != '_') {
    if (strcmp(c,"NULL") == 0) {
      memset(ptr,0,sz);
      return name;
    } else {
      return 0;
    }
  }
  return SWIG_UnpackData(++c,ptr,sz);
}

#ifdef __cplusplus
}
#endif

/*  Errors in SWIG */
#define  SWIG_UnknownError    	   -1
#define  SWIG_IOError        	   -2
#define  SWIG_RuntimeError   	   -3
#define  SWIG_IndexError     	   -4
#define  SWIG_TypeError      	   -5
#define  SWIG_DivisionByZero 	   -6
#define  SWIG_OverflowError  	   -7
#define  SWIG_SyntaxError    	   -8
#define  SWIG_ValueError     	   -9
#define  SWIG_SystemError    	   -10
#define  SWIG_AttributeError 	   -11
#define  SWIG_MemoryError    	   -12
#define  SWIG_NullReferenceError   -13



/* Compatibility macros for Python 3 */
#if PY_VERSION_HEX >= 0x03000000

#define PyClass_Check(obj) PyObject_IsInstance(obj, (PyObject *)&PyType_Type)
#define PyInt_Check(x) PyLong_Check(x)
#define PyInt_AsLong(x) PyLong_AsLong(x)
#define PyInt_FromLong(x) PyLong_FromLong(x)
#define PyInt_FromSize_t(x) PyLong_FromSize_t(x)
#define PyString_Check(name) PyBytes_Check(name)
#define PyString_FromString(x) PyUnicode_FromString(x)
#define PyString_Format(fmt, args)  PyUnicode_Format(fmt, args)
#define PyString_AsString(str) PyBytes_AsString(str)
#define PyString_Size(str) PyBytes_Size(str)	
#define PyString_InternFromString(key) PyUnicode_InternFromString(key)
#define Py_TPFLAGS_HAVE_CLASS Py_TPFLAGS_BASETYPE
#define PyString_AS_STRING(x) PyUnicode_AS_STRING(x)
#define _PyLong_FromSsize_t(x) PyLong_FromSsize_t(x)

#endif

#ifndef Py_TYPE
#  define Py_TYPE(op) ((op)->ob_type)
#endif

/* SWIG APIs for compatibility of both Python 2 & 3 */

#if PY_VERSION_HEX >= 0x03000000
#  define SWIG_Python_str_FromFormat PyUnicode_FromFormat
#else
#  define SWIG_Python_str_FromFormat PyString_FromFormat
#endif


/* Warning: This function will allocate a new string in Python 3,
 * so please call SWIG_Python_str_DelForPy3(x) to free the space.
 */
SWIGINTERN char*
SWIG_Python_str_AsChar(PyObject *str)
{
#if PY_VERSION_HEX >= 0x03030000
  return (char *)PyUnicode_AsUTF8(str);
#elif PY_VERSION_HEX >= 0x03000000
  char *newstr = 0;
  str = PyUnicode_AsUTF8String(str);
  if (str) {
    char *cstr;
    Py_ssize_t len;
    if (PyBytes_AsStringAndSize(str, &cstr, &len) != -1) {
      newstr = (char *) malloc(len+1);
      if (newstr)
        memcpy(newstr, cstr, len+1);
    }
    Py_XDECREF(str);
  }
  return newstr;
#else
  return PyString_AsString(str);
#endif
}

#if PY_VERSION_HEX >= 0x03030000 || PY_VERSION_HEX < 0x03000000
#  define SWIG_Python_str_DelForPy3(x)
#else
#  define SWIG_Python_str_DelForPy3(x) free( (void*) (x) )
#endif


SWIGINTERN PyObject*
SWIG_Python_str_FromChar(const char *c)
{
#if PY_VERSION_HEX >= 0x03000000
  return PyUnicode_FromString(c); 
#else
  return PyString_FromString(c);
#endif
}

#ifndef PyObject_DEL
# define PyObject_DEL PyObject_Del
#endif

// SWIGPY_USE_CAPSULE is no longer used within SWIG itself, but some user
// interface files check for it.
# define SWIGPY_USE_CAPSULE
# define SWIGPY_CAPSULE_NAME ("swig_runtime_data" SWIG_RUNTIME_VERSION ".type_pointer_capsule" SWIG_TYPE_TABLE_NAME)

#if PY_VERSION_HEX < 0x03020000
#define PyDescr_TYPE(x) (((PyDescrObject *)(x))->d_type)
#define PyDescr_NAME(x) (((PyDescrObject *)(x))->d_name)
#define Py_hash_t long
#endif

/* -----------------------------------------------------------------------------
 * error manipulation
 * ----------------------------------------------------------------------------- */

SWIGRUNTIME PyObject*
SWIG_Python_ErrorType(int code) {
  PyObject* type = 0;
  switch(code) {
  case SWIG_MemoryError:
    type = PyExc_MemoryError;
    break;
  case SWIG_IOError:
    type = PyExc_IOError;
    break;
  case SWIG_RuntimeError:
    type = PyExc_RuntimeError;
    break;
  case SWIG_IndexError:
    type = PyExc_IndexError;
    break;
  case SWIG_TypeError:
    type = PyExc_TypeError;
    break;
  case SWIG_DivisionByZero:
    type = PyExc_ZeroDivisionError;
    break;
  case SWIG_OverflowError:
    type = PyExc_OverflowError;
    break;
  case SWIG_SyntaxError:
    type = PyExc_SyntaxError;
    break;
  case SWIG_ValueError:
    type = PyExc_ValueError;
    break;
  case SWIG_SystemError:
    type = PyExc_SystemError;
    break;
  case SWIG_AttributeError:
    type = PyExc_AttributeError;
    break;
  default:
    type = PyExc_RuntimeError;
  }
  return type;
}


SWIGRUNTIME void
SWIG_Python_AddErrorMsg(const char* mesg)
{
  PyObject *type = 0;
  PyObject *value = 0;
  PyObject *traceback = 0;

  if (PyErr_Occurred())
    PyErr_Fetch(&type, &value, &traceback);
  if (value) {
    PyObject *old_str = PyObject_Str(value);
    const char *tmp = SWIG_Python_str_AsChar(old_str);
    PyErr_Clear();
    Py_XINCREF(type);
    if (tmp)
      PyErr_Format(type, "%s %s", tmp, mesg);
    else
      PyErr_Format(type, "%s", mesg);
    SWIG_Python_str_DelForPy3(tmp);
    Py_DECREF(old_str);
    Py_DECREF(value);
  } else {
    PyErr_SetString(PyExc_RuntimeError, mesg);
  }
}

SWIGRUNTIME int
SWIG_Python_TypeErrorOccurred(PyObject *obj)
{
  PyObject *error;
  if (obj)
    return 0;
  error = PyErr_Occurred();
  return error && PyErr_GivenExceptionMatches(error, PyExc_TypeError);
}

SWIGRUNTIME void
SWIG_Python_RaiseOrModifyTypeError(const char *message)
{
  if (SWIG_Python_TypeErrorOccurred(NULL)) {
    /* Use existing TypeError to preserve stacktrace and enhance with given message */
    PyObject *newvalue;
    PyObject *type = NULL, *value = NULL, *traceback = NULL;
    PyErr_Fetch(&type, &value, &traceback);
#if PY_VERSION_HEX >= 0x03000000
    newvalue = PyUnicode_FromFormat("%S\nAdditional information:\n%s", value, message);
#else
    newvalue = PyString_FromFormat("%s\nAdditional information:\n%s", PyString_AsString(value), message);
#endif
    Py_XDECREF(value);
    PyErr_Restore(type, newvalue, traceback);
  } else {
    /* Raise TypeError using given message */
    PyErr_SetString(PyExc_TypeError, message);
  }
}

#if defined(SWIG_PYTHON_NO_THREADS)
#  if defined(SWIG_PYTHON_THREADS)
#    undef SWIG_PYTHON_THREADS
#  endif
#endif
#if defined(SWIG_PYTHON_THREADS) /* Threading support is enabled */
#  if !defined(SWIG_PYTHON_USE_GIL) && !defined(SWIG_PYTHON_NO_USE_GIL)
#    define SWIG_PYTHON_USE_GIL
#  endif
#  if defined(SWIG_PYTHON_USE_GIL) /* Use PyGILState threads calls */
#    ifndef SWIG_PYTHON_INITIALIZE_THREADS
#     define SWIG_PYTHON_INITIALIZE_THREADS  PyEval_InitThreads() 
#    endif
#    ifdef __cplusplus /* C++ code */
       class SWIG_Python_Thread_Block {
         bool status;
         PyGILState_STATE state;
       public:
         void end() { if (status) { PyGILState_Release(state); status = false;} }
         SWIG_Python_Thread_Block() : status(true), state(PyGILState_Ensure()) {}
         ~SWIG_Python_Thread_Block() { end(); }
       };
       class SWIG_Python_Thread_Allow {
         bool status;
         PyThreadState *save;
       public:
         void end() { if (status) { PyEval_RestoreThread(save); status = false; }}
         SWIG_Python_Thread_Allow() : status(true), save(PyEval_SaveThread()) {}
         ~SWIG_Python_Thread_Allow() { end(); }
       };
#      define SWIG_PYTHON_THREAD_BEGIN_BLOCK   SWIG_Python_Thread_Block _swig_thread_block
#      define SWIG_PYTHON_THREAD_END_BLOCK     _swig_thread_block.end()
#      define SWIG_PYTHON_THREAD_BEGIN_ALLOW   SWIG_Python_Thread_Allow _swig_thread_allow
#      define SWIG_PYTHON_THREAD_END_ALLOW     _swig_thread_allow.end()
#    else /* C code */
#      define SWIG_PYTHON_THREAD_BEGIN_BLOCK   PyGILState_STATE _swig_thread_block = PyGILState_Ensure()
#      define SWIG_PYTHON_THREAD_END_BLOCK     PyGILState_Release(_swig_thread_block)
#      define SWIG_PYTHON_THREAD_BEGIN_ALLOW   PyThreadState *_swig_thread_allow = PyEval_SaveThread()
#      define SWIG_PYTHON_THREAD_END_ALLOW     PyEval_RestoreThread(_swig_thread_allow)
#    endif
#  else /* Old thread way, not implemented, user must provide it */
#    if !defined(SWIG_PYTHON_INITIALIZE_THREADS)
#      define SWIG_PYTHON_INITIALIZE_THREADS
#    endif
#    if !defined(SWIG_PYTHON_THREAD_BEGIN_BLOCK)
#      define SWIG_PYTHON_THREAD_BEGIN_BLOCK
#    endif
#    if !defined(SWIG_PYTHON_THREAD_END_BLOCK)
#      define SWIG_PYTHON_THREAD_END_BLOCK
#    endif
#    if !defined(SWIG_PYTHON_THREAD_BEGIN_ALLOW)
#      define SWIG_PYTHON_THREAD_BEGIN_ALLOW
#    endif
#    if !defined(SWIG_PYTHON_THREAD_END_ALLOW)
#      define SWIG_PYTHON_THREAD_END_ALLOW
#    endif
#  endif
#else /* No thread support */
#  define SWIG_PYTHON_INITIALIZE_THREADS
#  define SWIG_PYTHON_THREAD_BEGIN_BLOCK
#  define SWIG_PYTHON_THREAD_END_BLOCK
#  define SWIG_PYTHON_THREAD_BEGIN_ALLOW
#  define SWIG_PYTHON_THREAD_END_ALLOW
#endif

/* -----------------------------------------------------------------------------
 * Python API portion that goes into the runtime
 * ----------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------------------
 * Constant declarations
 * ----------------------------------------------------------------------------- */

/* Constant Types */
#define SWIG_PY_POINTER 4
#define SWIG_PY_BINARY  5

/* Constant information structure */
typedef struct swig_const_info {
  int type;
  const char *name;
  long lvalue;
  double dvalue;
  void   *pvalue;
  swig_type_info **ptype;
} swig_const_info;

#ifdef __cplusplus
}
#endif


/* -----------------------------------------------------------------------------
 * pyrun.swg
 *
 * This file contains the runtime support for Python modules
 * and includes code for managing global variables and pointer
 * type checking.
 *
 * ----------------------------------------------------------------------------- */

#if PY_VERSION_HEX < 0x02070000 /* 2.7.0 */
# error "This version of SWIG only supports Python >= 2.7"
#endif

#if PY_VERSION_HEX >= 0x03000000 && PY_VERSION_HEX < 0x03020000
# error "This version of SWIG only supports Python 3 >= 3.2"
#endif

/* Common SWIG API */

/* for raw pointers */
#define SWIG_Python_ConvertPtr(obj, pptr, type, flags)  SWIG_Python_ConvertPtrAndOwn(obj, pptr, type, flags, 0)
#define SWIG_ConvertPtr(obj, pptr, type, flags)         SWIG_Python_ConvertPtr(obj, pptr, type, flags)
#define SWIG_ConvertPtrAndOwn(obj,pptr,type,flags,own)  SWIG_Python_ConvertPtrAndOwn(obj, pptr, type, flags, own)

#ifdef SWIGPYTHON_BUILTIN
#define SWIG_NewPointerObj(ptr, type, flags)            SWIG_Python_NewPointerObj(self, ptr, type, flags)
#else
#define SWIG_NewPointerObj(ptr, type, flags)            SWIG_Python_NewPointerObj(NULL, ptr, type, flags)
#endif

#define SWIG_InternalNewPointerObj(ptr, type, flags)	SWIG_Python_NewPointerObj(NULL, ptr, type, flags)

#define SWIG_CheckImplicit(ty)                          SWIG_Python_CheckImplicit(ty) 
#define SWIG_AcquirePtr(ptr, src)                       SWIG_Python_AcquirePtr(ptr, src)
#define swig_owntype                                    int

/* for raw packed data */
#define SWIG_ConvertPacked(obj, ptr, sz, ty)            SWIG_Python_ConvertPacked(obj, ptr, sz, ty)
#define SWIG_NewPackedObj(ptr, sz, type)                SWIG_Python_NewPackedObj(ptr, sz, type)

/* for class or struct pointers */
#define SWIG_ConvertInstance(obj, pptr, type, flags)    SWIG_ConvertPtr(obj, pptr, type, flags)
#define SWIG_NewInstanceObj(ptr, type, flags)           SWIG_NewPointerObj(ptr, type, flags)

/* for C or C++ function pointers */
#define SWIG_ConvertFunctionPtr(obj, pptr, type)        SWIG_Python_ConvertFunctionPtr(obj, pptr, type)
#define SWIG_NewFunctionPtrObj(ptr, type)               SWIG_Python_NewPointerObj(NULL, ptr, type, 0)

/* for C++ member pointers, ie, member methods */
#define SWIG_ConvertMember(obj, ptr, sz, ty)            SWIG_Python_ConvertPacked(obj, ptr, sz, ty)
#define SWIG_NewMemberObj(ptr, sz, type)                SWIG_Python_NewPackedObj(ptr, sz, type)


/* Runtime API */

#define SWIG_GetModule(clientdata)                      SWIG_Python_GetModule(clientdata)
#define SWIG_SetModule(clientdata, pointer)             SWIG_Python_SetModule(pointer)
#define SWIG_NewClientData(obj)                         SwigPyClientData_New(obj)

#define SWIG_SetErrorObj                                SWIG_Python_SetErrorObj                            
#define SWIG_SetErrorMsg                        	SWIG_Python_SetErrorMsg				   
#define SWIG_ErrorType(code)                    	SWIG_Python_ErrorType(code)                        
#define SWIG_Error(code, msg)            		SWIG_Python_SetErrorMsg(SWIG_ErrorType(code), msg) 
#define SWIG_fail                        		goto fail					   


/* Runtime API implementation */

/* Error manipulation */

SWIGINTERN void 
SWIG_Python_SetErrorObj(PyObject *errtype, PyObject *obj) {
  SWIG_PYTHON_THREAD_BEGIN_BLOCK; 
  PyErr_SetObject(errtype, obj);
  Py_DECREF(obj);
  SWIG_PYTHON_THREAD_END_BLOCK;
}

SWIGINTERN void 
SWIG_Python_SetErrorMsg(PyObject *errtype, const char *msg) {
  SWIG_PYTHON_THREAD_BEGIN_BLOCK;
  PyErr_SetString(errtype, msg);
  SWIG_PYTHON_THREAD_END_BLOCK;
}

#define SWIG_Python_Raise(obj, type, desc)  SWIG_Python_SetErrorObj(SWIG_Python_ExceptionType(desc), obj)

/* Set a constant value */

#if defined(SWIGPYTHON_BUILTIN)

SWIGINTERN void
SwigPyBuiltin_AddPublicSymbol(PyObject *seq, const char *key) {
  PyObject *s = PyString_InternFromString(key);
  PyList_Append(seq, s);
  Py_DECREF(s);
}

SWIGINTERN void
SWIG_Python_SetConstant(PyObject *d, PyObject *public_interface, const char *name, PyObject *obj) {   
  PyDict_SetItemString(d, name, obj);
  Py_DECREF(obj);
  if (public_interface)
    SwigPyBuiltin_AddPublicSymbol(public_interface, name);
}

#else

SWIGINTERN void
SWIG_Python_SetConstant(PyObject *d, const char *name, PyObject *obj) {   
  PyDict_SetItemString(d, name, obj);
  Py_DECREF(obj);                            
}

#endif

/* Append a value to the result obj */

SWIGINTERN PyObject*
SWIG_Python_AppendOutput(PyObject* result, PyObject* obj) {
  if (!result) {
    result = obj;
  } else if (result == Py_None) {
    Py_DECREF(result);
    result = obj;
  } else {
    if (!PyList_Check(result)) {
      PyObject *o2 = result;
      result = PyList_New(1);
      PyList_SetItem(result, 0, o2);
    }
    PyList_Append(result,obj);
    Py_DECREF(obj);
  }
  return result;
}

/* Unpack the argument tuple */

SWIGINTERN Py_ssize_t
SWIG_Python_UnpackTuple(PyObject *args, const char *name, Py_ssize_t min, Py_ssize_t max, PyObject **objs)
{
  if (!args) {
    if (!min && !max) {
      return 1;
    } else {
      PyErr_Format(PyExc_TypeError, "%s expected %s%d arguments, got none", 
		   name, (min == max ? "" : "at least "), (int)min);
      return 0;
    }
  }  
  if (!PyTuple_Check(args)) {
    if (min <= 1 && max >= 1) {
      Py_ssize_t i;
      objs[0] = args;
      for (i = 1; i < max; ++i) {
	objs[i] = 0;
      }
      return 2;
    }
    PyErr_SetString(PyExc_SystemError, "UnpackTuple() argument list is not a tuple");
    return 0;
  } else {
    Py_ssize_t l = PyTuple_GET_SIZE(args);
    if (l < min) {
      PyErr_Format(PyExc_TypeError, "%s expected %s%d arguments, got %d", 
		   name, (min == max ? "" : "at least "), (int)min, (int)l);
      return 0;
    } else if (l > max) {
      PyErr_Format(PyExc_TypeError, "%s expected %s%d arguments, got %d", 
		   name, (min == max ? "" : "at most "), (int)max, (int)l);
      return 0;
    } else {
      Py_ssize_t i;
      for (i = 0; i < l; ++i) {
	objs[i] = PyTuple_GET_ITEM(args, i);
      }
      for (; l < max; ++l) {
	objs[l] = 0;
      }
      return i + 1;
    }    
  }
}

SWIGINTERN int
SWIG_Python_CheckNoKeywords(PyObject *kwargs, const char *name) {
  int no_kwargs = 1;
  if (kwargs) {
    assert(PyDict_Check(kwargs));
    if (PyDict_Size(kwargs) > 0) {
      PyErr_Format(PyExc_TypeError, "%s() does not take keyword arguments", name);
      no_kwargs = 0;
    }
  }
  return no_kwargs;
}

/* A functor is a function object with one single object argument */
#define SWIG_Python_CallFunctor(functor, obj)	        PyObject_CallFunctionObjArgs(functor, obj, NULL);

/*
  Helper for static pointer initialization for both C and C++ code, for example
  static PyObject *SWIG_STATIC_POINTER(MyVar) = NewSomething(...);
*/
#ifdef __cplusplus
#define SWIG_STATIC_POINTER(var)  var
#else
#define SWIG_STATIC_POINTER(var)  var = 0; if (!var) var
#endif

/* -----------------------------------------------------------------------------
 * Pointer declarations
 * ----------------------------------------------------------------------------- */

/* Flags for new pointer objects */
#define SWIG_POINTER_NOSHADOW       (SWIG_POINTER_OWN      << 1)
#define SWIG_POINTER_NEW            (SWIG_POINTER_NOSHADOW | SWIG_POINTER_OWN)

#define SWIG_POINTER_IMPLICIT_CONV  (SWIG_POINTER_DISOWN   << 1)

#define SWIG_BUILTIN_TP_INIT	    (SWIG_POINTER_OWN << 2)
#define SWIG_BUILTIN_INIT	    (SWIG_BUILTIN_TP_INIT | SWIG_POINTER_OWN)

#ifdef __cplusplus
extern "C" {
#endif

/* The python void return value */

SWIGRUNTIMEINLINE PyObject * 
SWIG_Py_Void(void)
{
  PyObject *none = Py_None;
  Py_INCREF(none);
  return none;
}

/* SwigPyClientData */

typedef struct {
  PyObject *klass;
  PyObject *newraw;
  PyObject *newargs;
  PyObject *destroy;
  int delargs;
  int implicitconv;
  PyTypeObject *pytype;
} SwigPyClientData;

SWIGRUNTIMEINLINE int 
SWIG_Python_CheckImplicit(swig_type_info *ty)
{
  SwigPyClientData *data = (SwigPyClientData *)ty->clientdata;
  int fail = data ? data->implicitconv : 0;
  if (fail)
    PyErr_SetString(PyExc_TypeError, "Implicit conversion is prohibited for explicit constructors.");
  return fail;
}

SWIGRUNTIMEINLINE PyObject *
SWIG_Python_ExceptionType(swig_type_info *desc) {
  SwigPyClientData *data = desc ? (SwigPyClientData *) desc->clientdata : 0;
  PyObject *klass = data ? data->klass : 0;
  return (klass ? klass : PyExc_RuntimeError);
}


SWIGRUNTIME SwigPyClientData * 
SwigPyClientData_New(PyObject* obj)
{
  if (!obj) {
    return 0;
  } else {
    SwigPyClientData *data = (SwigPyClientData *)malloc(sizeof(SwigPyClientData));
    /* the klass element */
    data->klass = obj;
    Py_INCREF(data->klass);
    /* the newraw method and newargs arguments used to create a new raw instance */
    if (PyClass_Check(obj)) {
      data->newraw = 0;
      data->newargs = obj;
      Py_INCREF(obj);
    } else {
      data->newraw = PyObject_GetAttrString(data->klass, "__new__");
      if (data->newraw) {
	Py_INCREF(data->newraw);
	data->newargs = PyTuple_New(1);
	PyTuple_SetItem(data->newargs, 0, obj);
      } else {
	data->newargs = obj;
      }
      Py_INCREF(data->newargs);
    }
    /* the destroy method, aka as the C++ delete method */
    data->destroy = PyObject_GetAttrString(data->klass, "__swig_destroy__");
    if (PyErr_Occurred()) {
      PyErr_Clear();
      data->destroy = 0;
    }
    if (data->destroy) {
      int flags;
      Py_INCREF(data->destroy);
      flags = PyCFunction_GET_FLAGS(data->destroy);
      data->delargs = !(flags & (METH_O));
    } else {
      data->delargs = 0;
    }
    data->implicitconv = 0;
    data->pytype = 0;
    return data;
  }
}

SWIGRUNTIME void 
SwigPyClientData_Del(SwigPyClientData *data) {
  Py_XDECREF(data->newraw);
  Py_XDECREF(data->newargs);
  Py_XDECREF(data->destroy);
}

/* =============== SwigPyObject =====================*/

typedef struct {
  PyObject_HEAD
  void *ptr;
  swig_type_info *ty;
  int own;
  PyObject *next;
#ifdef SWIGPYTHON_BUILTIN
  PyObject *dict;
#endif
} SwigPyObject;


#ifdef SWIGPYTHON_BUILTIN

SWIGRUNTIME PyObject *
SwigPyObject_get___dict__(PyObject *v, PyObject *SWIGUNUSEDPARM(args))
{
  SwigPyObject *sobj = (SwigPyObject *)v;

  if (!sobj->dict)
    sobj->dict = PyDict_New();

  Py_INCREF(sobj->dict);
  return sobj->dict;
}

#endif

SWIGRUNTIME PyObject *
SwigPyObject_long(SwigPyObject *v)
{
  return PyLong_FromVoidPtr(v->ptr);
}

SWIGRUNTIME PyObject *
SwigPyObject_format(const char* fmt, SwigPyObject *v)
{
  PyObject *res = NULL;
  PyObject *args = PyTuple_New(1);
  if (args) {
    if (PyTuple_SetItem(args, 0, SwigPyObject_long(v)) == 0) {
      PyObject *ofmt = SWIG_Python_str_FromChar(fmt);
      if (ofmt) {
#if PY_VERSION_HEX >= 0x03000000
	res = PyUnicode_Format(ofmt,args);
#else
	res = PyString_Format(ofmt,args);
#endif
	Py_DECREF(ofmt);
      }
      Py_DECREF(args);
    }
  }
  return res;
}

SWIGRUNTIME PyObject *
SwigPyObject_oct(SwigPyObject *v)
{
  return SwigPyObject_format("%o",v);
}

SWIGRUNTIME PyObject *
SwigPyObject_hex(SwigPyObject *v)
{
  return SwigPyObject_format("%x",v);
}

SWIGRUNTIME PyObject *
SwigPyObject_repr(SwigPyObject *v)
{
  const char *name = SWIG_TypePrettyName(v->ty);
  PyObject *repr = SWIG_Python_str_FromFormat("<Swig Object of type '%s' at %p>", (name ? name : "unknown"), (void *)v);
  if (v->next) {
    PyObject *nrep = SwigPyObject_repr((SwigPyObject *)v->next);
# if PY_VERSION_HEX >= 0x03000000
    PyObject *joined = PyUnicode_Concat(repr, nrep);
    Py_DecRef(repr);
    Py_DecRef(nrep);
    repr = joined;
# else
    PyString_ConcatAndDel(&repr,nrep);
# endif
  }
  return repr;  
}

/* We need a version taking two PyObject* parameters so it's a valid
 * PyCFunction to use in swigobject_methods[]. */
SWIGRUNTIME PyObject *
SwigPyObject_repr2(PyObject *v, PyObject *SWIGUNUSEDPARM(args))
{
  return SwigPyObject_repr((SwigPyObject*)v);
}

SWIGRUNTIME int
SwigPyObject_compare(SwigPyObject *v, SwigPyObject *w)
{
  void *i = v->ptr;
  void *j = w->ptr;
  return (i < j) ? -1 : ((i > j) ? 1 : 0);
}

/* Added for Python 3.x, would it also be useful for Python 2.x? */
SWIGRUNTIME PyObject*
SwigPyObject_richcompare(SwigPyObject *v, SwigPyObject *w, int op)
{
  PyObject* res;
  if( op != Py_EQ && op != Py_NE ) {
    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
  }
  res = PyBool_FromLong( (SwigPyObject_compare(v, w)==0) == (op == Py_EQ) ? 1 : 0);
  return res;  
}


SWIGRUNTIME PyTypeObject* SwigPyObject_TypeOnce(void);

#ifdef SWIGPYTHON_BUILTIN
static swig_type_info *SwigPyObject_stype = 0;
SWIGRUNTIME PyTypeObject*
SwigPyObject_type(void) {
    SwigPyClientData *cd;
    assert(SwigPyObject_stype);
    cd = (SwigPyClientData*) SwigPyObject_stype->clientdata;
    assert(cd);
    assert(cd->pytype);
    return cd->pytype;
}
#else
SWIGRUNTIME PyTypeObject*
SwigPyObject_type(void) {
  static PyTypeObject *SWIG_STATIC_POINTER(type) = SwigPyObject_TypeOnce();
  return type;
}
#endif

SWIGRUNTIMEINLINE int
SwigPyObject_Check(PyObject *op) {
#ifdef SWIGPYTHON_BUILTIN
  PyTypeObject *target_tp = SwigPyObject_type();
  if (PyType_IsSubtype(op->ob_type, target_tp))
    return 1;
  return (strcmp(op->ob_type->tp_name, "SwigPyObject") == 0);
#else
  return (Py_TYPE(op) == SwigPyObject_type())
    || (strcmp(Py_TYPE(op)->tp_name,"SwigPyObject") == 0);
#endif
}

SWIGRUNTIME PyObject *
SwigPyObject_New(void *ptr, swig_type_info *ty, int own);

SWIGRUNTIME void
SwigPyObject_dealloc(PyObject *v)
{
  SwigPyObject *sobj = (SwigPyObject *) v;
  PyObject *next = sobj->next;
  if (sobj->own == SWIG_POINTER_OWN) {
    swig_type_info *ty = sobj->ty;
    SwigPyClientData *data = ty ? (SwigPyClientData *) ty->clientdata : 0;
    PyObject *destroy = data ? data->destroy : 0;
    if (destroy) {
      /* destroy is always a VARARGS method */
      PyObject *res;

      /* PyObject_CallFunction() has the potential to silently drop
         the active exception.  In cases of unnamed temporary
         variable or where we just finished iterating over a generator
         StopIteration will be active right now, and this needs to
         remain true upon return from SwigPyObject_dealloc.  So save
         and restore. */
      
      PyObject *type = NULL, *value = NULL, *traceback = NULL;
      PyErr_Fetch(&type, &value, &traceback);

      if (data->delargs) {
        /* we need to create a temporary object to carry the destroy operation */
        PyObject *tmp = SwigPyObject_New(sobj->ptr, ty, 0);
        res = SWIG_Python_CallFunctor(destroy, tmp);
        Py_DECREF(tmp);
      } else {
        PyCFunction meth = PyCFunction_GET_FUNCTION(destroy);
        PyObject *mself = PyCFunction_GET_SELF(destroy);
        res = ((*meth)(mself, v));
      }
      if (!res)
        PyErr_WriteUnraisable(destroy);

      PyErr_Restore(type, value, traceback);

      Py_XDECREF(res);
    } 
#if !defined(SWIG_PYTHON_SILENT_MEMLEAK)
    else {
      const char *name = SWIG_TypePrettyName(ty);
      printf("swig/python detected a memory leak of type '%s', no destructor found.\n", (name ? name : "unknown"));
    }
#endif
  } 
  Py_XDECREF(next);
  PyObject_DEL(v);
}

SWIGRUNTIME PyObject* 
SwigPyObject_append(PyObject* v, PyObject* next)
{
  SwigPyObject *sobj = (SwigPyObject *) v;
  if (!SwigPyObject_Check(next)) {
    PyErr_SetString(PyExc_TypeError, "Attempt to append a non SwigPyObject");
    return NULL;
  }
  sobj->next = next;
  Py_INCREF(next);
  return SWIG_Py_Void();
}

SWIGRUNTIME PyObject* 
SwigPyObject_next(PyObject* v, PyObject *SWIGUNUSEDPARM(args))
{
  SwigPyObject *sobj = (SwigPyObject *) v;
  if (sobj->next) {    
    Py_INCREF(sobj->next);
    return sobj->next;
  } else {
    return SWIG_Py_Void();
  }
}

SWIGINTERN PyObject*
SwigPyObject_disown(PyObject* v, PyObject *SWIGUNUSEDPARM(args))
{
  SwigPyObject *sobj = (SwigPyObject *)v;
  sobj->own = 0;
  return SWIG_Py_Void();
}

SWIGINTERN PyObject*
SwigPyObject_acquire(PyObject* v, PyObject *SWIGUNUSEDPARM(args))
{
  SwigPyObject *sobj = (SwigPyObject *)v;
  sobj->own = SWIG_POINTER_OWN;
  return SWIG_Py_Void();
}

SWIGINTERN PyObject*
SwigPyObject_own(PyObject *v, PyObject *args)
{
  PyObject *val = 0;
  if (!PyArg_UnpackTuple(args, "own", 0, 1, &val)) {
    return NULL;
  } else {
    SwigPyObject *sobj = (SwigPyObject *)v;
    PyObject *obj = PyBool_FromLong(sobj->own);
    if (val) {
      if (PyObject_IsTrue(val)) {
        SwigPyObject_acquire(v,args);
      } else {
        SwigPyObject_disown(v,args);
      }
    } 
    return obj;
  }
}

static PyMethodDef
swigobject_methods[] = {
  {"disown",  SwigPyObject_disown,  METH_NOARGS,  "releases ownership of the pointer"},
  {"acquire", SwigPyObject_acquire, METH_NOARGS,  "acquires ownership of the pointer"},
  {"own",     SwigPyObject_own,     METH_VARARGS, "returns/sets ownership of the pointer"},
  {"append",  SwigPyObject_append,  METH_O,       "appends another 'this' object"},
  {"next",    SwigPyObject_next,    METH_NOARGS,  "returns the next 'this' object"},
  {"__repr__",SwigPyObject_repr2,   METH_NOARGS,  "returns object representation"},
  {0, 0, 0, 0}  
};

SWIGRUNTIME PyTypeObject*
SwigPyObject_TypeOnce(void) {
  static char swigobject_doc[] = "Swig object carries a C/C++ instance pointer";

  static PyNumberMethods SwigPyObject_as_number = {
    (binaryfunc)0, /*nb_add*/
    (binaryfunc)0, /*nb_subtract*/
    (binaryfunc)0, /*nb_multiply*/
    /* nb_divide removed in Python 3 */
#if PY_VERSION_HEX < 0x03000000
    (binaryfunc)0, /*nb_divide*/
#endif
    (binaryfunc)0, /*nb_remainder*/
    (binaryfunc)0, /*nb_divmod*/
    (ternaryfunc)0,/*nb_power*/
    (unaryfunc)0,  /*nb_negative*/
    (unaryfunc)0,  /*nb_positive*/
    (unaryfunc)0,  /*nb_absolute*/
    (inquiry)0,    /*nb_nonzero*/
    0,		   /*nb_invert*/
    0,		   /*nb_lshift*/
    0,		   /*nb_rshift*/
    0,		   /*nb_and*/
    0,		   /*nb_xor*/
    0,		   /*nb_or*/
#if PY_VERSION_HEX < 0x03000000
    0,   /*nb_coerce*/
#endif
    (unaryfunc)SwigPyObject_long, /*nb_int*/
#if PY_VERSION_HEX < 0x03000000
    (unaryfunc)SwigPyObject_long, /*nb_long*/
#else
    0, /*nb_reserved*/
#endif
    (unaryfunc)0,                 /*nb_float*/
#if PY_VERSION_HEX < 0x03000000
    (unaryfunc)SwigPyObject_oct,  /*nb_oct*/
    (unaryfunc)SwigPyObject_hex,  /*nb_hex*/
#endif
#if PY_VERSION_HEX >= 0x03050000 /* 3.5 */
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /* nb_inplace_add -> nb_inplace_matrix_multiply */
#elif PY_VERSION_HEX >= 0x03000000 /* 3.0 */
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /* nb_inplace_add -> nb_index, nb_inplace_divide removed */
#else
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /* nb_inplace_add -> nb_index */
#endif
  };

  static PyTypeObject swigpyobject_type;
  static int type_init = 0;
  if (!type_init) {
    const PyTypeObject tmp = {
#if PY_VERSION_HEX >= 0x03000000
      PyVarObject_HEAD_INIT(NULL, 0)
#else
      PyObject_HEAD_INIT(NULL)
      0,                                    /* ob_size */
#endif
      "SwigPyObject",                       /* tp_name */
      sizeof(SwigPyObject),                 /* tp_basicsize */
      0,                                    /* tp_itemsize */
      (destructor)SwigPyObject_dealloc,     /* tp_dealloc */
      0,                                    /* tp_print */
      (getattrfunc)0,                       /* tp_getattr */
      (setattrfunc)0,                       /* tp_setattr */
#if PY_VERSION_HEX >= 0x03000000
      0, /* tp_reserved in 3.0.1, tp_compare in 3.0.0 but not used */
#else
      (cmpfunc)SwigPyObject_compare,        /* tp_compare */
#endif
      (reprfunc)SwigPyObject_repr,          /* tp_repr */
      &SwigPyObject_as_number,              /* tp_as_number */
      0,                                    /* tp_as_sequence */
      0,                                    /* tp_as_mapping */
      (hashfunc)0,                          /* tp_hash */
      (ternaryfunc)0,                       /* tp_call */
      0,                                    /* tp_str */
      PyObject_GenericGetAttr,              /* tp_getattro */
      0,                                    /* tp_setattro */
      0,                                    /* tp_as_buffer */
      Py_TPFLAGS_DEFAULT,                   /* tp_flags */
      swigobject_doc,                       /* tp_doc */
      0,                                    /* tp_traverse */
      0,                                    /* tp_clear */
      (richcmpfunc)SwigPyObject_richcompare,/* tp_richcompare */
      0,                                    /* tp_weaklistoffset */
      0,                                    /* tp_iter */
      0,                                    /* tp_iternext */
      swigobject_methods,                   /* tp_methods */
      0,                                    /* tp_members */
      0,                                    /* tp_getset */
      0,                                    /* tp_base */
      0,                                    /* tp_dict */
      0,                                    /* tp_descr_get */
      0,                                    /* tp_descr_set */
      0,                                    /* tp_dictoffset */
      0,                                    /* tp_init */
      0,                                    /* tp_alloc */
      0,                                    /* tp_new */
      0,                                    /* tp_free */
      0,                                    /* tp_is_gc */
      0,                                    /* tp_bases */
      0,                                    /* tp_mro */
      0,                                    /* tp_cache */
      0,                                    /* tp_subclasses */
      0,                                    /* tp_weaklist */
      0,                                    /* tp_del */
      0,                                    /* tp_version_tag */
#if PY_VERSION_HEX >= 0x03040000
      0,                                    /* tp_finalize */
#endif
#if PY_VERSION_HEX >= 0x03080000
      0,                                    /* tp_vectorcall */
#endif
#if (PY_VERSION_HEX >= 0x03080000) && (PY_VERSION_HEX < 0x03090000)
      0,                                    /* tp_print */
#endif
#ifdef COUNT_ALLOCS
      0,                                    /* tp_allocs */
      0,                                    /* tp_frees */
      0,                                    /* tp_maxalloc */
      0,                                    /* tp_prev */
      0                                     /* tp_next */
#endif
    };
    swigpyobject_type = tmp;
    type_init = 1;
    if (PyType_Ready(&swigpyobject_type) < 0)
      return NULL;
  }
  return &swigpyobject_type;
}

SWIGRUNTIME PyObject *
SwigPyObject_New(void *ptr, swig_type_info *ty, int own)
{
  SwigPyObject *sobj = PyObject_NEW(SwigPyObject, SwigPyObject_type());
  if (sobj) {
    sobj->ptr  = ptr;
    sobj->ty   = ty;
    sobj->own  = own;
    sobj->next = 0;
  }
  return (PyObject *)sobj;
}

/* -----------------------------------------------------------------------------
 * Implements a simple Swig Packed type, and use it instead of string
 * ----------------------------------------------------------------------------- */

typedef struct {
  PyObject_HEAD
  void *pack;
  swig_type_info *ty;
  size_t size;
} SwigPyPacked;

SWIGRUNTIME PyObject *
SwigPyPacked_repr(SwigPyPacked *v)
{
  char result[SWIG_BUFFER_SIZE];
  if (SWIG_PackDataName(result, v->pack, v->size, 0, sizeof(result))) {
    return SWIG_Python_str_FromFormat("<Swig Packed at %s%s>", result, v->ty->name);
  } else {
    return SWIG_Python_str_FromFormat("<Swig Packed %s>", v->ty->name);
  }  
}

SWIGRUNTIME PyObject *
SwigPyPacked_str(SwigPyPacked *v)
{
  char result[SWIG_BUFFER_SIZE];
  if (SWIG_PackDataName(result, v->pack, v->size, 0, sizeof(result))){
    return SWIG_Python_str_FromFormat("%s%s", result, v->ty->name);
  } else {
    return SWIG_Python_str_FromChar(v->ty->name);
  }  
}

SWIGRUNTIME int
SwigPyPacked_compare(SwigPyPacked *v, SwigPyPacked *w)
{
  size_t i = v->size;
  size_t j = w->size;
  int s = (i < j) ? -1 : ((i > j) ? 1 : 0);
  return s ? s : strncmp((const char *)v->pack, (const char *)w->pack, 2*v->size);
}

SWIGRUNTIME PyTypeObject* SwigPyPacked_TypeOnce(void);

SWIGRUNTIME PyTypeObject*
SwigPyPacked_type(void) {
  static PyTypeObject *SWIG_STATIC_POINTER(type) = SwigPyPacked_TypeOnce();
  return type;
}

SWIGRUNTIMEINLINE int
SwigPyPacked_Check(PyObject *op) {
  return ((op)->ob_type == SwigPyPacked_TypeOnce()) 
    || (strcmp((op)->ob_type->tp_name,"SwigPyPacked") == 0);
}

SWIGRUNTIME void
SwigPyPacked_dealloc(PyObject *v)
{
  if (SwigPyPacked_Check(v)) {
    SwigPyPacked *sobj = (SwigPyPacked *) v;
    free(sobj->pack);
  }
  PyObject_DEL(v);
}

SWIGRUNTIME PyTypeObject*
SwigPyPacked_TypeOnce(void) {
  static char swigpacked_doc[] = "Swig object carries a C/C++ instance pointer";
  static PyTypeObject swigpypacked_type;
  static int type_init = 0;
  if (!type_init) {
    const PyTypeObject tmp = {
#if PY_VERSION_HEX>=0x03000000
      PyVarObject_HEAD_INIT(NULL, 0)
#else
      PyObject_HEAD_INIT(NULL)
      0,                                    /* ob_size */
#endif
      "SwigPyPacked",                       /* tp_name */
      sizeof(SwigPyPacked),                 /* tp_basicsize */
      0,                                    /* tp_itemsize */
      (destructor)SwigPyPacked_dealloc,     /* tp_dealloc */
      0,                                    /* tp_print */
      (getattrfunc)0,                       /* tp_getattr */
      (setattrfunc)0,                       /* tp_setattr */
#if PY_VERSION_HEX>=0x03000000
      0, /* tp_reserved in 3.0.1 */
#else
      (cmpfunc)SwigPyPacked_compare,        /* tp_compare */
#endif
      (reprfunc)SwigPyPacked_repr,          /* tp_repr */
      0,                                    /* tp_as_number */
      0,                                    /* tp_as_sequence */
      0,                                    /* tp_as_mapping */
      (hashfunc)0,                          /* tp_hash */
      (ternaryfunc)0,                       /* tp_call */
      (reprfunc)SwigPyPacked_str,           /* tp_str */
      PyObject_GenericGetAttr,              /* tp_getattro */
      0,                                    /* tp_setattro */
      0,                                    /* tp_as_buffer */
      Py_TPFLAGS_DEFAULT,                   /* tp_flags */
      swigpacked_doc,                       /* tp_doc */
      0,                                    /* tp_traverse */
      0,                                    /* tp_clear */
      0,                                    /* tp_richcompare */
      0,                                    /* tp_weaklistoffset */
      0,                                    /* tp_iter */
      0,                                    /* tp_iternext */
      0,                                    /* tp_methods */
      0,                                    /* tp_members */
      0,                                    /* tp_getset */
      0,                                    /* tp_base */
      0,                                    /* tp_dict */
      0,                                    /* tp_descr_get */
      0,                                    /* tp_descr_set */
      0,                                    /* tp_dictoffset */
      0,                                    /* tp_init */
      0,                                    /* tp_alloc */
      0,                                    /* tp_new */
      0,                                    /* tp_free */
      0,                                    /* tp_is_gc */
      0,                                    /* tp_bases */
      0,                                    /* tp_mro */
      0,                                    /* tp_cache */
      0,                                    /* tp_subclasses */
      0,                                    /* tp_weaklist */
      0,                                    /* tp_del */
      0,                                    /* tp_version_tag */
#if PY_VERSION_HEX >= 0x03040000
      0,                                    /* tp_finalize */
#endif
#if PY_VERSION_HEX >= 0x03080000
      0,                                    /* tp_vectorcall */
#endif
#if (PY_VERSION_HEX >= 0x03080000) && (PY_VERSION_HEX < 0x03090000)
      0,                                    /* tp_print */
#endif
#ifdef COUNT_ALLOCS
      0,                                    /* tp_allocs */
      0,                                    /* tp_frees */
      0,                                    /* tp_maxalloc */
      0,                                    /* tp_prev */
      0                                     /* tp_next */
#endif
    };
    swigpypacked_type = tmp;
    type_init = 1;
    if (PyType_Ready(&swigpypacked_type) < 0)
      return NULL;
  }
  return &swigpypacked_type;
}

SWIGRUNTIME PyObject *
SwigPyPacked_New(void *ptr, size_t size, swig_type_info *ty)
{
  SwigPyPacked *sobj = PyObject_NEW(SwigPyPacked, SwigPyPacked_type());
  if (sobj) {
    void *pack = malloc(size);
    if (pack) {
      memcpy(pack, ptr, size);
      sobj->pack = pack;
      sobj->ty   = ty;
      sobj->size = size;
    } else {
      PyObject_DEL((PyObject *) sobj);
      sobj = 0;
    }
  }
  return (PyObject *) sobj;
}

SWIGRUNTIME swig_type_info *
SwigPyPacked_UnpackData(PyObject *obj, void *ptr, size_t size)
{
  if (SwigPyPacked_Check(obj)) {
    SwigPyPacked *sobj = (SwigPyPacked *)obj;
    if (sobj->size != size) return 0;
    memcpy(ptr, sobj->pack, size);
    return sobj->ty;
  } else {
    return 0;
  }
}

/* -----------------------------------------------------------------------------
 * pointers/data manipulation
 * ----------------------------------------------------------------------------- */

static PyObject *Swig_This_global = NULL;

SWIGRUNTIME PyObject *
SWIG_This(void)
{
  if (Swig_This_global == NULL)
    Swig_This_global = SWIG_Python_str_FromChar("this");
  return Swig_This_global;
}

/* #define SWIG_PYTHON_SLOW_GETSET_THIS */

/* TODO: I don't know how to implement the fast getset in Python 3 right now */
#if PY_VERSION_HEX>=0x03000000
#define SWIG_PYTHON_SLOW_GETSET_THIS 
#endif

SWIGRUNTIME SwigPyObject *
SWIG_Python_GetSwigThis(PyObject *pyobj) 
{
  PyObject *obj;

  if (SwigPyObject_Check(pyobj))
    return (SwigPyObject *) pyobj;

#ifdef SWIGPYTHON_BUILTIN
  (void)obj;
# ifdef PyWeakref_CheckProxy
  if (PyWeakref_CheckProxy(pyobj)) {
    pyobj = PyWeakref_GET_OBJECT(pyobj);
    if (pyobj && SwigPyObject_Check(pyobj))
      return (SwigPyObject*) pyobj;
  }
# endif
  return NULL;
#else

  obj = 0;

#if !defined(SWIG_PYTHON_SLOW_GETSET_THIS)
  if (PyInstance_Check(pyobj)) {
    obj = _PyInstance_Lookup(pyobj, SWIG_This());      
  } else {
    PyObject **dictptr = _PyObject_GetDictPtr(pyobj);
    if (dictptr != NULL) {
      PyObject *dict = *dictptr;
      obj = dict ? PyDict_GetItem(dict, SWIG_This()) : 0;
    } else {
#ifdef PyWeakref_CheckProxy
      if (PyWeakref_CheckProxy(pyobj)) {
	PyObject *wobj = PyWeakref_GET_OBJECT(pyobj);
	return wobj ? SWIG_Python_GetSwigThis(wobj) : 0;
      }
#endif
      obj = PyObject_GetAttr(pyobj,SWIG_This());
      if (obj) {
	Py_DECREF(obj);
      } else {
	if (PyErr_Occurred()) PyErr_Clear();
	return 0;
      }
    }
  }
#else
  obj = PyObject_GetAttr(pyobj,SWIG_This());
  if (obj) {
    Py_DECREF(obj);
  } else {
    if (PyErr_Occurred()) PyErr_Clear();
    return 0;
  }
#endif
  if (obj && !SwigPyObject_Check(obj)) {
    /* a PyObject is called 'this', try to get the 'real this'
       SwigPyObject from it */ 
    return SWIG_Python_GetSwigThis(obj);
  }
  return (SwigPyObject *)obj;
#endif
}

/* Acquire a pointer value */

SWIGRUNTIME int
SWIG_Python_AcquirePtr(PyObject *obj, int own) {
  if (own == SWIG_POINTER_OWN) {
    SwigPyObject *sobj = SWIG_Python_GetSwigThis(obj);
    if (sobj) {
      int oldown = sobj->own;
      sobj->own = own;
      return oldown;
    }
  }
  return 0;
}

/* Convert a pointer value */

SWIGRUNTIME int
SWIG_Python_ConvertPtrAndOwn(PyObject *obj, void **ptr, swig_type_info *ty, int flags, int *own) {
  int res;
  SwigPyObject *sobj;
  int implicit_conv = (flags & SWIG_POINTER_IMPLICIT_CONV) != 0;

  if (!obj)
    return SWIG_ERROR;
  if (obj == Py_None && !implicit_conv) {
    if (ptr)
      *ptr = 0;
    return (flags & SWIG_POINTER_NO_NULL) ? SWIG_NullReferenceError : SWIG_OK;
  }

  res = SWIG_ERROR;

  sobj = SWIG_Python_GetSwigThis(obj);
  if (own)
    *own = 0;
  while (sobj) {
    void *vptr = sobj->ptr;
    if (ty) {
      swig_type_info *to = sobj->ty;
      if (to == ty) {
        /* no type cast needed */
        if (ptr) *ptr = vptr;
        break;
      } else {
        swig_cast_info *tc = SWIG_TypeCheck(to->name,ty);
        if (!tc) {
          sobj = (SwigPyObject *)sobj->next;
        } else {
          if (ptr) {
            int newmemory = 0;
            *ptr = SWIG_TypeCast(tc,vptr,&newmemory);
            if (newmemory == SWIG_CAST_NEW_MEMORY) {
              assert(own); /* badly formed typemap which will lead to a memory leak - it must set and use own to delete *ptr */
              if (own)
                *own = *own | SWIG_CAST_NEW_MEMORY;
            }
          }
          break;
        }
      }
    } else {
      if (ptr) *ptr = vptr;
      break;
    }
  }
  if (sobj) {
    if (own)
      *own = *own | sobj->own;
    if (flags & SWIG_POINTER_DISOWN) {
      sobj->own = 0;
    }
    res = SWIG_OK;
  } else {
    if (implicit_conv) {
      SwigPyClientData *data = ty ? (SwigPyClientData *) ty->clientdata : 0;
      if (data && !data->implicitconv) {
        PyObject *klass = data->klass;
        if (klass) {
          PyObject *impconv;
          data->implicitconv = 1; /* avoid recursion and call 'explicit' constructors*/
          impconv = SWIG_Python_CallFunctor(klass, obj);
          data->implicitconv = 0;
          if (PyErr_Occurred()) {
            PyErr_Clear();
            impconv = 0;
          }
          if (impconv) {
            SwigPyObject *iobj = SWIG_Python_GetSwigThis(impconv);
            if (iobj) {
              void *vptr;
              res = SWIG_Python_ConvertPtrAndOwn((PyObject*)iobj, &vptr, ty, 0, 0);
              if (SWIG_IsOK(res)) {
                if (ptr) {
                  *ptr = vptr;
                  /* transfer the ownership to 'ptr' */
                  iobj->own = 0;
                  res = SWIG_AddCast(res);
                  res = SWIG_AddNewMask(res);
                } else {
                  res = SWIG_AddCast(res);		    
                }
              }
            }
            Py_DECREF(impconv);
          }
        }
      }
      if (!SWIG_IsOK(res) && obj == Py_None) {
        if (ptr)
          *ptr = 0;
        if (PyErr_Occurred())
          PyErr_Clear();
        res = SWIG_OK;
      }
    }
  }
  return res;
}

/* Convert a function ptr value */

SWIGRUNTIME int
SWIG_Python_ConvertFunctionPtr(PyObject *obj, void **ptr, swig_type_info *ty) {
  if (!PyCFunction_Check(obj)) {
    return SWIG_ConvertPtr(obj, ptr, ty, 0);
  } else {
    void *vptr = 0;
    swig_cast_info *tc;

    /* here we get the method pointer for callbacks */
    const char *doc = (((PyCFunctionObject *)obj) -> m_ml -> ml_doc);
    const char *desc = doc ? strstr(doc, "swig_ptr: ") : 0;
    if (desc)
      desc = ty ? SWIG_UnpackVoidPtr(desc + 10, &vptr, ty->name) : 0;
    if (!desc)
      return SWIG_ERROR;
    tc = SWIG_TypeCheck(desc,ty);
    if (tc) {
      int newmemory = 0;
      *ptr = SWIG_TypeCast(tc,vptr,&newmemory);
      assert(!newmemory); /* newmemory handling not yet implemented */
    } else {
      return SWIG_ERROR;
    }
    return SWIG_OK;
  }
}

/* Convert a packed pointer value */

SWIGRUNTIME int
SWIG_Python_ConvertPacked(PyObject *obj, void *ptr, size_t sz, swig_type_info *ty) {
  swig_type_info *to = SwigPyPacked_UnpackData(obj, ptr, sz);
  if (!to) return SWIG_ERROR;
  if (ty) {
    if (to != ty) {
      /* check type cast? */
      swig_cast_info *tc = SWIG_TypeCheck(to->name,ty);
      if (!tc) return SWIG_ERROR;
    }
  }
  return SWIG_OK;
}  

/* -----------------------------------------------------------------------------
 * Create a new pointer object
 * ----------------------------------------------------------------------------- */

/*
  Create a new instance object, without calling __init__, and set the
  'this' attribute.
*/

SWIGRUNTIME PyObject* 
SWIG_Python_NewShadowInstance(SwigPyClientData *data, PyObject *swig_this)
{
  PyObject *inst = 0;
  PyObject *newraw = data->newraw;
  if (newraw) {
    inst = PyObject_Call(newraw, data->newargs, NULL);
    if (inst) {
#if !defined(SWIG_PYTHON_SLOW_GETSET_THIS)
      PyObject **dictptr = _PyObject_GetDictPtr(inst);
      if (dictptr != NULL) {
	PyObject *dict = *dictptr;
	if (dict == NULL) {
	  dict = PyDict_New();
	  *dictptr = dict;
	  PyDict_SetItem(dict, SWIG_This(), swig_this);
	}
      }
#else
      if (PyObject_SetAttr(inst, SWIG_This(), swig_this) == -1) {
        Py_DECREF(inst);
        inst = 0;
      }
#endif
    }
  } else {
#if PY_VERSION_HEX >= 0x03000000
    PyObject *empty_args = PyTuple_New(0);
    if (empty_args) {
      PyObject *empty_kwargs = PyDict_New();
      if (empty_kwargs) {
        inst = ((PyTypeObject *)data->newargs)->tp_new((PyTypeObject *)data->newargs, empty_args, empty_kwargs);
        Py_DECREF(empty_kwargs);
        if (inst) {
          if (PyObject_SetAttr(inst, SWIG_This(), swig_this) == -1) {
            Py_DECREF(inst);
            inst = 0;
          } else {
            Py_TYPE(inst)->tp_flags &= ~Py_TPFLAGS_VALID_VERSION_TAG;
          }
        }
      }
      Py_DECREF(empty_args);
    }
#else
    PyObject *dict = PyDict_New();
    if (dict) {
      PyDict_SetItem(dict, SWIG_This(), swig_this);
      inst = PyInstance_NewRaw(data->newargs, dict);
      Py_DECREF(dict);
    }
#endif
  }
  return inst;
}

SWIGRUNTIME int
SWIG_Python_SetSwigThis(PyObject *inst, PyObject *swig_this)
{
#if !defined(SWIG_PYTHON_SLOW_GETSET_THIS)
  PyObject **dictptr = _PyObject_GetDictPtr(inst);
  if (dictptr != NULL) {
    PyObject *dict = *dictptr;
    if (dict == NULL) {
      dict = PyDict_New();
      *dictptr = dict;
    }
    return PyDict_SetItem(dict, SWIG_This(), swig_this);
  }
#endif
  return PyObject_SetAttr(inst, SWIG_This(), swig_this);
} 


SWIGINTERN PyObject *
SWIG_Python_InitShadowInstance(PyObject *args) {
  PyObject *obj[2];
  if (!SWIG_Python_UnpackTuple(args, "swiginit", 2, 2, obj)) {
    return NULL;
  } else {
    SwigPyObject *sthis = SWIG_Python_GetSwigThis(obj[0]);
    if (sthis) {
      SwigPyObject_append((PyObject*) sthis, obj[1]);
    } else {
      if (SWIG_Python_SetSwigThis(obj[0], obj[1]) != 0)
        return NULL;
    }
    return SWIG_Py_Void();
  }
}

/* Create a new pointer object */

SWIGRUNTIME PyObject *
SWIG_Python_NewPointerObj(PyObject *self, void *ptr, swig_type_info *type, int flags) {
  SwigPyClientData *clientdata;
  PyObject * robj;
  int own;

  if (!ptr)
    return SWIG_Py_Void();

  clientdata = type ? (SwigPyClientData *)(type->clientdata) : 0;
  own = (flags & SWIG_POINTER_OWN) ? SWIG_POINTER_OWN : 0;
  if (clientdata && clientdata->pytype) {
    SwigPyObject *newobj;
    if (flags & SWIG_BUILTIN_TP_INIT) {
      newobj = (SwigPyObject*) self;
      if (newobj->ptr) {
        PyObject *next_self = clientdata->pytype->tp_alloc(clientdata->pytype, 0);
        while (newobj->next)
	  newobj = (SwigPyObject *) newobj->next;
        newobj->next = next_self;
        newobj = (SwigPyObject *)next_self;
#ifdef SWIGPYTHON_BUILTIN
        newobj->dict = 0;
#endif
      }
    } else {
      newobj = PyObject_New(SwigPyObject, clientdata->pytype);
#ifdef SWIGPYTHON_BUILTIN
      newobj->dict = 0;
#endif
    }
    if (newobj) {
      newobj->ptr = ptr;
      newobj->ty = type;
      newobj->own = own;
      newobj->next = 0;
      return (PyObject*) newobj;
    }
    return SWIG_Py_Void();
  }

  assert(!(flags & SWIG_BUILTIN_TP_INIT));

  robj = SwigPyObject_New(ptr, type, own);
  if (robj && clientdata && !(flags & SWIG_POINTER_NOSHADOW)) {
    PyObject *inst = SWIG_Python_NewShadowInstance(clientdata, robj);
    Py_DECREF(robj);
    robj = inst;
  }
  return robj;
}

/* Create a new packed object */

SWIGRUNTIMEINLINE PyObject *
SWIG_Python_NewPackedObj(void *ptr, size_t sz, swig_type_info *type) {
  return ptr ? SwigPyPacked_New((void *) ptr, sz, type) : SWIG_Py_Void();
}

/* -----------------------------------------------------------------------------*
 *  Get type list 
 * -----------------------------------------------------------------------------*/

#ifdef SWIG_LINK_RUNTIME
void *SWIG_ReturnGlobalTypeList(void *);
#endif

SWIGRUNTIME swig_module_info *
SWIG_Python_GetModule(void *SWIGUNUSEDPARM(clientdata)) {
  static void *type_pointer = (void *)0;
  /* first check if module already created */
  if (!type_pointer) {
#ifdef SWIG_LINK_RUNTIME
    type_pointer = SWIG_ReturnGlobalTypeList((void *)0);
#else
    type_pointer = PyCapsule_Import(SWIGPY_CAPSULE_NAME, 0);
    if (PyErr_Occurred()) {
      PyErr_Clear();
      type_pointer = (void *)0;
    }
#endif
  }
  return (swig_module_info *) type_pointer;
}

SWIGRUNTIME void
SWIG_Python_DestroyModule(PyObject *obj)
{
  swig_module_info *swig_module = (swig_module_info *) PyCapsule_GetPointer(obj, SWIGPY_CAPSULE_NAME);
  swig_type_info **types = swig_module->types;
  size_t i;
  for (i =0; i < swig_module->size; ++i) {
    swig_type_info *ty = types[i];
    if (ty->owndata) {
      SwigPyClientData *data = (SwigPyClientData *) ty->clientdata;
      if (data) SwigPyClientData_Del(data);
    }
  }
  Py_DECREF(SWIG_This());
  Swig_This_global = NULL;
}

SWIGRUNTIME void
SWIG_Python_SetModule(swig_module_info *swig_module) {
#if PY_VERSION_HEX >= 0x03000000
 /* Add a dummy module object into sys.modules */
  PyObject *module = PyImport_AddModule("swig_runtime_data" SWIG_RUNTIME_VERSION);
#else
  static PyMethodDef swig_empty_runtime_method_table[] = { {NULL, NULL, 0, NULL} }; /* Sentinel */
  PyObject *module = Py_InitModule("swig_runtime_data" SWIG_RUNTIME_VERSION, swig_empty_runtime_method_table);
#endif
  PyObject *pointer = PyCapsule_New((void *) swig_module, SWIGPY_CAPSULE_NAME, SWIG_Python_DestroyModule);
  if (pointer && module) {
    PyModule_AddObject(module, "type_pointer_capsule" SWIG_TYPE_TABLE_NAME, pointer);
  } else {
    Py_XDECREF(pointer);
  }
}

/* The python cached type query */
SWIGRUNTIME PyObject *
SWIG_Python_TypeCache(void) {
  static PyObject *SWIG_STATIC_POINTER(cache) = PyDict_New();
  return cache;
}

SWIGRUNTIME swig_type_info *
SWIG_Python_TypeQuery(const char *type)
{
  PyObject *cache = SWIG_Python_TypeCache();
  PyObject *key = SWIG_Python_str_FromChar(type); 
  PyObject *obj = PyDict_GetItem(cache, key);
  swig_type_info *descriptor;
  if (obj) {
    descriptor = (swig_type_info *) PyCapsule_GetPointer(obj, NULL);
  } else {
    swig_module_info *swig_module = SWIG_GetModule(0);
    descriptor = SWIG_TypeQueryModule(swig_module, swig_module, type);
    if (descriptor) {
      obj = PyCapsule_New((void*) descriptor, NULL, NULL);
      PyDict_SetItem(cache, key, obj);
      Py_DECREF(obj);
    }
  }
  Py_DECREF(key);
  return descriptor;
}

/* 
   For backward compatibility only
*/
#define SWIG_POINTER_EXCEPTION  0
#define SWIG_arg_fail(arg)      SWIG_Python_ArgFail(arg)
#define SWIG_MustGetPtr(p, type, argnum, flags)  SWIG_Python_MustGetPtr(p, type, argnum, flags)

SWIGRUNTIME int
SWIG_Python_AddErrMesg(const char* mesg, int infront)
{  
  if (PyErr_Occurred()) {
    PyObject *type = 0;
    PyObject *value = 0;
    PyObject *traceback = 0;
    PyErr_Fetch(&type, &value, &traceback);
    if (value) {
      PyObject *old_str = PyObject_Str(value);
      const char *tmp = SWIG_Python_str_AsChar(old_str);
      const char *errmesg = tmp ? tmp : "Invalid error message";
      Py_XINCREF(type);
      PyErr_Clear();
      if (infront) {
	PyErr_Format(type, "%s %s", mesg, errmesg);
      } else {
	PyErr_Format(type, "%s %s", errmesg, mesg);
      }
      SWIG_Python_str_DelForPy3(tmp);
      Py_DECREF(old_str);
    }
    return 1;
  } else {
    return 0;
  }
}
  
SWIGRUNTIME int
SWIG_Python_ArgFail(int argnum)
{
  if (PyErr_Occurred()) {
    /* add information about failing argument */
    char mesg[256];
    PyOS_snprintf(mesg, sizeof(mesg), "argument number %d:", argnum);
    return SWIG_Python_AddErrMesg(mesg, 1);
  } else {
    return 0;
  }
}

SWIGRUNTIMEINLINE const char *
SwigPyObject_GetDesc(PyObject *self)
{
  SwigPyObject *v = (SwigPyObject *)self;
  swig_type_info *ty = v ? v->ty : 0;
  return ty ? ty->str : "";
}

SWIGRUNTIME void
SWIG_Python_TypeError(const char *type, PyObject *obj)
{
  if (type) {
#if defined(SWIG_COBJECT_TYPES)
    if (obj && SwigPyObject_Check(obj)) {
      const char *otype = (const char *) SwigPyObject_GetDesc(obj);
      if (otype) {
	PyErr_Format(PyExc_TypeError, "a '%s' is expected, 'SwigPyObject(%s)' is received",
		     type, otype);
	return;
      }
    } else 
#endif      
    {
      const char *otype = (obj ? obj->ob_type->tp_name : 0); 
      if (otype) {
	PyObject *str = PyObject_Str(obj);
	const char *cstr = str ? SWIG_Python_str_AsChar(str) : 0;
	if (cstr) {
	  PyErr_Format(PyExc_TypeError, "a '%s' is expected, '%s(%s)' is received",
		       type, otype, cstr);
          SWIG_Python_str_DelForPy3(cstr);
	} else {
	  PyErr_Format(PyExc_TypeError, "a '%s' is expected, '%s' is received",
		       type, otype);
	}
	Py_XDECREF(str);
	return;
      }
    }   
    PyErr_Format(PyExc_TypeError, "a '%s' is expected", type);
  } else {
    PyErr_Format(PyExc_TypeError, "unexpected type is received");
  }
}


/* Convert a pointer value, signal an exception on a type mismatch */
SWIGRUNTIME void *
SWIG_Python_MustGetPtr(PyObject *obj, swig_type_info *ty, int SWIGUNUSEDPARM(argnum), int flags) {
  void *result;
  if (SWIG_Python_ConvertPtr(obj, &result, ty, flags) == -1) {
    PyErr_Clear();
#if SWIG_POINTER_EXCEPTION
    if (flags) {
      SWIG_Python_TypeError(SWIG_TypePrettyName(ty), obj);
      SWIG_Python_ArgFail(argnum);
    }
#endif
  }
  return result;
}

#ifdef SWIGPYTHON_BUILTIN
SWIGRUNTIME int
SWIG_Python_NonDynamicSetAttr(PyObject *obj, PyObject *name, PyObject *value) {
  PyTypeObject *tp = obj->ob_type;
  PyObject *descr;
  PyObject *encoded_name;
  descrsetfunc f;
  int res = -1;

# ifdef Py_USING_UNICODE
  if (PyString_Check(name)) {
    name = PyUnicode_Decode(PyString_AsString(name), PyString_Size(name), NULL, NULL);
    if (!name)
      return -1;
  } else if (!PyUnicode_Check(name))
# else
  if (!PyString_Check(name))
# endif
  {
    PyErr_Format(PyExc_TypeError, "attribute name must be string, not '%.200s'", name->ob_type->tp_name);
    return -1;
  } else {
    Py_INCREF(name);
  }

  if (!tp->tp_dict) {
    if (PyType_Ready(tp) < 0)
      goto done;
  }

  descr = _PyType_Lookup(tp, name);
  f = NULL;
  if (descr != NULL)
    f = descr->ob_type->tp_descr_set;
  if (!f) {
    if (PyString_Check(name)) {
      encoded_name = name;
      Py_INCREF(name);
    } else {
      encoded_name = PyUnicode_AsUTF8String(name);
      if (!encoded_name)
        return -1;
    }
    PyErr_Format(PyExc_AttributeError, "'%.100s' object has no attribute '%.200s'", tp->tp_name, PyString_AsString(encoded_name));
    Py_DECREF(encoded_name);
  } else {
    res = f(descr, obj, value);
  }
  
  done:
  Py_DECREF(name);
  return res;
}
#endif


#ifdef __cplusplus
}
#endif



#define SWIG_exception_fail(code, msg) do { SWIG_Error(code, msg); SWIG_fail; } while(0) 

#define SWIG_contract_assert(expr, msg) if (!(expr)) { SWIG_Error(SWIG_RuntimeError, msg); SWIG_fail; } else 



#ifdef __cplusplus
extern "C" {
#endif

/* Method creation and docstring support functions */

SWIGINTERN PyMethodDef *SWIG_PythonGetProxyDoc(const char *name);
SWIGINTERN PyObject *SWIG_PyInstanceMethod_New(PyObject *SWIGUNUSEDPARM(self), PyObject *func);
SWIGINTERN PyObject *SWIG_PyStaticMethod_New(PyObject *SWIGUNUSEDPARM(self), PyObject *func);

#ifdef __cplusplus
}
#endif


  #define SWIG_exception(code, msg) do { SWIG_Error(code, msg); SWIG_fail;; } while(0) 


/* -------- TYPES TABLE (BEGIN) -------- */

#define SWIGTYPE_p_allocator_type swig_types[0]
#define SWIGTYPE_p_char swig_types[1]
#define SWIGTYPE_p_difference_type swig_types[2]
#define SWIGTYPE_p_p_PyObject swig_types[3]
#define SWIGTYPE_p_size_type swig_types[4]
#define SWIGTYPE_p_std__allocatorT_double_t swig_types[5]
#define SWIGTYPE_p_std__allocatorT_std__string_t swig_types[6]
#define SWIGTYPE_p_std__invalid_argument swig_types[7]
#define SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t swig_types[8]
#define SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t swig_types[9]
#define SWIGTYPE_p_std__vectorT_std__vectorT_double_std__allocatorT_double_t_t_std__allocatorT_std__vectorT_double_std__allocatorT_double_t_t_t_t swig_types[10]
#define SWIGTYPE_p_swig__SwigPyIterator swig_types[11]
#define SWIGTYPE_p_value_type swig_types[12]
static swig_type_info *swig_types[14];
static swig_module_info swig_module = {swig_types, 13, 0, 0, 0, 0};
#define SWIG_TypeQuery(name) SWIG_TypeQueryModule(&swig_module, &swig_module, name)
#define SWIG_MangledTypeQuery(name) SWIG_MangledTypeQueryModule(&swig_module, &swig_module, name)

/* -------- TYPES TABLE (END) -------- */

#ifdef SWIG_TypeQuery
# undef SWIG_TypeQuery
#endif
#define SWIG_TypeQuery SWIG_Python_TypeQuery

/*-----------------------------------------------
              @(target):= _fractional_bs_model.so
  ------------------------------------------------*/
#if PY_VERSION_HEX >= 0x03000000
#  define SWIG_init    PyInit__fractional_bs_model

#else
#  define SWIG_init    init_fractional_bs_model

#endif
#define SWIG_name    "_fractional_bs_model"

#define SWIGVERSION 0x040002 
#define SWIG_VERSION SWIGVERSION


#define SWIG_as_voidptr(a) const_cast< void * >(static_cast< const void * >(a)) 
#define SWIG_as_voidptrptr(a) ((void)SWIG_as_voidptr(*a),reinterpret_cast< void** >(a)) 


#include <stdexcept>


namespace swig {
  class SwigPtr_PyObject {
  protected:
    PyObject *_obj;

  public:
    SwigPtr_PyObject() :_obj(0)
    {
    }

    SwigPtr_PyObject(const SwigPtr_PyObject& item) : _obj(item._obj)
    {
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
      Py_XINCREF(_obj);      
      SWIG_PYTHON_THREAD_END_BLOCK;
    }
    
    SwigPtr_PyObject(PyObject *obj, bool initial_ref = true) :_obj(obj)
    {
      if (initial_ref) {
        SWIG_PYTHON_THREAD_BEGIN_BLOCK;
        Py_XINCREF(_obj);
        SWIG_PYTHON_THREAD_END_BLOCK;
      }
    }
    
    SwigPtr_PyObject & operator=(const SwigPtr_PyObject& item) 
    {
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
      Py_XINCREF(item._obj);
      Py_XDECREF(_obj);
      _obj = item._obj;
      SWIG_PYTHON_THREAD_END_BLOCK;
      return *this;      
    }
    
    ~SwigPtr_PyObject() 
    {
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
      Py_XDECREF(_obj);
      SWIG_PYTHON_THREAD_END_BLOCK;
    }
    
    operator PyObject *() const
    {
      return _obj;
    }

    PyObject *operator->() const
    {
      return _obj;
    }
  };
}


namespace swig {
  struct SwigVar_PyObject : SwigPtr_PyObject {
    SwigVar_PyObject(PyObject* obj = 0) : SwigPtr_PyObject(obj, false) { }
    
    SwigVar_PyObject & operator = (PyObject* obj)
    {
      Py_XDECREF(_obj);
      _obj = obj;
      return *this;      
    }
  };
}


// Include necessary headers
#include "fractional_bs_model.h"
#include <vector>
#include <string>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>  // NumPy C API header

// Function to initialize NumPy's C API
int init_numpy() {
    if (_import_array() < 0) {
        PyErr_Print();  // Print the error to help diagnose the issue
        return -1;  // Indicate failure
    }
    return 0;  // Indicate success
}


#include <iostream>

#if PY_VERSION_HEX >= 0x03020000
# define SWIGPY_SLICE_ARG(obj) ((PyObject*) (obj))
#else
# define SWIGPY_SLICE_ARG(obj) ((PySliceObject*) (obj))
#endif


#include <typeinfo>
#include <stdexcept>


#if defined(__GNUC__)
#  if __GNUC__ == 2 && __GNUC_MINOR <= 96
#     define SWIG_STD_NOMODERN_STL
#  endif
#endif


#include <string>


#include <stddef.h>


namespace swig {
  struct stop_iteration {
  };

  struct SwigPyIterator {
  private:
    SwigPtr_PyObject _seq;

  protected:
    SwigPyIterator(PyObject *seq) : _seq(seq)
    {
    }
      
  public:
    virtual ~SwigPyIterator() {}

    // Access iterator method, required by Python
    virtual PyObject *value() const = 0;

    // Forward iterator method, required by Python
    virtual SwigPyIterator *incr(size_t n = 1) = 0;
    
    // Backward iterator method, very common in C++, but not required in Python
    virtual SwigPyIterator *decr(size_t /*n*/ = 1)
    {
      throw stop_iteration();
    }

    // Random access iterator methods, but not required in Python
    virtual ptrdiff_t distance(const SwigPyIterator &/*x*/) const
    {
      throw std::invalid_argument("operation not supported");
    }

    virtual bool equal (const SwigPyIterator &/*x*/) const
    {
      throw std::invalid_argument("operation not supported");
    }
    
    // C++ common/needed methods
    virtual SwigPyIterator *copy() const = 0;

    PyObject *next()     
    {
      SWIG_PYTHON_THREAD_BEGIN_BLOCK; // disable threads       
      PyObject *obj = value();
      incr();       
      SWIG_PYTHON_THREAD_END_BLOCK; // re-enable threads
      return obj;     
    }

    /* Make an alias for Python 3.x */
    PyObject *__next__()
    {
      return next();
    }

    PyObject *previous()
    {
      SWIG_PYTHON_THREAD_BEGIN_BLOCK; // disable threads       
      decr();
      PyObject *obj = value();
      SWIG_PYTHON_THREAD_END_BLOCK; // re-enable threads       
      return obj;
    }

    SwigPyIterator *advance(ptrdiff_t n)
    {
      return  (n > 0) ?  incr(n) : decr(-n);
    }
      
    bool operator == (const SwigPyIterator& x)  const
    {
      return equal(x);
    }
      
    bool operator != (const SwigPyIterator& x) const
    {
      return ! operator==(x);
    }
      
    SwigPyIterator& operator += (ptrdiff_t n)
    {
      return *advance(n);
    }

    SwigPyIterator& operator -= (ptrdiff_t n)
    {
      return *advance(-n);
    }
      
    SwigPyIterator* operator + (ptrdiff_t n) const
    {
      return copy()->advance(n);
    }

    SwigPyIterator* operator - (ptrdiff_t n) const
    {
      return copy()->advance(-n);
    }
      
    ptrdiff_t operator - (const SwigPyIterator& x) const
    {
      return x.distance(*this);
    }
      
    static swig_type_info* descriptor() {
      static int init = 0;
      static swig_type_info* desc = 0;
      if (!init) {
	desc = SWIG_TypeQuery("swig::SwigPyIterator *");
	init = 1;
      }	
      return desc;
    }    
  };

#if defined(SWIGPYTHON_BUILTIN)
  inline PyObject* make_output_iterator_builtin (PyObject *pyself)
  {
    Py_INCREF(pyself);
    return pyself;
  }
#endif
}


SWIGINTERN int
SWIG_AsVal_double (PyObject *obj, double *val)
{
  int res = SWIG_TypeError;
  if (PyFloat_Check(obj)) {
    if (val) *val = PyFloat_AsDouble(obj);
    return SWIG_OK;
#if PY_VERSION_HEX < 0x03000000
  } else if (PyInt_Check(obj)) {
    if (val) *val = (double) PyInt_AsLong(obj);
    return SWIG_OK;
#endif
  } else if (PyLong_Check(obj)) {
    double v = PyLong_AsDouble(obj);
    if (!PyErr_Occurred()) {
      if (val) *val = v;
      return SWIG_OK;
    } else {
      PyErr_Clear();
    }
  }
#ifdef SWIG_PYTHON_CAST_MODE
  {
    int dispatch = 0;
    double d = PyFloat_AsDouble(obj);
    if (!PyErr_Occurred()) {
      if (val) *val = d;
      return SWIG_AddCast(SWIG_OK);
    } else {
      PyErr_Clear();
    }
    if (!dispatch) {
      long v = PyLong_AsLong(obj);
      if (!PyErr_Occurred()) {
	if (val) *val = v;
	return SWIG_AddCast(SWIG_AddCast(SWIG_OK));
      } else {
	PyErr_Clear();
      }
    }
  }
#endif
  return res;
}


#include <float.h>


#include <math.h>


SWIGINTERNINLINE int
SWIG_CanCastAsInteger(double *d, double min, double max) {
  double x = *d;
  if ((min <= x && x <= max)) {
   double fx = floor(x);
   double cx = ceil(x);
   double rd =  ((x - fx) < 0.5) ? fx : cx; /* simple rint */
   if ((errno == EDOM) || (errno == ERANGE)) {
     errno = 0;
   } else {
     double summ, reps, diff;
     if (rd < x) {
       diff = x - rd;
     } else if (rd > x) {
       diff = rd - x;
     } else {
       return 1;
     }
     summ = rd + x;
     reps = diff/summ;
     if (reps < 8*DBL_EPSILON) {
       *d = rd;
       return 1;
     }
   }
  }
  return 0;
}


SWIGINTERN int
SWIG_AsVal_unsigned_SS_long (PyObject *obj, unsigned long *val) 
{
#if PY_VERSION_HEX < 0x03000000
  if (PyInt_Check(obj)) {
    long v = PyInt_AsLong(obj);
    if (v >= 0) {
      if (val) *val = v;
      return SWIG_OK;
    } else {
      return SWIG_OverflowError;
    }
  } else
#endif
  if (PyLong_Check(obj)) {
    unsigned long v = PyLong_AsUnsignedLong(obj);
    if (!PyErr_Occurred()) {
      if (val) *val = v;
      return SWIG_OK;
    } else {
      PyErr_Clear();
      return SWIG_OverflowError;
    }
  }
#ifdef SWIG_PYTHON_CAST_MODE
  {
    int dispatch = 0;
    unsigned long v = PyLong_AsUnsignedLong(obj);
    if (!PyErr_Occurred()) {
      if (val) *val = v;
      return SWIG_AddCast(SWIG_OK);
    } else {
      PyErr_Clear();
    }
    if (!dispatch) {
      double d;
      int res = SWIG_AddCast(SWIG_AsVal_double (obj,&d));
      if (SWIG_IsOK(res) && SWIG_CanCastAsInteger(&d, 0, ULONG_MAX)) {
	if (val) *val = (unsigned long)(d);
	return res;
      }
    }
  }
#endif
  return SWIG_TypeError;
}


#include <limits.h>
#if !defined(SWIG_NO_LLONG_MAX)
# if !defined(LLONG_MAX) && defined(__GNUC__) && defined (__LONG_LONG_MAX__)
#   define LLONG_MAX __LONG_LONG_MAX__
#   define LLONG_MIN (-LLONG_MAX - 1LL)
#   define ULLONG_MAX (LLONG_MAX * 2ULL + 1ULL)
# endif
#endif


#if defined(LLONG_MAX) && !defined(SWIG_LONG_LONG_AVAILABLE)
#  define SWIG_LONG_LONG_AVAILABLE
#endif


#ifdef SWIG_LONG_LONG_AVAILABLE
SWIGINTERN int
SWIG_AsVal_unsigned_SS_long_SS_long (PyObject *obj, unsigned long long *val)
{
  int res = SWIG_TypeError;
  if (PyLong_Check(obj)) {
    unsigned long long v = PyLong_AsUnsignedLongLong(obj);
    if (!PyErr_Occurred()) {
      if (val) *val = v;
      return SWIG_OK;
    } else {
      PyErr_Clear();
      res = SWIG_OverflowError;
    }
  } else {
    unsigned long v;
    res = SWIG_AsVal_unsigned_SS_long (obj,&v);
    if (SWIG_IsOK(res)) {
      if (val) *val = v;
      return res;
    }
  }
#ifdef SWIG_PYTHON_CAST_MODE
  {
    const double mant_max = 1LL << DBL_MANT_DIG;
    double d;
    res = SWIG_AsVal_double (obj,&d);
    if (SWIG_IsOK(res) && !SWIG_CanCastAsInteger(&d, 0, mant_max))
      return SWIG_OverflowError;
    if (SWIG_IsOK(res) && SWIG_CanCastAsInteger(&d, 0, mant_max)) {
      if (val) *val = (unsigned long long)(d);
      return SWIG_AddCast(res);
    }
    res = SWIG_TypeError;
  }
#endif
  return res;
}
#endif


SWIGINTERNINLINE int
SWIG_AsVal_size_t (PyObject * obj, size_t *val)
{
  int res = SWIG_TypeError;
#ifdef SWIG_LONG_LONG_AVAILABLE
  if (sizeof(size_t) <= sizeof(unsigned long)) {
#endif
    unsigned long v;
    res = SWIG_AsVal_unsigned_SS_long (obj, val ? &v : 0);
    if (SWIG_IsOK(res) && val) *val = static_cast< size_t >(v);
#ifdef SWIG_LONG_LONG_AVAILABLE
  } else if (sizeof(size_t) <= sizeof(unsigned long long)) {
    unsigned long long v;
    res = SWIG_AsVal_unsigned_SS_long_SS_long (obj, val ? &v : 0);
    if (SWIG_IsOK(res) && val) *val = static_cast< size_t >(v);
  }
#endif
  return res;
}


  #define SWIG_From_long   PyInt_FromLong 


#ifdef SWIG_LONG_LONG_AVAILABLE
SWIGINTERNINLINE PyObject* 
SWIG_From_long_SS_long  (long long value)
{
  return ((value < LONG_MIN) || (value > LONG_MAX)) ?
    PyLong_FromLongLong(value) : PyInt_FromLong(static_cast< long >(value));
}
#endif


SWIGINTERNINLINE PyObject *
SWIG_From_ptrdiff_t  (ptrdiff_t value)
{    
#ifdef SWIG_LONG_LONG_AVAILABLE
  if (sizeof(ptrdiff_t) <= sizeof(long)) {
#endif
    return SWIG_From_long  (static_cast< long >(value));
#ifdef SWIG_LONG_LONG_AVAILABLE
  } else {
    /* assume sizeof(ptrdiff_t) <= sizeof(long long) */
    return SWIG_From_long_SS_long  (static_cast< long long >(value));
  }
#endif
}


SWIGINTERNINLINE PyObject*
  SWIG_From_bool  (bool value)
{
  return PyBool_FromLong(value ? 1 : 0);
}


SWIGINTERN int
SWIG_AsVal_long (PyObject *obj, long* val)
{
#if PY_VERSION_HEX < 0x03000000
  if (PyInt_Check(obj)) {
    if (val) *val = PyInt_AsLong(obj);
    return SWIG_OK;
  } else
#endif
  if (PyLong_Check(obj)) {
    long v = PyLong_AsLong(obj);
    if (!PyErr_Occurred()) {
      if (val) *val = v;
      return SWIG_OK;
    } else {
      PyErr_Clear();
      return SWIG_OverflowError;
    }
  }
#ifdef SWIG_PYTHON_CAST_MODE
  {
    int dispatch = 0;
    long v = PyInt_AsLong(obj);
    if (!PyErr_Occurred()) {
      if (val) *val = v;
      return SWIG_AddCast(SWIG_OK);
    } else {
      PyErr_Clear();
    }
    if (!dispatch) {
      double d;
      int res = SWIG_AddCast(SWIG_AsVal_double (obj,&d));
      if (SWIG_IsOK(res) && SWIG_CanCastAsInteger(&d, LONG_MIN, LONG_MAX)) {
	if (val) *val = (long)(d);
	return res;
      }
    }
  }
#endif
  return SWIG_TypeError;
}


#ifdef SWIG_LONG_LONG_AVAILABLE
SWIGINTERN int
SWIG_AsVal_long_SS_long (PyObject *obj, long long *val)
{
  int res = SWIG_TypeError;
  if (PyLong_Check(obj)) {
    long long v = PyLong_AsLongLong(obj);
    if (!PyErr_Occurred()) {
      if (val) *val = v;
      return SWIG_OK;
    } else {
      PyErr_Clear();
      res = SWIG_OverflowError;
    }
  } else {
    long v;
    res = SWIG_AsVal_long (obj,&v);
    if (SWIG_IsOK(res)) {
      if (val) *val = v;
      return res;
    }
  }
#ifdef SWIG_PYTHON_CAST_MODE
  {
    const double mant_max = 1LL << DBL_MANT_DIG;
    const double mant_min = -mant_max;
    double d;
    res = SWIG_AsVal_double (obj,&d);
    if (SWIG_IsOK(res) && !SWIG_CanCastAsInteger(&d, mant_min, mant_max))
      return SWIG_OverflowError;
    if (SWIG_IsOK(res) && SWIG_CanCastAsInteger(&d, mant_min, mant_max)) {
      if (val) *val = (long long)(d);
      return SWIG_AddCast(res);
    }
    res = SWIG_TypeError;
  }
#endif
  return res;
}
#endif


SWIGINTERNINLINE int
SWIG_AsVal_ptrdiff_t (PyObject * obj, ptrdiff_t *val)
{
  int res = SWIG_TypeError;
#ifdef SWIG_LONG_LONG_AVAILABLE
  if (sizeof(ptrdiff_t) <= sizeof(long)) {
#endif
    long v;
    res = SWIG_AsVal_long (obj, val ? &v : 0);
    if (SWIG_IsOK(res) && val) *val = static_cast< ptrdiff_t >(v);
#ifdef SWIG_LONG_LONG_AVAILABLE
  } else if (sizeof(ptrdiff_t) <= sizeof(long long)) {
    long long v;
    res = SWIG_AsVal_long_SS_long (obj, val ? &v : 0);
    if (SWIG_IsOK(res) && val) *val = static_cast< ptrdiff_t >(v);
  }
#endif
  return res;
}


#include <algorithm>


#include <vector>


#ifndef SWIG_FILE_WITH_INIT
#define NO_IMPORT_ARRAY
#endif
#include "stdio.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>


#include <complex> 


namespace swig {
  template <class Type>
  struct noconst_traits {
    typedef Type noconst_type;
  };

  template <class Type>
  struct noconst_traits<const Type> {
    typedef Type noconst_type;
  };

  /*
    type categories
  */
  struct pointer_category { };
  struct value_category { };

  /*
    General traits that provides type_name and type_info
  */
  template <class Type> struct traits { };

  template <class Type>
  inline const char* type_name() {
    return traits<typename noconst_traits<Type >::noconst_type >::type_name();
  }

  template <class Type> struct traits_info {
    static swig_type_info *type_query(std::string name) {
      name += " *";
      return SWIG_TypeQuery(name.c_str());
    }
    static swig_type_info *type_info() {
      static swig_type_info *info = type_query(type_name<Type>());
      return info;
    }
  };

  /*
    Partial specialization for pointers (traits_info)
  */
  template <class Type> struct traits_info<Type *> {
    static swig_type_info *type_query(std::string name) {
      name += " *";
      return SWIG_TypeQuery(name.c_str());
    }
    static swig_type_info *type_info() {
      static swig_type_info *info = type_query(type_name<Type>());
      return info;
    }
  };

  template <class Type>
  inline swig_type_info *type_info() {
    return traits_info<Type>::type_info();
  }

  /*
    Partial specialization for pointers (traits)
  */
  template <class Type> struct traits <Type *> {
    typedef pointer_category category;
    static std::string make_ptr_name(const char* name) {
      std::string ptrname = name;
      ptrname += " *";
      return ptrname;
    }
    static const char* type_name() {
      static std::string name = make_ptr_name(swig::type_name<Type>());
      return name.c_str();
    }
  };

  template <class Type, class Category>
  struct traits_as { };

  template <class Type, class Category>
  struct traits_check { };

}


namespace swig {  
  /*
    Traits that provides the from method
  */
  template <class Type> struct traits_from_ptr {
    static PyObject *from(Type *val, int owner = 0) {
      return SWIG_InternalNewPointerObj(val, type_info<Type>(), owner);
    }
  };

  template <class Type> struct traits_from {
    static PyObject *from(const Type& val) {
      return traits_from_ptr<Type>::from(new Type(val), 1);
    }
  };

  template <class Type> struct traits_from<Type *> {
    static PyObject *from(Type* val) {
      return traits_from_ptr<Type>::from(val, 0);
    }
  };

  template <class Type> struct traits_from<const Type *> {
    static PyObject *from(const Type* val) {
      return traits_from_ptr<Type>::from(const_cast<Type*>(val), 0);
    }
  };


  template <class Type>
  inline PyObject *from(const Type& val) {
    return traits_from<Type>::from(val);
  }

  template <class Type>
  inline PyObject *from_ptr(Type* val, int owner) {
    return traits_from_ptr<Type>::from(val, owner);
  }

  /*
    Traits that provides the asval/as/check method
  */
  template <class Type>
  struct traits_asptr {   
    static int asptr(PyObject *obj, Type **val) {
      int res = SWIG_ERROR;
      swig_type_info *descriptor = type_info<Type>();
      if (val) {
        Type *p = 0;
        int newmem = 0;
        res = descriptor ? SWIG_ConvertPtrAndOwn(obj, (void **)&p, descriptor, 0, &newmem) : SWIG_ERROR;
        if (SWIG_IsOK(res)) {
          if (newmem & SWIG_CAST_NEW_MEMORY) {
            res |= SWIG_NEWOBJMASK;
          }
          *val = p;
        }
      } else {
        res = descriptor ? SWIG_ConvertPtr(obj, 0, descriptor, 0) : SWIG_ERROR;
      }
      return res;
    }
  }; 

  template <class Type>
  inline int asptr(PyObject *obj, Type **vptr) {
    return traits_asptr<Type>::asptr(obj, vptr);
  }

  template <class Type> 
  struct traits_asval {
    static int asval(PyObject *obj, Type *val) {
      if (val) {
	Type *p = 0;
	int res = traits_asptr<Type>::asptr(obj, &p);
	if (!SWIG_IsOK(res)) return res;	
	if (p) {
	  typedef typename noconst_traits<Type>::noconst_type noconst_type;
	  *(const_cast<noconst_type*>(val)) = *p;
	  if (SWIG_IsNewObj(res)){
	    delete p;
	    res = SWIG_DelNewMask(res);
	  }
	  return res;
	} else {
	  return SWIG_ERROR;
	}
      } else {
	return traits_asptr<Type>::asptr(obj, (Type **)(0));
      }
    }
  };

  template <class Type> struct traits_asval<Type*> {
    static int asval(PyObject *obj, Type **val) {
      if (val) {
        typedef typename noconst_traits<Type>::noconst_type noconst_type;
        noconst_type *p = 0;
        int res = traits_asptr<noconst_type>::asptr(obj,  &p);
        if (SWIG_IsOK(res)) {
          *(const_cast<noconst_type**>(val)) = p;
	}
	return res;
      } else {
	return traits_asptr<Type>::asptr(obj, (Type **)(0));
      }
    }
  };
  
  template <class Type>
  inline int asval(PyObject *obj, Type *val) {
    return traits_asval<Type>::asval(obj, val);
  }

  template <class Type> 
  struct traits_as<Type, value_category> {
    static Type as(PyObject *obj) {
      Type v;
      int res = asval(obj, &v);
      if (!obj || !SWIG_IsOK(res)) {
	if (!PyErr_Occurred()) {
	  ::SWIG_Error(SWIG_TypeError,  swig::type_name<Type>());
	}
	throw std::invalid_argument("bad type");
      }
      return v;
    }
  };

  template <class Type> 
  struct traits_as<Type, pointer_category> {
    static Type as(PyObject *obj) {
      Type *v = 0;      
      int res = (obj ? traits_asptr<Type>::asptr(obj, &v) : SWIG_ERROR);
      if (SWIG_IsOK(res) && v) {
	if (SWIG_IsNewObj(res)) {
	  Type r(*v);
	  delete v;
	  return r;
	} else {
	  return *v;
	}
      } else {
	if (!PyErr_Occurred()) {
	  SWIG_Error(SWIG_TypeError,  swig::type_name<Type>());
	}
	throw std::invalid_argument("bad type");
      }
    }
  };

  template <class Type> 
  struct traits_as<Type*, pointer_category> {
    static Type* as(PyObject *obj) {
      Type *v = 0;      
      int res = (obj ? traits_asptr<Type>::asptr(obj, &v) : SWIG_ERROR);
      if (SWIG_IsOK(res)) {
	return v;
      } else {
	if (!PyErr_Occurred()) {
	  SWIG_Error(SWIG_TypeError,  swig::type_name<Type>());
	}
	throw std::invalid_argument("bad type");
      }
    }
  };
    
  template <class Type>
  inline Type as(PyObject *obj) {
    return traits_as<Type, typename traits<Type>::category>::as(obj);
  }

  template <class Type> 
  struct traits_check<Type, value_category> {
    static bool check(PyObject *obj) {
      int res = obj ? asval(obj, (Type *)(0)) : SWIG_ERROR;
      return SWIG_IsOK(res) ? true : false;
    }
  };

  template <class Type> 
  struct traits_check<Type, pointer_category> {
    static bool check(PyObject *obj) {
      int res = obj ? asptr(obj, (Type **)(0)) : SWIG_ERROR;
      return SWIG_IsOK(res) ? true : false;
    }
  };

  template <class Type>
  inline bool check(PyObject *obj) {
    return traits_check<Type, typename traits<Type>::category>::check(obj);
  }
}


#include <functional>

namespace std {
  template <>
  struct less <PyObject *>
  {
    bool
    operator()(PyObject * v, PyObject *w) const
    { 
      bool res;
      SWIG_PYTHON_THREAD_BEGIN_BLOCK;
      res = PyObject_RichCompareBool(v, w, Py_LT) ? true : false;
      /* This may fall into a case of inconsistent
               eg. ObjA > ObjX > ObjB
               but ObjA < ObjB
      */
      if( PyErr_Occurred() && PyErr_ExceptionMatches(PyExc_TypeError) )
      {
        /* Objects can't be compared, this mostly occurred in Python 3.0 */
        /* Compare their ptr directly for a workaround */
        res = (v < w);
        PyErr_Clear();
      }
      SWIG_PYTHON_THREAD_END_BLOCK;
      return res;
    }
  };

  template <>
  struct less <swig::SwigPtr_PyObject>
  {
    bool
    operator()(const swig::SwigPtr_PyObject& v, const swig::SwigPtr_PyObject& w) const
    {
      return std::less<PyObject *>()(v, w);
    }
  };

  template <>
  struct less <swig::SwigVar_PyObject>
  {
    bool
    operator()(const swig::SwigVar_PyObject& v, const swig::SwigVar_PyObject& w) const
    {
      return std::less<PyObject *>()(v, w);
    }
  };

}

namespace swig {
  template <> struct traits<PyObject *> {
    typedef value_category category;
    static const char* type_name() { return "PyObject *"; }
  };  

  template <>  struct traits_asval<PyObject * > {   
    typedef PyObject * value_type;
    static int asval(PyObject *obj, value_type *val) {
      if (val) *val = obj;
      return SWIG_OK;
    }
  };

  template <> 
  struct traits_check<PyObject *, value_category> {
    static bool check(PyObject *) {
      return true;
    }
  };

  template <>  struct traits_from<PyObject *> {
    typedef PyObject * value_type;
    static PyObject *from(const value_type& val) {
      Py_XINCREF(val);
      return val;
    }
  };
  
}

namespace swig {
  template <class Difference>
  inline size_t
  check_index(Difference i, size_t size, bool insert = false) {
    if ( i < 0 ) {
      if ((size_t) (-i) <= size)
	return (size_t) (i + size);
    } else if ( (size_t) i < size ) {
      return (size_t) i;
    } else if (insert && ((size_t) i == size)) {
      return size;
    }
    throw std::out_of_range("index out of range");
  }

  template <class Difference>
  void
  slice_adjust(Difference i, Difference j, Py_ssize_t step, size_t size, Difference &ii, Difference &jj, bool insert = false) {
    if (step == 0) {
      throw std::invalid_argument("slice step cannot be zero");
    } else if (step > 0) {
      // Required range: 0 <= i < size, 0 <= j < size, i <= j
      if (i < 0) {
        ii = 0;
      } else if (i < (Difference)size) {
        ii = i;
      } else if (insert && (i >= (Difference)size)) {
        ii = (Difference)size;
      }
      if (j < 0) {
        jj = 0;
      } else {
        jj = (j < (Difference)size) ? j : (Difference)size;
      }
      if (jj < ii)
        jj = ii;
    } else {
      // Required range: -1 <= i < size-1, -1 <= j < size-1, i >= j
      if (i < -1) {
        ii = -1;
      } else if (i < (Difference) size) {
        ii = i;
      } else if (i >= (Difference)(size-1)) {
        ii = (Difference)(size-1);
      }
      if (j < -1) {
        jj = -1;
      } else {
        jj = (j < (Difference)size ) ? j : (Difference)(size-1);
      }
      if (ii < jj)
        ii = jj;
    }
  }

  template <class Sequence, class Difference>
  inline typename Sequence::iterator
  getpos(Sequence* self, Difference i)  {
    typename Sequence::iterator pos = self->begin();
    std::advance(pos, check_index(i,self->size()));
    return pos;
  }

  template <class Sequence, class Difference>
  inline typename Sequence::const_iterator
  cgetpos(const Sequence* self, Difference i)  {
    typename Sequence::const_iterator pos = self->begin();
    std::advance(pos, check_index(i,self->size()));
    return pos;
  }

  template <class Sequence>
  inline void
  erase(Sequence* seq, const typename Sequence::iterator& position) {
    seq->erase(position);
  }

  template <class Sequence>
  struct traits_reserve {
    static void reserve(Sequence & /*seq*/, typename Sequence::size_type /*n*/) {
      // This should be specialized for types that support reserve
    }
  };

  template <class Sequence, class Difference>
  inline Sequence*
  getslice(const Sequence* self, Difference i, Difference j, Py_ssize_t step) {
    typename Sequence::size_type size = self->size();
    Difference ii = 0;
    Difference jj = 0;
    swig::slice_adjust(i, j, step, size, ii, jj);

    if (step > 0) {
      typename Sequence::const_iterator sb = self->begin();
      typename Sequence::const_iterator se = self->begin();
      std::advance(sb,ii);
      std::advance(se,jj);
      if (step == 1) {
        return new Sequence(sb, se);
      } else {
        Sequence *sequence = new Sequence();
        swig::traits_reserve<Sequence>::reserve(*sequence, (jj - ii + step - 1) / step);
        typename Sequence::const_iterator it = sb;
        while (it!=se) {
          sequence->push_back(*it);
          for (Py_ssize_t c=0; c<step && it!=se; ++c)
            it++;
        }
        return sequence;
      } 
    } else {
      Sequence *sequence = new Sequence();
      swig::traits_reserve<Sequence>::reserve(*sequence, (ii - jj - step - 1) / -step);
      typename Sequence::const_reverse_iterator sb = self->rbegin();
      typename Sequence::const_reverse_iterator se = self->rbegin();
      std::advance(sb,size-ii-1);
      std::advance(se,size-jj-1);
      typename Sequence::const_reverse_iterator it = sb;
      while (it!=se) {
        sequence->push_back(*it);
        for (Py_ssize_t c=0; c<-step && it!=se; ++c)
          it++;
      }
      return sequence;
    }
  }

  template <class Sequence, class Difference, class InputSeq>
  inline void
  setslice(Sequence* self, Difference i, Difference j, Py_ssize_t step, const InputSeq& is = InputSeq()) {
    typename Sequence::size_type size = self->size();
    Difference ii = 0;
    Difference jj = 0;
    swig::slice_adjust(i, j, step, size, ii, jj, true);
    if (step > 0) {
      if (step == 1) {
        size_t ssize = jj - ii;
        if (ssize <= is.size()) {
          // expanding/staying the same size
          swig::traits_reserve<Sequence>::reserve(*self, self->size() - ssize + is.size());
          typename Sequence::iterator sb = self->begin();
          typename InputSeq::const_iterator isit = is.begin();
          std::advance(sb,ii);
          std::advance(isit, jj - ii);
          self->insert(std::copy(is.begin(), isit, sb), isit, is.end());
        } else {
          // shrinking
          typename Sequence::iterator sb = self->begin();
          typename Sequence::iterator se = self->begin();
          std::advance(sb,ii);
          std::advance(se,jj);
          self->erase(sb,se);
          sb = self->begin();
          std::advance(sb,ii);
          self->insert(sb, is.begin(), is.end());
        }
      } else {
        size_t replacecount = (jj - ii + step - 1) / step;
        if (is.size() != replacecount) {
          char msg[1024];
          sprintf(msg, "attempt to assign sequence of size %lu to extended slice of size %lu", (unsigned long)is.size(), (unsigned long)replacecount);
          throw std::invalid_argument(msg);
        }
        typename Sequence::const_iterator isit = is.begin();
        typename Sequence::iterator it = self->begin();
        std::advance(it,ii);
        for (size_t rc=0; rc<replacecount && it != self->end(); ++rc) {
          *it++ = *isit++;
          for (Py_ssize_t c=0; c<(step-1) && it != self->end(); ++c)
            it++;
        }
      }
    } else {
      size_t replacecount = (ii - jj - step - 1) / -step;
      if (is.size() != replacecount) {
        char msg[1024];
        sprintf(msg, "attempt to assign sequence of size %lu to extended slice of size %lu", (unsigned long)is.size(), (unsigned long)replacecount);
        throw std::invalid_argument(msg);
      }
      typename Sequence::const_iterator isit = is.begin();
      typename Sequence::reverse_iterator it = self->rbegin();
      std::advance(it,size-ii-1);
      for (size_t rc=0; rc<replacecount && it != self->rend(); ++rc) {
        *it++ = *isit++;
        for (Py_ssize_t c=0; c<(-step-1) && it != self->rend(); ++c)
          it++;
      }
    }
  }

  template <class Sequence, class Difference>
  inline void
  delslice(Sequence* self, Difference i, Difference j, Py_ssize_t step) {
    typename Sequence::size_type size = self->size();
    Difference ii = 0;
    Difference jj = 0;
    swig::slice_adjust(i, j, step, size, ii, jj, true);
    if (step > 0) {
      typename Sequence::iterator sb = self->begin();
      std::advance(sb,ii);
      if (step == 1) {
        typename Sequence::iterator se = self->begin();
        std::advance(se,jj);
        self->erase(sb,se);
      } else {
        typename Sequence::iterator it = sb;
        size_t delcount = (jj - ii + step - 1) / step;
        while (delcount) {
          it = self->erase(it);
          for (Py_ssize_t c=0; c<(step-1) && it != self->end(); ++c)
            it++;
          delcount--;
        }
      }
    } else {
      typename Sequence::reverse_iterator sb = self->rbegin();
      std::advance(sb,size-ii-1);
      typename Sequence::reverse_iterator it = sb;
      size_t delcount = (ii - jj - step - 1) / -step;
      while (delcount) {
        it = typename Sequence::reverse_iterator(self->erase((++it).base()));
        for (Py_ssize_t c=0; c<(-step-1) && it != self->rend(); ++c)
          it++;
        delcount--;
      }
    }
  }
}


#if defined(__SUNPRO_CC) && defined(_RWSTD_VER)
#  if !defined(SWIG_NO_STD_NOITERATOR_TRAITS_STL)
#    define SWIG_STD_NOITERATOR_TRAITS_STL
#  endif
#endif

#if !defined(SWIG_STD_NOITERATOR_TRAITS_STL)
#include <iterator>
#else
namespace std {
  template <class Iterator>
  struct iterator_traits {
    typedef ptrdiff_t difference_type;
    typedef typename Iterator::value_type value_type;
  };

  template <class Iterator, class Category,class T, class Reference, class Pointer, class Distance>
  struct iterator_traits<__reverse_bi_iterator<Iterator,Category,T,Reference,Pointer,Distance> > {
    typedef Distance difference_type;
    typedef T value_type;
  };

  template <class T>
  struct iterator_traits<T*> {
    typedef T value_type;
    typedef ptrdiff_t difference_type;
  };

  template<typename _InputIterator>
  inline typename iterator_traits<_InputIterator>::difference_type
  distance(_InputIterator __first, _InputIterator __last)
  {
    typename iterator_traits<_InputIterator>::difference_type __n = 0;
    while (__first != __last) {
      ++__first; ++__n;
    }
    return __n;
  }
}
#endif


namespace swig {
  template<typename OutIterator>
  class SwigPyIterator_T :  public SwigPyIterator
  {
  public:
    typedef OutIterator out_iterator;
    typedef typename std::iterator_traits<out_iterator>::value_type value_type;    
    typedef SwigPyIterator_T<out_iterator> self_type;

    SwigPyIterator_T(out_iterator curr, PyObject *seq)
      : SwigPyIterator(seq), current(curr)
    {
    }

    const out_iterator& get_current() const
    {
      return current;
    }

    
    bool equal (const SwigPyIterator &iter) const
    {
      const self_type *iters = dynamic_cast<const self_type *>(&iter);
      if (iters) {
	return (current == iters->get_current());
      } else {
	throw std::invalid_argument("bad iterator type");
      }
    }
    
    ptrdiff_t distance(const SwigPyIterator &iter) const
    {
      const self_type *iters = dynamic_cast<const self_type *>(&iter);
      if (iters) {
	return std::distance(current, iters->get_current());
      } else {
	throw std::invalid_argument("bad iterator type");
      }
    }    
    
  protected:
    out_iterator current;
  };
  
  template <class ValueType>
  struct from_oper 
  {
    typedef const ValueType& argument_type;
    typedef PyObject *result_type;
    result_type operator()(argument_type v) const
    {
      return swig::from(v);
    }
  };

  template<typename OutIterator, 
	   typename ValueType = typename std::iterator_traits<OutIterator>::value_type,
	   typename FromOper = from_oper<ValueType> >
  class SwigPyForwardIteratorOpen_T :  public SwigPyIterator_T<OutIterator>
  {
  public:
    FromOper from;
    typedef OutIterator out_iterator;
    typedef ValueType value_type;
    typedef SwigPyIterator_T<out_iterator>  base;
    typedef SwigPyForwardIteratorOpen_T<OutIterator, ValueType, FromOper> self_type;
    
    SwigPyForwardIteratorOpen_T(out_iterator curr, PyObject *seq)
      : SwigPyIterator_T<OutIterator>(curr, seq)
    {
    }
    
    PyObject *value() const {
      return from(static_cast<const value_type&>(*(base::current)));
    }
    
    SwigPyIterator *copy() const
    {
      return new self_type(*this);
    }

    SwigPyIterator *incr(size_t n = 1)
    {
      while (n--) {
	++base::current;
      }
      return this;
    }

  };

  template<typename OutIterator, 
	   typename ValueType = typename std::iterator_traits<OutIterator>::value_type,
	   typename FromOper = from_oper<ValueType> >
  class SwigPyIteratorOpen_T :  public SwigPyForwardIteratorOpen_T<OutIterator, ValueType, FromOper>
  {
  public:
    FromOper from;
    typedef OutIterator out_iterator;
    typedef ValueType value_type;
    typedef SwigPyIterator_T<out_iterator>  base;
    typedef SwigPyIteratorOpen_T<OutIterator, ValueType, FromOper> self_type;
    
    SwigPyIteratorOpen_T(out_iterator curr, PyObject *seq)
      : SwigPyForwardIteratorOpen_T<OutIterator>(curr, seq)
    {
    }

    SwigPyIterator *decr(size_t n = 1)
    {
      while (n--) {
	--base::current;
      }
      return this;
    }
  };

  template<typename OutIterator, 
	   typename ValueType = typename std::iterator_traits<OutIterator>::value_type,
	   typename FromOper = from_oper<ValueType> >
  class SwigPyForwardIteratorClosed_T :  public SwigPyIterator_T<OutIterator>
  {
  public:
    FromOper from;
    typedef OutIterator out_iterator;
    typedef ValueType value_type;
    typedef SwigPyIterator_T<out_iterator>  base;    
    typedef SwigPyForwardIteratorClosed_T<OutIterator, ValueType, FromOper> self_type;
    
    SwigPyForwardIteratorClosed_T(out_iterator curr, out_iterator first, out_iterator last, PyObject *seq)
      : SwigPyIterator_T<OutIterator>(curr, seq), begin(first), end(last)
    {
    }
    
    PyObject *value() const {
      if (base::current == end) {
	throw stop_iteration();
      } else {
	return from(static_cast<const value_type&>(*(base::current)));
      }
    }
    
    SwigPyIterator *copy() const
    {
      return new self_type(*this);
    }

    SwigPyIterator *incr(size_t n = 1)
    {
      while (n--) {
	if (base::current == end) {
	  throw stop_iteration();
	} else {
	  ++base::current;
	}
      }
      return this;
    }

  protected:
    out_iterator begin;
    out_iterator end;
  };

  template<typename OutIterator, 
	   typename ValueType = typename std::iterator_traits<OutIterator>::value_type,
	   typename FromOper = from_oper<ValueType> >
  class SwigPyIteratorClosed_T :  public SwigPyForwardIteratorClosed_T<OutIterator,ValueType,FromOper>
  {
  public:
    FromOper from;
    typedef OutIterator out_iterator;
    typedef ValueType value_type;
    typedef SwigPyIterator_T<out_iterator>  base;
    typedef SwigPyForwardIteratorClosed_T<OutIterator, ValueType, FromOper> base0;
    typedef SwigPyIteratorClosed_T<OutIterator, ValueType, FromOper> self_type;
    
    SwigPyIteratorClosed_T(out_iterator curr, out_iterator first, out_iterator last, PyObject *seq)
      : SwigPyForwardIteratorClosed_T<OutIterator,ValueType,FromOper>(curr, first, last, seq)
    {
    }

    SwigPyIterator *decr(size_t n = 1)
    {
      while (n--) {
	if (base::current == base0::begin) {
	  throw stop_iteration();
	} else {
	  --base::current;
	}
      }
      return this;
    }
  };


  template<typename OutIter>
  inline SwigPyIterator*
  make_output_forward_iterator(const OutIter& current, const OutIter& begin,const OutIter& end, PyObject *seq = 0)
  {
    return new SwigPyForwardIteratorClosed_T<OutIter>(current, begin, end, seq);
  }

  template<typename OutIter>
  inline SwigPyIterator*
  make_output_iterator(const OutIter& current, const OutIter& begin,const OutIter& end, PyObject *seq = 0)
  {
    return new SwigPyIteratorClosed_T<OutIter>(current, begin, end, seq);
  }

  template<typename OutIter>
  inline SwigPyIterator*
  make_output_forward_iterator(const OutIter& current, PyObject *seq = 0)
  {
    return new SwigPyForwardIteratorOpen_T<OutIter>(current, seq);
  }

  template<typename OutIter>
  inline SwigPyIterator*
  make_output_iterator(const OutIter& current, PyObject *seq = 0)
  {
    return new SwigPyIteratorOpen_T<OutIter>(current, seq);
  }

}


namespace swig
{
  template <class T>
  struct SwigPySequence_Ref
  {
    SwigPySequence_Ref(PyObject* seq, Py_ssize_t index)
      : _seq(seq), _index(index)
    {
    }
    
    operator T () const
    {
      swig::SwigVar_PyObject item = PySequence_GetItem(_seq, _index);
      try {
	return swig::as<T>(item);
      } catch (const std::invalid_argument& e) {
	char msg[1024];
	sprintf(msg, "in sequence element %d ", (int)_index);
	if (!PyErr_Occurred()) {
	  ::SWIG_Error(SWIG_TypeError,  swig::type_name<T>());
	}
	SWIG_Python_AddErrorMsg(msg);
	SWIG_Python_AddErrorMsg(e.what());
	throw;
      }
    }

    SwigPySequence_Ref& operator=(const T& v)
    {
      PySequence_SetItem(_seq, _index, swig::from<T>(v));
      return *this;
    }

  private:
    PyObject* _seq;
    Py_ssize_t _index;
  };

  template <class T>
  struct SwigPySequence_ArrowProxy
  {
    SwigPySequence_ArrowProxy(const T& x): m_value(x) {}
    const T* operator->() const { return &m_value; }
    operator const T*() const { return &m_value; }
    T m_value;
  };

  template <class T, class Reference >
  struct SwigPySequence_InputIterator
  {
    typedef SwigPySequence_InputIterator<T, Reference > self;

    typedef std::random_access_iterator_tag iterator_category;
    typedef Reference reference;
    typedef T value_type;
    typedef T* pointer;
    typedef Py_ssize_t difference_type;

    SwigPySequence_InputIterator()
    {
    }

    SwigPySequence_InputIterator(PyObject* seq, Py_ssize_t index)
      : _seq(seq), _index(index)
    {
    }

    reference operator*() const
    {
      return reference(_seq, _index);
    }

    SwigPySequence_ArrowProxy<T>
    operator->() const {
      return SwigPySequence_ArrowProxy<T>(operator*());
    }

    bool operator==(const self& ri) const
    {
      return (_index == ri._index) && (_seq == ri._seq);
    }

    bool operator!=(const self& ri) const
    {
      return !(operator==(ri));
    }

    self& operator ++ ()
    {
      ++_index;
      return *this;
    }

    self& operator -- ()
    {
      --_index;
      return *this;
    }

    self& operator += (difference_type n)
    {
      _index += n;
      return *this;
    }

    self operator +(difference_type n) const
    {
      return self(_seq, _index + n);
    }

    self& operator -= (difference_type n)
    {
      _index -= n;
      return *this;
    }

    self operator -(difference_type n) const
    {
      return self(_seq, _index - n);
    }

    difference_type operator - (const self& ri) const
    {
      return _index - ri._index;
    }

    bool operator < (const self& ri) const
    {
      return _index < ri._index;
    }

    reference
    operator[](difference_type n) const
    {
      return reference(_seq, _index + n);
    }

  private:
    PyObject* _seq;
    difference_type _index;
  };

  // STL container wrapper around a Python sequence
  template <class T>
  struct SwigPySequence_Cont
  {
    typedef SwigPySequence_Ref<T> reference;
    typedef const SwigPySequence_Ref<T> const_reference;
    typedef T value_type;
    typedef T* pointer;
    typedef Py_ssize_t difference_type;
    typedef size_t size_type;
    typedef const pointer const_pointer;
    typedef SwigPySequence_InputIterator<T, reference> iterator;
    typedef SwigPySequence_InputIterator<T, const_reference> const_iterator;

    SwigPySequence_Cont(PyObject* seq) : _seq(0)
    {
      if (!PySequence_Check(seq)) {
	throw std::invalid_argument("a sequence is expected");
      }
      _seq = seq;
      Py_INCREF(_seq);
    }

    ~SwigPySequence_Cont()
    {
      Py_XDECREF(_seq);
    }

    size_type size() const
    {
      return static_cast<size_type>(PySequence_Size(_seq));
    }

    bool empty() const
    {
      return size() == 0;
    }

    iterator begin()
    {
      return iterator(_seq, 0);
    }

    const_iterator begin() const
    {
      return const_iterator(_seq, 0);
    }

    iterator end()
    {
      return iterator(_seq, size());
    }

    const_iterator end() const
    {
      return const_iterator(_seq, size());
    }

    reference operator[](difference_type n)
    {
      return reference(_seq, n);
    }

    const_reference operator[](difference_type n)  const
    {
      return const_reference(_seq, n);
    }

    bool check() const
    {
      Py_ssize_t s = size();
      for (Py_ssize_t i = 0; i < s; ++i) {
	swig::SwigVar_PyObject item = PySequence_GetItem(_seq, i);
	if (!swig::check<value_type>(item))
	  return false;
      }
      return true;
    }

  private:
    PyObject* _seq;
  };

}


  #define SWIG_From_double   PyFloat_FromDouble 


namespace swig {
  template <> struct traits< double > {
    typedef value_category category;
    static const char* type_name() { return"double"; }
  };
  template <>  struct traits_asval< double > {
    typedef double value_type;
    static int asval(PyObject *obj, value_type *val) {
      return SWIG_AsVal_double (obj, val);
    }
  };
  template <>  struct traits_from< double > {
    typedef double value_type;
    static PyObject *from(const value_type& val) {
      return SWIG_From_double  (val);
    }
  };
}


namespace swig {
  template <class SwigPySeq, class Seq>
  inline void
  assign(const SwigPySeq& swigpyseq, Seq* seq) {
    // seq->assign(swigpyseq.begin(), swigpyseq.end()); // not used as not always implemented
    typedef typename SwigPySeq::value_type value_type;
    typename SwigPySeq::const_iterator it = swigpyseq.begin();
    for (;it != swigpyseq.end(); ++it) {
      seq->insert(seq->end(),(value_type)(*it));
    }
  }

  template <class Seq, class T = typename Seq::value_type >
  struct traits_asptr_stdseq {
    typedef Seq sequence;
    typedef T value_type;

    static int asptr(PyObject *obj, sequence **seq) {
      if (obj == Py_None || SWIG_Python_GetSwigThis(obj)) {
	sequence *p;
	swig_type_info *descriptor = swig::type_info<sequence>();
	if (descriptor && SWIG_IsOK(::SWIG_ConvertPtr(obj, (void **)&p, descriptor, 0))) {
	  if (seq) *seq = p;
	  return SWIG_OLDOBJ;
	}
      } else if (PySequence_Check(obj)) {
	try {
	  SwigPySequence_Cont<value_type> swigpyseq(obj);
	  if (seq) {
	    sequence *pseq = new sequence();
	    assign(swigpyseq, pseq);
	    *seq = pseq;
	    return SWIG_NEWOBJ;
	  } else {
	    return swigpyseq.check() ? SWIG_OK : SWIG_ERROR;
	  }
	} catch (std::exception& e) {
	  if (seq) {
	    if (!PyErr_Occurred()) {
	      PyErr_SetString(PyExc_TypeError, e.what());
	    }
	  }
	  return SWIG_ERROR;
	}
      }
      return SWIG_ERROR;
    }
  };

  template <class Seq, class T = typename Seq::value_type >
  struct traits_from_stdseq {
    typedef Seq sequence;
    typedef T value_type;
    typedef typename Seq::size_type size_type;
    typedef typename sequence::const_iterator const_iterator;

    static PyObject *from(const sequence& seq) {
#ifdef SWIG_PYTHON_EXTRA_NATIVE_CONTAINERS
      swig_type_info *desc = swig::type_info<sequence>();
      if (desc && desc->clientdata) {
	return SWIG_InternalNewPointerObj(new sequence(seq), desc, SWIG_POINTER_OWN);
      }
#endif
      size_type size = seq.size();
      if (size <= (size_type)INT_MAX) {
	PyObject *obj = PyTuple_New((Py_ssize_t)size);
	Py_ssize_t i = 0;
	for (const_iterator it = seq.begin(); it != seq.end(); ++it, ++i) {
	  PyTuple_SetItem(obj,i,swig::from<value_type>(*it));
	}
	return obj;
      } else {
	PyErr_SetString(PyExc_OverflowError,"sequence size not valid in python");
	return NULL;
      }
    }
  };
}


  namespace swig {
    template <class T>
    struct traits_reserve<std::vector<T> > {
      static void reserve(std::vector<T> &seq, typename std::vector<T>::size_type n) {
        seq.reserve(n);
      }
    };

    template <class T>
    struct traits_asptr<std::vector<T> >  {
      static int asptr(PyObject *obj, std::vector<T> **vec) {
	return traits_asptr_stdseq<std::vector<T> >::asptr(obj, vec);
      }
    };
    
    template <class T>
    struct traits_from<std::vector<T> > {
      static PyObject *from(const std::vector<T>& vec) {
	return traits_from_stdseq<std::vector<T> >::from(vec);
      }
    };
  }


      namespace swig {
	template <>  struct traits<std::vector< double, std::allocator< double > > > {
	  typedef pointer_category category;
	  static const char* type_name() {
	    return "std::vector<" "double" "," "std::allocator< double >" " >";
	  }
	};
      }
    
SWIGINTERN swig::SwigPyIterator *std_vector_Sl_double_Sg__iterator(std::vector< double > *self,PyObject **PYTHON_SELF){
      return swig::make_output_iterator(self->begin(), self->begin(), self->end(), *PYTHON_SELF);
    }
SWIGINTERN bool std_vector_Sl_double_Sg____nonzero__(std::vector< double > const *self){
      return !(self->empty());
    }
SWIGINTERN bool std_vector_Sl_double_Sg____bool__(std::vector< double > const *self){
      return !(self->empty());
    }
SWIGINTERN std::vector< double >::size_type std_vector_Sl_double_Sg____len__(std::vector< double > const *self){
      return self->size();
    }

SWIGINTERNINLINE PyObject* 
SWIG_From_unsigned_SS_long  (unsigned long value)
{
  return (value > LONG_MAX) ?
    PyLong_FromUnsignedLong(value) : PyInt_FromLong(static_cast< long >(value));
}


#ifdef SWIG_LONG_LONG_AVAILABLE
SWIGINTERNINLINE PyObject* 
SWIG_From_unsigned_SS_long_SS_long  (unsigned long long value)
{
  return (value > LONG_MAX) ?
    PyLong_FromUnsignedLongLong(value) : PyInt_FromLong(static_cast< long >(value));
}
#endif


SWIGINTERNINLINE PyObject *
SWIG_From_size_t  (size_t value)
{    
#ifdef SWIG_LONG_LONG_AVAILABLE
  if (sizeof(size_t) <= sizeof(unsigned long)) {
#endif
    return SWIG_From_unsigned_SS_long  (static_cast< unsigned long >(value));
#ifdef SWIG_LONG_LONG_AVAILABLE
  } else {
    /* assume sizeof(size_t) <= sizeof(unsigned long long) */
    return SWIG_From_unsigned_SS_long_SS_long  (static_cast< unsigned long long >(value));
  }
#endif
}

SWIGINTERN std::vector< double,std::allocator< double > > *std_vector_Sl_double_Sg____getslice__(std::vector< double > *self,std::vector< double >::difference_type i,std::vector< double >::difference_type j){
      return swig::getslice(self, i, j, 1);
    }
SWIGINTERN void std_vector_Sl_double_Sg____setslice____SWIG_0(std::vector< double > *self,std::vector< double >::difference_type i,std::vector< double >::difference_type j){
      swig::setslice(self, i, j, 1, std::vector< double,std::allocator< double > >());
    }
SWIGINTERN void std_vector_Sl_double_Sg____setslice____SWIG_1(std::vector< double > *self,std::vector< double >::difference_type i,std::vector< double >::difference_type j,std::vector< double,std::allocator< double > > const &v){
      swig::setslice(self, i, j, 1, v);
    }
SWIGINTERN void std_vector_Sl_double_Sg____delslice__(std::vector< double > *self,std::vector< double >::difference_type i,std::vector< double >::difference_type j){
      swig::delslice(self, i, j, 1);
    }
SWIGINTERN void std_vector_Sl_double_Sg____delitem____SWIG_0(std::vector< double > *self,std::vector< double >::difference_type i){
      swig::erase(self, swig::getpos(self, i));
    }
SWIGINTERN std::vector< double,std::allocator< double > > *std_vector_Sl_double_Sg____getitem____SWIG_0(std::vector< double > *self,PySliceObject *slice){
      Py_ssize_t i, j, step;
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
        return NULL;
      }
      PySlice_GetIndices(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &i, &j, &step);
      std::vector< double,std::allocator< double > >::difference_type id = i;
      std::vector< double,std::allocator< double > >::difference_type jd = j;
      return swig::getslice(self, id, jd, step);
    }
SWIGINTERN void std_vector_Sl_double_Sg____setitem____SWIG_0(std::vector< double > *self,PySliceObject *slice,std::vector< double,std::allocator< double > > const &v){
      Py_ssize_t i, j, step;
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
        return;
      }
      PySlice_GetIndices(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &i, &j, &step);
      std::vector< double,std::allocator< double > >::difference_type id = i;
      std::vector< double,std::allocator< double > >::difference_type jd = j;
      swig::setslice(self, id, jd, step, v);
    }
SWIGINTERN void std_vector_Sl_double_Sg____setitem____SWIG_1(std::vector< double > *self,PySliceObject *slice){
      Py_ssize_t i, j, step;
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
        return;
      }
      PySlice_GetIndices(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &i, &j, &step);
      std::vector< double,std::allocator< double > >::difference_type id = i;
      std::vector< double,std::allocator< double > >::difference_type jd = j;
      swig::delslice(self, id, jd, step);
    }
SWIGINTERN void std_vector_Sl_double_Sg____delitem____SWIG_1(std::vector< double > *self,PySliceObject *slice){
      Py_ssize_t i, j, step;
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
        return;
      }
      PySlice_GetIndices(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &i, &j, &step);
      std::vector< double,std::allocator< double > >::difference_type id = i;
      std::vector< double,std::allocator< double > >::difference_type jd = j;
      swig::delslice(self, id, jd, step);
    }
SWIGINTERN std::vector< double >::value_type const &std_vector_Sl_double_Sg____getitem____SWIG_1(std::vector< double > const *self,std::vector< double >::difference_type i){
      return *(swig::cgetpos(self, i));
    }

namespace swig {
  static PyObject* container_owner_attribute() {
    static PyObject* attr = SWIG_Python_str_FromChar("__swig_container");
    return attr;
  }

  template <typename T>
  struct container_owner {
    // By default, do not add the back-reference (for value types)
    // Specialization below will check the reference for pointer types.
    static bool back_reference(PyObject* /*child*/, PyObject* /*owner*/) {
      return false;
    }
  };

  template <>
  struct container_owner<swig::pointer_category> {  
    /*
     * Call to add a back-reference to the owning object when returning a 
     * reference from a container.  Will only set the reference if child
     * is a SWIG wrapper object that does not own the pointer.
     *
     * returns whether the reference was set or not
     */
    static bool back_reference(PyObject* child, PyObject* owner) {
      SwigPyObject* swigThis = SWIG_Python_GetSwigThis(child);
      if (swigThis && (swigThis->own & SWIG_POINTER_OWN) != SWIG_POINTER_OWN) {
        return PyObject_SetAttr(child, container_owner_attribute(), owner) != -1;
      }
      return false;
    }
  };
}

SWIGINTERN void std_vector_Sl_double_Sg____setitem____SWIG_2(std::vector< double > *self,std::vector< double >::difference_type i,std::vector< double >::value_type const &x){
      *(swig::getpos(self,i)) = x;
    }
SWIGINTERN std::vector< double >::value_type std_vector_Sl_double_Sg__pop(std::vector< double > *self){
      if (self->size() == 0)
	throw std::out_of_range("pop from empty container");
      std::vector< double,std::allocator< double > >::value_type x = self->back();
      self->pop_back();
      return x;
    }
SWIGINTERN void std_vector_Sl_double_Sg__append(std::vector< double > *self,std::vector< double >::value_type const &x){
      self->push_back(x);
    }
SWIGINTERN std::vector< double >::iterator std_vector_Sl_double_Sg__erase__SWIG_0(std::vector< double > *self,std::vector< double >::iterator pos){ return self->erase(pos); }
SWIGINTERN std::vector< double >::iterator std_vector_Sl_double_Sg__erase__SWIG_1(std::vector< double > *self,std::vector< double >::iterator first,std::vector< double >::iterator last){ return self->erase(first, last); }
SWIGINTERN std::vector< double >::iterator std_vector_Sl_double_Sg__insert__SWIG_0(std::vector< double > *self,std::vector< double >::iterator pos,std::vector< double >::value_type const &x){ return self->insert(pos, x); }
SWIGINTERN void std_vector_Sl_double_Sg__insert__SWIG_1(std::vector< double > *self,std::vector< double >::iterator pos,std::vector< double >::size_type n,std::vector< double >::value_type const &x){ self->insert(pos, n, x); }

SWIGINTERN swig_type_info*
SWIG_pchar_descriptor(void)
{
  static int init = 0;
  static swig_type_info* info = 0;
  if (!init) {
    info = SWIG_TypeQuery("_p_char");
    init = 1;
  }
  return info;
}


SWIGINTERN int
SWIG_AsCharPtrAndSize(PyObject *obj, char** cptr, size_t* psize, int *alloc)
{
#if PY_VERSION_HEX>=0x03000000
#if defined(SWIG_PYTHON_STRICT_BYTE_CHAR)
  if (PyBytes_Check(obj))
#else
  if (PyUnicode_Check(obj))
#endif
#else  
  if (PyString_Check(obj))
#endif
  {
    char *cstr; Py_ssize_t len;
    int ret = SWIG_OK;
#if PY_VERSION_HEX>=0x03000000
#if !defined(SWIG_PYTHON_STRICT_BYTE_CHAR)
    if (!alloc && cptr) {
        /* We can't allow converting without allocation, since the internal
           representation of string in Python 3 is UCS-2/UCS-4 but we require
           a UTF-8 representation.
           TODO(bhy) More detailed explanation */
        return SWIG_RuntimeError;
    }
    obj = PyUnicode_AsUTF8String(obj);
    if (!obj)
      return SWIG_TypeError;
    if (alloc)
      *alloc = SWIG_NEWOBJ;
#endif
    if (PyBytes_AsStringAndSize(obj, &cstr, &len) == -1)
      return SWIG_TypeError;
#else
    if (PyString_AsStringAndSize(obj, &cstr, &len) == -1)
      return SWIG_TypeError;
#endif
    if (cptr) {
      if (alloc) {
	if (*alloc == SWIG_NEWOBJ) {
	  *cptr = reinterpret_cast< char* >(memcpy(new char[len + 1], cstr, sizeof(char)*(len + 1)));
	  *alloc = SWIG_NEWOBJ;
	} else {
	  *cptr = cstr;
	  *alloc = SWIG_OLDOBJ;
	}
      } else {
#if PY_VERSION_HEX>=0x03000000
#if defined(SWIG_PYTHON_STRICT_BYTE_CHAR)
	*cptr = PyBytes_AsString(obj);
#else
	assert(0); /* Should never reach here with Unicode strings in Python 3 */
#endif
#else
	*cptr = SWIG_Python_str_AsChar(obj);
        if (!*cptr)
          ret = SWIG_TypeError;
#endif
      }
    }
    if (psize) *psize = len + 1;
#if PY_VERSION_HEX>=0x03000000 && !defined(SWIG_PYTHON_STRICT_BYTE_CHAR)
    Py_XDECREF(obj);
#endif
    return ret;
  } else {
#if defined(SWIG_PYTHON_2_UNICODE)
#if defined(SWIG_PYTHON_STRICT_BYTE_CHAR)
#error "Cannot use both SWIG_PYTHON_2_UNICODE and SWIG_PYTHON_STRICT_BYTE_CHAR at once"
#endif
#if PY_VERSION_HEX<0x03000000
    if (PyUnicode_Check(obj)) {
      char *cstr; Py_ssize_t len;
      if (!alloc && cptr) {
        return SWIG_RuntimeError;
      }
      obj = PyUnicode_AsUTF8String(obj);
      if (!obj)
        return SWIG_TypeError;
      if (PyString_AsStringAndSize(obj, &cstr, &len) != -1) {
        if (cptr) {
          if (alloc) *alloc = SWIG_NEWOBJ;
          *cptr = reinterpret_cast< char* >(memcpy(new char[len + 1], cstr, sizeof(char)*(len + 1)));
        }
        if (psize) *psize = len + 1;

        Py_XDECREF(obj);
        return SWIG_OK;
      } else {
        Py_XDECREF(obj);
      }
    }
#endif
#endif

    swig_type_info* pchar_descriptor = SWIG_pchar_descriptor();
    if (pchar_descriptor) {
      void* vptr = 0;
      if (SWIG_ConvertPtr(obj, &vptr, pchar_descriptor, 0) == SWIG_OK) {
	if (cptr) *cptr = (char *) vptr;
	if (psize) *psize = vptr ? (strlen((char *)vptr) + 1) : 0;
	if (alloc) *alloc = SWIG_OLDOBJ;
	return SWIG_OK;
      }
    }
  }
  return SWIG_TypeError;
}


SWIGINTERN int
SWIG_AsPtr_std_string (PyObject * obj, std::string **val) 
{
  char* buf = 0 ; size_t size = 0; int alloc = SWIG_OLDOBJ;
  if (SWIG_IsOK((SWIG_AsCharPtrAndSize(obj, &buf, &size, &alloc)))) {
    if (buf) {
      if (val) *val = new std::string(buf, size - 1);
      if (alloc == SWIG_NEWOBJ) delete[] buf;
      return SWIG_NEWOBJ;
    } else {
      if (val) *val = 0;
      return SWIG_OLDOBJ;
    }
  } else {
    static int init = 0;
    static swig_type_info* descriptor = 0;
    if (!init) {
      descriptor = SWIG_TypeQuery("std::string" " *");
      init = 1;
    }
    if (descriptor) {
      std::string *vptr;
      int res = SWIG_ConvertPtr(obj, (void**)&vptr, descriptor, 0);
      if (SWIG_IsOK(res) && val) *val = vptr;
      return res;
    }
  }
  return SWIG_ERROR;
}


SWIGINTERN int
SWIG_AsVal_std_string (PyObject * obj, std::string *val)
{
  std::string* v = (std::string *) 0;
  int res = SWIG_AsPtr_std_string (obj, &v);
  if (!SWIG_IsOK(res)) return res;
  if (v) {
    if (val) *val = *v;
    if (SWIG_IsNewObj(res)) {
      delete v;
      res = SWIG_DelNewMask(res);
    }
    return res;
  }
  return SWIG_ERROR;
}


SWIGINTERNINLINE PyObject *
SWIG_FromCharPtrAndSize(const char* carray, size_t size)
{
  if (carray) {
    if (size > INT_MAX) {
      swig_type_info* pchar_descriptor = SWIG_pchar_descriptor();
      return pchar_descriptor ? 
	SWIG_InternalNewPointerObj(const_cast< char * >(carray), pchar_descriptor, 0) : SWIG_Py_Void();
    } else {
#if PY_VERSION_HEX >= 0x03000000
#if defined(SWIG_PYTHON_STRICT_BYTE_CHAR)
      return PyBytes_FromStringAndSize(carray, static_cast< Py_ssize_t >(size));
#else
      return PyUnicode_DecodeUTF8(carray, static_cast< Py_ssize_t >(size), "surrogateescape");
#endif
#else
      return PyString_FromStringAndSize(carray, static_cast< Py_ssize_t >(size));
#endif
    }
  } else {
    return SWIG_Py_Void();
  }
}


SWIGINTERNINLINE PyObject *
SWIG_From_std_string  (const std::string& s)
{
  return SWIG_FromCharPtrAndSize(s.data(), s.size());
}


namespace swig {
  template <> struct traits< std::string > {
    typedef value_category category;
    static const char* type_name() { return"std::string"; }
  };
  template <>  struct traits_asval< std::string > {
    typedef std::string value_type;
    static int asval(PyObject *obj, value_type *val) {
      return SWIG_AsVal_std_string (obj, val);
    }
  };
  template <>  struct traits_from< std::string > {
    typedef std::string value_type;
    static PyObject *from(const value_type& val) {
      return SWIG_From_std_string  (val);
    }
  };
}


      namespace swig {
	template <>  struct traits<std::vector< std::string, std::allocator< std::string > > > {
	  typedef pointer_category category;
	  static const char* type_name() {
	    return "std::vector<" "std::string" "," "std::allocator< std::string >" " >";
	  }
	};
      }
    
SWIGINTERN swig::SwigPyIterator *std_vector_Sl_std_string_Sg__iterator(std::vector< std::string > *self,PyObject **PYTHON_SELF){
      return swig::make_output_iterator(self->begin(), self->begin(), self->end(), *PYTHON_SELF);
    }
SWIGINTERN bool std_vector_Sl_std_string_Sg____nonzero__(std::vector< std::string > const *self){
      return !(self->empty());
    }
SWIGINTERN bool std_vector_Sl_std_string_Sg____bool__(std::vector< std::string > const *self){
      return !(self->empty());
    }
SWIGINTERN std::vector< std::string >::size_type std_vector_Sl_std_string_Sg____len__(std::vector< std::string > const *self){
      return self->size();
    }
SWIGINTERN std::vector< std::string,std::allocator< std::string > > *std_vector_Sl_std_string_Sg____getslice__(std::vector< std::string > *self,std::vector< std::string >::difference_type i,std::vector< std::string >::difference_type j){
      return swig::getslice(self, i, j, 1);
    }
SWIGINTERN void std_vector_Sl_std_string_Sg____setslice____SWIG_0(std::vector< std::string > *self,std::vector< std::string >::difference_type i,std::vector< std::string >::difference_type j){
      swig::setslice(self, i, j, 1, std::vector< std::string,std::allocator< std::string > >());
    }
SWIGINTERN void std_vector_Sl_std_string_Sg____setslice____SWIG_1(std::vector< std::string > *self,std::vector< std::string >::difference_type i,std::vector< std::string >::difference_type j,std::vector< std::string,std::allocator< std::string > > const &v){
      swig::setslice(self, i, j, 1, v);
    }
SWIGINTERN void std_vector_Sl_std_string_Sg____delslice__(std::vector< std::string > *self,std::vector< std::string >::difference_type i,std::vector< std::string >::difference_type j){
      swig::delslice(self, i, j, 1);
    }
SWIGINTERN void std_vector_Sl_std_string_Sg____delitem____SWIG_0(std::vector< std::string > *self,std::vector< std::string >::difference_type i){
      swig::erase(self, swig::getpos(self, i));
    }
SWIGINTERN std::vector< std::string,std::allocator< std::string > > *std_vector_Sl_std_string_Sg____getitem____SWIG_0(std::vector< std::string > *self,PySliceObject *slice){
      Py_ssize_t i, j, step;
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
        return NULL;
      }
      PySlice_GetIndices(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &i, &j, &step);
      std::vector< std::string,std::allocator< std::string > >::difference_type id = i;
      std::vector< std::string,std::allocator< std::string > >::difference_type jd = j;
      return swig::getslice(self, id, jd, step);
    }
SWIGINTERN void std_vector_Sl_std_string_Sg____setitem____SWIG_0(std::vector< std::string > *self,PySliceObject *slice,std::vector< std::string,std::allocator< std::string > > const &v){
      Py_ssize_t i, j, step;
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
        return;
      }
      PySlice_GetIndices(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &i, &j, &step);
      std::vector< std::string,std::allocator< std::string > >::difference_type id = i;
      std::vector< std::string,std::allocator< std::string > >::difference_type jd = j;
      swig::setslice(self, id, jd, step, v);
    }
SWIGINTERN void std_vector_Sl_std_string_Sg____setitem____SWIG_1(std::vector< std::string > *self,PySliceObject *slice){
      Py_ssize_t i, j, step;
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
        return;
      }
      PySlice_GetIndices(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &i, &j, &step);
      std::vector< std::string,std::allocator< std::string > >::difference_type id = i;
      std::vector< std::string,std::allocator< std::string > >::difference_type jd = j;
      swig::delslice(self, id, jd, step);
    }
SWIGINTERN void std_vector_Sl_std_string_Sg____delitem____SWIG_1(std::vector< std::string > *self,PySliceObject *slice){
      Py_ssize_t i, j, step;
      if( !PySlice_Check(slice) ) {
        SWIG_Error(SWIG_TypeError, "Slice object expected.");
        return;
      }
      PySlice_GetIndices(SWIGPY_SLICE_ARG(slice), (Py_ssize_t)self->size(), &i, &j, &step);
      std::vector< std::string,std::allocator< std::string > >::difference_type id = i;
      std::vector< std::string,std::allocator< std::string > >::difference_type jd = j;
      swig::delslice(self, id, jd, step);
    }
SWIGINTERN std::vector< std::string >::value_type const &std_vector_Sl_std_string_Sg____getitem____SWIG_1(std::vector< std::string > const *self,std::vector< std::string >::difference_type i){
      return *(swig::cgetpos(self, i));
    }
SWIGINTERN void std_vector_Sl_std_string_Sg____setitem____SWIG_2(std::vector< std::string > *self,std::vector< std::string >::difference_type i,std::vector< std::string >::value_type const &x){
      *(swig::getpos(self,i)) = x;
    }
SWIGINTERN std::vector< std::string >::value_type std_vector_Sl_std_string_Sg__pop(std::vector< std::string > *self){
      if (self->size() == 0)
	throw std::out_of_range("pop from empty container");
      std::vector< std::string,std::allocator< std::string > >::value_type x = self->back();
      self->pop_back();
      return x;
    }
SWIGINTERN void std_vector_Sl_std_string_Sg__append(std::vector< std::string > *self,std::vector< std::string >::value_type const &x){
      self->push_back(x);
    }
SWIGINTERN std::vector< std::string >::iterator std_vector_Sl_std_string_Sg__erase__SWIG_0(std::vector< std::string > *self,std::vector< std::string >::iterator pos){ return self->erase(pos); }
SWIGINTERN std::vector< std::string >::iterator std_vector_Sl_std_string_Sg__erase__SWIG_1(std::vector< std::string > *self,std::vector< std::string >::iterator first,std::vector< std::string >::iterator last){ return self->erase(first, last); }
SWIGINTERN std::vector< std::string >::iterator std_vector_Sl_std_string_Sg__insert__SWIG_0(std::vector< std::string > *self,std::vector< std::string >::iterator pos,std::vector< std::string >::value_type const &x){ return self->insert(pos, x); }
SWIGINTERN void std_vector_Sl_std_string_Sg__insert__SWIG_1(std::vector< std::string > *self,std::vector< std::string >::iterator pos,std::vector< std::string >::size_type n,std::vector< std::string >::value_type const &x){ self->insert(pos, n, x); }

SWIGINTERN int
SWIG_AsVal_int (PyObject * obj, int *val)
{
  long v;
  int res = SWIG_AsVal_long (obj, &v);
  if (SWIG_IsOK(res)) {
    if ((v < INT_MIN || v > INT_MAX)) {
      return SWIG_OverflowError;
    } else {
      if (val) *val = static_cast< int >(v);
    }
  }  
  return res;
}


std::vector<double> fractional_black_scholes_alpha(const std::vector<double>& S, double alpha, double sigma, double r,
                                                   double dx, int N, int L, double K, double T, const std::string& option_type) {
    return ::fractional_black_scholes_alpha(S, alpha, sigma, r, dx, N, L, K, T, option_type);
}

// Wrap the calibration objective function
double calibration_objective(double alpha, const std::vector<std::vector<double>>& options_data,
                             const std::string& option_type, const std::vector<double>& S, double sigma,
                             double r, double dx, int N, int L, double S0) {
    return ::calibration_objective(alpha, options_data, option_type, S, sigma, r, dx, N, L, S0);
}

#ifdef __cplusplus
extern "C" {
#endif
SWIGINTERN PyObject *_wrap_delete_SwigPyIterator(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_DISOWN |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "delete_SwigPyIterator" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  delete arg1;
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_value(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  PyObject *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_value" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  try {
    result = (PyObject *)((swig::SwigPyIterator const *)arg1)->value();
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = result;
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_incr__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  size_t arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  swig::SwigPyIterator *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_incr" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "SwigPyIterator_incr" "', argument " "2"" of type '" "size_t""'");
  } 
  arg2 = static_cast< size_t >(val2);
  try {
    result = (swig::SwigPyIterator *)(arg1)->incr(arg2);
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_incr__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *result = 0 ;
  
  if ((nobjs < 1) || (nobjs > 1)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_incr" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  try {
    result = (swig::SwigPyIterator *)(arg1)->incr();
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_incr(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "SwigPyIterator_incr", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 1) {
    int _v;
    void *vptr = 0;
    int res = SWIG_ConvertPtr(argv[0], &vptr, SWIGTYPE_p_swig__SwigPyIterator, 0);
    _v = SWIG_CheckState(res);
    if (_v) {
      return _wrap_SwigPyIterator_incr__SWIG_1(self, argc, argv);
    }
  }
  if (argc == 2) {
    int _v;
    void *vptr = 0;
    int res = SWIG_ConvertPtr(argv[0], &vptr, SWIGTYPE_p_swig__SwigPyIterator, 0);
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_size_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_SwigPyIterator_incr__SWIG_0(self, argc, argv);
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'SwigPyIterator_incr'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    swig::SwigPyIterator::incr(size_t)\n"
    "    swig::SwigPyIterator::incr()\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_decr__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  size_t arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  swig::SwigPyIterator *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_decr" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "SwigPyIterator_decr" "', argument " "2"" of type '" "size_t""'");
  } 
  arg2 = static_cast< size_t >(val2);
  try {
    result = (swig::SwigPyIterator *)(arg1)->decr(arg2);
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_decr__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *result = 0 ;
  
  if ((nobjs < 1) || (nobjs > 1)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_decr" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  try {
    result = (swig::SwigPyIterator *)(arg1)->decr();
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_decr(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "SwigPyIterator_decr", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 1) {
    int _v;
    void *vptr = 0;
    int res = SWIG_ConvertPtr(argv[0], &vptr, SWIGTYPE_p_swig__SwigPyIterator, 0);
    _v = SWIG_CheckState(res);
    if (_v) {
      return _wrap_SwigPyIterator_decr__SWIG_1(self, argc, argv);
    }
  }
  if (argc == 2) {
    int _v;
    void *vptr = 0;
    int res = SWIG_ConvertPtr(argv[0], &vptr, SWIGTYPE_p_swig__SwigPyIterator, 0);
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_size_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_SwigPyIterator_decr__SWIG_0(self, argc, argv);
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'SwigPyIterator_decr'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    swig::SwigPyIterator::decr(size_t)\n"
    "    swig::SwigPyIterator::decr()\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_distance(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  swig::SwigPyIterator *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  void *argp2 = 0 ;
  int res2 = 0 ;
  PyObject *swig_obj[2] ;
  ptrdiff_t result;
  
  if (!SWIG_Python_UnpackTuple(args, "SwigPyIterator_distance", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_distance" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], &argp2, SWIGTYPE_p_swig__SwigPyIterator,  0  | 0);
  if (!SWIG_IsOK(res2)) {
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "SwigPyIterator_distance" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  if (!argp2) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "SwigPyIterator_distance" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  arg2 = reinterpret_cast< swig::SwigPyIterator * >(argp2);
  try {
    result = ((swig::SwigPyIterator const *)arg1)->distance((swig::SwigPyIterator const &)*arg2);
  } catch(std::invalid_argument &_e) {
    SWIG_Python_Raise(SWIG_NewPointerObj((new std::invalid_argument(static_cast< const std::invalid_argument& >(_e))),SWIGTYPE_p_std__invalid_argument,SWIG_POINTER_OWN), "std::invalid_argument", SWIGTYPE_p_std__invalid_argument); SWIG_fail;
  }
  resultobj = SWIG_From_ptrdiff_t(static_cast< ptrdiff_t >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_equal(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  swig::SwigPyIterator *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  void *argp2 = 0 ;
  int res2 = 0 ;
  PyObject *swig_obj[2] ;
  bool result;
  
  if (!SWIG_Python_UnpackTuple(args, "SwigPyIterator_equal", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_equal" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], &argp2, SWIGTYPE_p_swig__SwigPyIterator,  0  | 0);
  if (!SWIG_IsOK(res2)) {
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "SwigPyIterator_equal" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  if (!argp2) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "SwigPyIterator_equal" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  arg2 = reinterpret_cast< swig::SwigPyIterator * >(argp2);
  try {
    result = (bool)((swig::SwigPyIterator const *)arg1)->equal((swig::SwigPyIterator const &)*arg2);
  } catch(std::invalid_argument &_e) {
    SWIG_Python_Raise(SWIG_NewPointerObj((new std::invalid_argument(static_cast< const std::invalid_argument& >(_e))),SWIGTYPE_p_std__invalid_argument,SWIG_POINTER_OWN), "std::invalid_argument", SWIGTYPE_p_std__invalid_argument); SWIG_fail;
  }
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_copy(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  swig::SwigPyIterator *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_copy" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  result = (swig::SwigPyIterator *)((swig::SwigPyIterator const *)arg1)->copy();
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_next(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  PyObject *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_next" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  try {
    result = (PyObject *)(arg1)->next();
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = result;
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___next__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  PyObject *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator___next__" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  try {
    result = (PyObject *)(arg1)->__next__();
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = result;
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_previous(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  PyObject *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_previous" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  try {
    result = (PyObject *)(arg1)->previous();
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = result;
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator_advance(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  ptrdiff_t arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  swig::SwigPyIterator *result = 0 ;
  
  if (!SWIG_Python_UnpackTuple(args, "SwigPyIterator_advance", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator_advance" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "SwigPyIterator_advance" "', argument " "2"" of type '" "ptrdiff_t""'");
  } 
  arg2 = static_cast< ptrdiff_t >(val2);
  try {
    result = (swig::SwigPyIterator *)(arg1)->advance(arg2);
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___eq__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  swig::SwigPyIterator *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  void *argp2 = 0 ;
  int res2 = 0 ;
  PyObject *swig_obj[2] ;
  bool result;
  
  if (!SWIG_Python_UnpackTuple(args, "SwigPyIterator___eq__", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator___eq__" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], &argp2, SWIGTYPE_p_swig__SwigPyIterator,  0  | 0);
  if (!SWIG_IsOK(res2)) {
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "SwigPyIterator___eq__" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  if (!argp2) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "SwigPyIterator___eq__" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  arg2 = reinterpret_cast< swig::SwigPyIterator * >(argp2);
  result = (bool)((swig::SwigPyIterator const *)arg1)->operator ==((swig::SwigPyIterator const &)*arg2);
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  PyErr_Clear();
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___ne__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  swig::SwigPyIterator *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  void *argp2 = 0 ;
  int res2 = 0 ;
  PyObject *swig_obj[2] ;
  bool result;
  
  if (!SWIG_Python_UnpackTuple(args, "SwigPyIterator___ne__", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator___ne__" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], &argp2, SWIGTYPE_p_swig__SwigPyIterator,  0  | 0);
  if (!SWIG_IsOK(res2)) {
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "SwigPyIterator___ne__" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  if (!argp2) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "SwigPyIterator___ne__" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  arg2 = reinterpret_cast< swig::SwigPyIterator * >(argp2);
  result = (bool)((swig::SwigPyIterator const *)arg1)->operator !=((swig::SwigPyIterator const &)*arg2);
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  PyErr_Clear();
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___iadd__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  ptrdiff_t arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  swig::SwigPyIterator *result = 0 ;
  
  if (!SWIG_Python_UnpackTuple(args, "SwigPyIterator___iadd__", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_DISOWN |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator___iadd__" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "SwigPyIterator___iadd__" "', argument " "2"" of type '" "ptrdiff_t""'");
  } 
  arg2 = static_cast< ptrdiff_t >(val2);
  try {
    result = (swig::SwigPyIterator *) &(arg1)->operator +=(arg2);
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___isub__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  ptrdiff_t arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  swig::SwigPyIterator *result = 0 ;
  
  if (!SWIG_Python_UnpackTuple(args, "SwigPyIterator___isub__", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_DISOWN |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator___isub__" "', argument " "1"" of type '" "swig::SwigPyIterator *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "SwigPyIterator___isub__" "', argument " "2"" of type '" "ptrdiff_t""'");
  } 
  arg2 = static_cast< ptrdiff_t >(val2);
  try {
    result = (swig::SwigPyIterator *) &(arg1)->operator -=(arg2);
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___add__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  ptrdiff_t arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  swig::SwigPyIterator *result = 0 ;
  
  if (!SWIG_Python_UnpackTuple(args, "SwigPyIterator___add__", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator___add__" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "SwigPyIterator___add__" "', argument " "2"" of type '" "ptrdiff_t""'");
  } 
  arg2 = static_cast< ptrdiff_t >(val2);
  try {
    result = (swig::SwigPyIterator *)((swig::SwigPyIterator const *)arg1)->operator +(arg2);
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  PyErr_Clear();
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___sub____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  ptrdiff_t arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  swig::SwigPyIterator *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator___sub__" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "SwigPyIterator___sub__" "', argument " "2"" of type '" "ptrdiff_t""'");
  } 
  arg2 = static_cast< ptrdiff_t >(val2);
  try {
    result = (swig::SwigPyIterator *)((swig::SwigPyIterator const *)arg1)->operator -(arg2);
  } catch(swig::stop_iteration &_e) {
    {
      (void)_e;
      SWIG_SetErrorObj(PyExc_StopIteration, SWIG_Py_Void());
      SWIG_fail;
    }
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  PyErr_Clear();
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___sub____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  swig::SwigPyIterator *arg1 = (swig::SwigPyIterator *) 0 ;
  swig::SwigPyIterator *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  void *argp2 = 0 ;
  int res2 = 0 ;
  ptrdiff_t result;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_swig__SwigPyIterator, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "SwigPyIterator___sub__" "', argument " "1"" of type '" "swig::SwigPyIterator const *""'"); 
  }
  arg1 = reinterpret_cast< swig::SwigPyIterator * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], &argp2, SWIGTYPE_p_swig__SwigPyIterator,  0  | 0);
  if (!SWIG_IsOK(res2)) {
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "SwigPyIterator___sub__" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  if (!argp2) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "SwigPyIterator___sub__" "', argument " "2"" of type '" "swig::SwigPyIterator const &""'"); 
  }
  arg2 = reinterpret_cast< swig::SwigPyIterator * >(argp2);
  result = ((swig::SwigPyIterator const *)arg1)->operator -((swig::SwigPyIterator const &)*arg2);
  resultobj = SWIG_From_ptrdiff_t(static_cast< ptrdiff_t >(result));
  return resultobj;
fail:
  PyErr_Clear();
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}


SWIGINTERN PyObject *_wrap_SwigPyIterator___sub__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "SwigPyIterator___sub__", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    void *vptr = 0;
    int res = SWIG_ConvertPtr(argv[0], &vptr, SWIGTYPE_p_swig__SwigPyIterator, 0);
    _v = SWIG_CheckState(res);
    if (_v) {
      int res = SWIG_ConvertPtr(argv[1], 0, SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_NO_NULL | 0);
      _v = SWIG_CheckState(res);
      if (_v) {
        return _wrap_SwigPyIterator___sub____SWIG_1(self, argc, argv);
      }
    }
  }
  if (argc == 2) {
    int _v;
    void *vptr = 0;
    int res = SWIG_ConvertPtr(argv[0], &vptr, SWIGTYPE_p_swig__SwigPyIterator, 0);
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_SwigPyIterator___sub____SWIG_0(self, argc, argv);
      }
    }
  }
  
fail:
  Py_INCREF(Py_NotImplemented);
  return Py_NotImplemented;
}


SWIGINTERN PyObject *SwigPyIterator_swigregister(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *obj;
  if (!SWIG_Python_UnpackTuple(args, "swigregister", 1, 1, &obj)) return NULL;
  SWIG_TypeNewClientData(SWIGTYPE_p_swig__SwigPyIterator, SWIG_NewClientData(obj));
  return SWIG_Py_Void();
}

SWIGINTERN PyObject *_wrap_VectorDouble_iterator(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  PyObject **arg2 = (PyObject **) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  swig::SwigPyIterator *result = 0 ;
  
  arg2 = &swig_obj[0];
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_iterator" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (swig::SwigPyIterator *)std_vector_Sl_double_Sg__iterator(arg1,arg2);
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___nonzero__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  bool result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___nonzero__" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (bool)std_vector_Sl_double_Sg____nonzero__((std::vector< double > const *)arg1);
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___bool__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  bool result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___bool__" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (bool)std_vector_Sl_double_Sg____bool__((std::vector< double > const *)arg1);
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___len__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::size_type result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___len__" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = std_vector_Sl_double_Sg____len__((std::vector< double > const *)arg1);
  resultobj = SWIG_From_size_t(static_cast< size_t >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___getslice__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::difference_type arg2 ;
  std::vector< double >::difference_type arg3 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  ptrdiff_t val3 ;
  int ecode3 = 0 ;
  PyObject *swig_obj[3] ;
  std::vector< double,std::allocator< double > > *result = 0 ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorDouble___getslice__", 3, 3, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___getslice__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble___getslice__" "', argument " "2"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< double >::difference_type >(val2);
  ecode3 = SWIG_AsVal_ptrdiff_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble___getslice__" "', argument " "3"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg3 = static_cast< std::vector< double >::difference_type >(val3);
  try {
    result = (std::vector< double,std::allocator< double > > *)std_vector_Sl_double_Sg____getslice__(arg1,arg2,arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___setslice____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::difference_type arg2 ;
  std::vector< double >::difference_type arg3 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  ptrdiff_t val3 ;
  int ecode3 = 0 ;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___setslice__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble___setslice__" "', argument " "2"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< double >::difference_type >(val2);
  ecode3 = SWIG_AsVal_ptrdiff_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble___setslice__" "', argument " "3"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg3 = static_cast< std::vector< double >::difference_type >(val3);
  try {
    std_vector_Sl_double_Sg____setslice____SWIG_0(arg1,arg2,arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___setslice____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::difference_type arg2 ;
  std::vector< double >::difference_type arg3 ;
  std::vector< double,std::allocator< double > > *arg4 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  ptrdiff_t val3 ;
  int ecode3 = 0 ;
  int res4 = SWIG_OLDOBJ ;
  
  if ((nobjs < 4) || (nobjs > 4)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___setslice__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble___setslice__" "', argument " "2"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< double >::difference_type >(val2);
  ecode3 = SWIG_AsVal_ptrdiff_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble___setslice__" "', argument " "3"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg3 = static_cast< std::vector< double >::difference_type >(val3);
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res4 = swig::asptr(swig_obj[3], &ptr);
    if (!SWIG_IsOK(res4)) {
      SWIG_exception_fail(SWIG_ArgError(res4), "in method '" "VectorDouble___setslice__" "', argument " "4"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorDouble___setslice__" "', argument " "4"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg4 = ptr;
  }
  try {
    std_vector_Sl_double_Sg____setslice____SWIG_1(arg1,arg2,arg3,(std::vector< double,std::allocator< double > > const &)*arg4);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res4)) delete arg4;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res4)) delete arg4;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___setslice__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[5] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorDouble___setslice__", 0, 4, argv))) SWIG_fail;
  --argc;
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        {
          int res = SWIG_AsVal_ptrdiff_t(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          return _wrap_VectorDouble___setslice____SWIG_0(self, argc, argv);
        }
      }
    }
  }
  if (argc == 4) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        {
          int res = SWIG_AsVal_ptrdiff_t(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          int res = swig::asptr(argv[3], (std::vector< double,std::allocator< double > >**)(0));
          _v = SWIG_CheckState(res);
          if (_v) {
            return _wrap_VectorDouble___setslice____SWIG_1(self, argc, argv);
          }
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorDouble___setslice__'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< double >::__setslice__(std::vector< double >::difference_type,std::vector< double >::difference_type)\n"
    "    std::vector< double >::__setslice__(std::vector< double >::difference_type,std::vector< double >::difference_type,std::vector< double,std::allocator< double > > const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorDouble___delslice__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::difference_type arg2 ;
  std::vector< double >::difference_type arg3 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  ptrdiff_t val3 ;
  int ecode3 = 0 ;
  PyObject *swig_obj[3] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorDouble___delslice__", 3, 3, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___delslice__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble___delslice__" "', argument " "2"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< double >::difference_type >(val2);
  ecode3 = SWIG_AsVal_ptrdiff_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble___delslice__" "', argument " "3"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg3 = static_cast< std::vector< double >::difference_type >(val3);
  try {
    std_vector_Sl_double_Sg____delslice__(arg1,arg2,arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___delitem____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::difference_type arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___delitem__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble___delitem__" "', argument " "2"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< double >::difference_type >(val2);
  try {
    std_vector_Sl_double_Sg____delitem____SWIG_0(arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___getitem____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  PySliceObject *arg2 = (PySliceObject *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  std::vector< double,std::allocator< double > > *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___getitem__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  {
    if (!PySlice_Check(swig_obj[1])) {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble___getitem__" "', argument " "2"" of type '" "PySliceObject *""'");
    }
    arg2 = (PySliceObject *) swig_obj[1];
  }
  try {
    result = (std::vector< double,std::allocator< double > > *)std_vector_Sl_double_Sg____getitem____SWIG_0(arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___setitem____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  PySliceObject *arg2 = (PySliceObject *) 0 ;
  std::vector< double,std::allocator< double > > *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  int res3 = SWIG_OLDOBJ ;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___setitem__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  {
    if (!PySlice_Check(swig_obj[1])) {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble___setitem__" "', argument " "2"" of type '" "PySliceObject *""'");
    }
    arg2 = (PySliceObject *) swig_obj[1];
  }
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res3 = swig::asptr(swig_obj[2], &ptr);
    if (!SWIG_IsOK(res3)) {
      SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "VectorDouble___setitem__" "', argument " "3"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorDouble___setitem__" "', argument " "3"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg3 = ptr;
  }
  try {
    std_vector_Sl_double_Sg____setitem____SWIG_0(arg1,arg2,(std::vector< double,std::allocator< double > > const &)*arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res3)) delete arg3;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res3)) delete arg3;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___setitem____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  PySliceObject *arg2 = (PySliceObject *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___setitem__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  {
    if (!PySlice_Check(swig_obj[1])) {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble___setitem__" "', argument " "2"" of type '" "PySliceObject *""'");
    }
    arg2 = (PySliceObject *) swig_obj[1];
  }
  try {
    std_vector_Sl_double_Sg____setitem____SWIG_1(arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___delitem____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  PySliceObject *arg2 = (PySliceObject *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___delitem__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  {
    if (!PySlice_Check(swig_obj[1])) {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble___delitem__" "', argument " "2"" of type '" "PySliceObject *""'");
    }
    arg2 = (PySliceObject *) swig_obj[1];
  }
  try {
    std_vector_Sl_double_Sg____delitem____SWIG_1(arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___delitem__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorDouble___delitem__", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        _v = PySlice_Check(argv[1]);
      }
      if (_v) {
        return _wrap_VectorDouble___delitem____SWIG_1(self, argc, argv);
      }
    }
  }
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_VectorDouble___delitem____SWIG_0(self, argc, argv);
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorDouble___delitem__'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< double >::__delitem__(std::vector< double >::difference_type)\n"
    "    std::vector< double >::__delitem__(PySliceObject *)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorDouble___getitem____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::difference_type arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  std::vector< double >::value_type *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___getitem__" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble___getitem__" "', argument " "2"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< double >::difference_type >(val2);
  try {
    result = (std::vector< double >::value_type *) &std_vector_Sl_double_Sg____getitem____SWIG_1((std::vector< double > const *)arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  }
  resultobj = SWIG_From_double(static_cast< double >(*result));
  (void)swig::container_owner<swig::traits<std::vector< double >::value_type>::category>::back_reference(resultobj, swig_obj[0]);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___getitem__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorDouble___getitem__", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        _v = PySlice_Check(argv[1]);
      }
      if (_v) {
        return _wrap_VectorDouble___getitem____SWIG_0(self, argc, argv);
      }
    }
  }
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_VectorDouble___getitem____SWIG_1(self, argc, argv);
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorDouble___getitem__'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< double >::__getitem__(PySliceObject *)\n"
    "    std::vector< double >::__getitem__(std::vector< double >::difference_type) const\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorDouble___setitem____SWIG_2(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::difference_type arg2 ;
  std::vector< double >::value_type *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  std::vector< double >::value_type temp3 ;
  double val3 ;
  int ecode3 = 0 ;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble___setitem__" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble___setitem__" "', argument " "2"" of type '" "std::vector< double >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< double >::difference_type >(val2);
  ecode3 = SWIG_AsVal_double(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble___setitem__" "', argument " "3"" of type '" "std::vector< double >::value_type""'");
  } 
  temp3 = static_cast< std::vector< double >::value_type >(val3);
  arg3 = &temp3;
  try {
    std_vector_Sl_double_Sg____setitem____SWIG_2(arg1,arg2,(double const &)*arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble___setitem__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[4] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorDouble___setitem__", 0, 3, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        _v = PySlice_Check(argv[1]);
      }
      if (_v) {
        return _wrap_VectorDouble___setitem____SWIG_1(self, argc, argv);
      }
    }
  }
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        _v = PySlice_Check(argv[1]);
      }
      if (_v) {
        int res = swig::asptr(argv[2], (std::vector< double,std::allocator< double > >**)(0));
        _v = SWIG_CheckState(res);
        if (_v) {
          return _wrap_VectorDouble___setitem____SWIG_0(self, argc, argv);
        }
      }
    }
  }
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        {
          int res = SWIG_AsVal_double(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          return _wrap_VectorDouble___setitem____SWIG_2(self, argc, argv);
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorDouble___setitem__'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< double >::__setitem__(PySliceObject *,std::vector< double,std::allocator< double > > const &)\n"
    "    std::vector< double >::__setitem__(PySliceObject *)\n"
    "    std::vector< double >::__setitem__(std::vector< double >::difference_type,std::vector< double >::value_type const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorDouble_pop(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::value_type result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_pop" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  try {
    result = (std::vector< double >::value_type)std_vector_Sl_double_Sg__pop(arg1);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  }
  resultobj = SWIG_From_double(static_cast< double >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_append(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::value_type *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  std::vector< double >::value_type temp2 ;
  double val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorDouble_append", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_append" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_double(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble_append" "', argument " "2"" of type '" "std::vector< double >::value_type""'");
  } 
  temp2 = static_cast< std::vector< double >::value_type >(val2);
  arg2 = &temp2;
  std_vector_Sl_double_Sg__append(arg1,(double const &)*arg2);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_new_VectorDouble__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **SWIGUNUSEDPARM(swig_obj)) {
  PyObject *resultobj = 0;
  std::vector< double > *result = 0 ;
  
  if ((nobjs < 0) || (nobjs > 0)) SWIG_fail;
  result = (std::vector< double > *)new std::vector< double >();
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, SWIG_POINTER_NEW |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_new_VectorDouble__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = 0 ;
  int res1 = SWIG_OLDOBJ ;
  std::vector< double > *result = 0 ;
  
  if ((nobjs < 1) || (nobjs > 1)) SWIG_fail;
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res1 = swig::asptr(swig_obj[0], &ptr);
    if (!SWIG_IsOK(res1)) {
      SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "new_VectorDouble" "', argument " "1"" of type '" "std::vector< double > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "new_VectorDouble" "', argument " "1"" of type '" "std::vector< double > const &""'"); 
    }
    arg1 = ptr;
  }
  result = (std::vector< double > *)new std::vector< double >((std::vector< double > const &)*arg1);
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, SWIG_POINTER_NEW |  0 );
  if (SWIG_IsNewObj(res1)) delete arg1;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res1)) delete arg1;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_empty(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  bool result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_empty" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (bool)((std::vector< double > const *)arg1)->empty();
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_size(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::size_type result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_size" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = ((std::vector< double > const *)arg1)->size();
  resultobj = SWIG_From_size_t(static_cast< size_t >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_swap(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double > *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  void *argp2 = 0 ;
  int res2 = 0 ;
  PyObject *swig_obj[2] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorDouble_swap", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_swap" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], &argp2, SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t,  0 );
  if (!SWIG_IsOK(res2)) {
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "VectorDouble_swap" "', argument " "2"" of type '" "std::vector< double > &""'"); 
  }
  if (!argp2) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorDouble_swap" "', argument " "2"" of type '" "std::vector< double > &""'"); 
  }
  arg2 = reinterpret_cast< std::vector< double > * >(argp2);
  (arg1)->swap(*arg2);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_begin(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::iterator result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_begin" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (arg1)->begin();
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< double >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_end(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::iterator result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_end" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (arg1)->end();
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< double >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_rbegin(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::reverse_iterator result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_rbegin" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (arg1)->rbegin();
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< double >::reverse_iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_rend(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::reverse_iterator result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_rend" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (arg1)->rend();
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< double >::reverse_iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_clear(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_clear" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  (arg1)->clear();
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_get_allocator(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  SwigValueWrapper< std::allocator< double > > result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_get_allocator" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = ((std::vector< double > const *)arg1)->get_allocator();
  resultobj = SWIG_NewPointerObj((new std::vector< double >::allocator_type(static_cast< const std::vector< double >::allocator_type& >(result))), SWIGTYPE_p_std__allocatorT_double_t, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_new_VectorDouble__SWIG_2(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double >::size_type arg1 ;
  size_t val1 ;
  int ecode1 = 0 ;
  std::vector< double > *result = 0 ;
  
  if ((nobjs < 1) || (nobjs > 1)) SWIG_fail;
  ecode1 = SWIG_AsVal_size_t(swig_obj[0], &val1);
  if (!SWIG_IsOK(ecode1)) {
    SWIG_exception_fail(SWIG_ArgError(ecode1), "in method '" "new_VectorDouble" "', argument " "1"" of type '" "std::vector< double >::size_type""'");
  } 
  arg1 = static_cast< std::vector< double >::size_type >(val1);
  result = (std::vector< double > *)new std::vector< double >(arg1);
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, SWIG_POINTER_NEW |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_pop_back(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_pop_back" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  (arg1)->pop_back();
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_resize__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::size_type arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_resize" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble_resize" "', argument " "2"" of type '" "std::vector< double >::size_type""'");
  } 
  arg2 = static_cast< std::vector< double >::size_type >(val2);
  (arg1)->resize(arg2);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_erase__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::iterator arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *iter2 = 0 ;
  int res2 ;
  std::vector< double >::iterator result;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_erase" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], SWIG_as_voidptrptr(&iter2), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res2) || !iter2) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_erase" "', argument " "2"" of type '" "std::vector< double >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< double >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter2);
    if (iter_t) {
      arg2 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_erase" "', argument " "2"" of type '" "std::vector< double >::iterator""'");
    }
  }
  result = std_vector_Sl_double_Sg__erase__SWIG_0(arg1,arg2);
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< double >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_erase__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::iterator arg2 ;
  std::vector< double >::iterator arg3 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *iter2 = 0 ;
  int res2 ;
  swig::SwigPyIterator *iter3 = 0 ;
  int res3 ;
  std::vector< double >::iterator result;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_erase" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], SWIG_as_voidptrptr(&iter2), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res2) || !iter2) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_erase" "', argument " "2"" of type '" "std::vector< double >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< double >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter2);
    if (iter_t) {
      arg2 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_erase" "', argument " "2"" of type '" "std::vector< double >::iterator""'");
    }
  }
  res3 = SWIG_ConvertPtr(swig_obj[2], SWIG_as_voidptrptr(&iter3), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res3) || !iter3) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_erase" "', argument " "3"" of type '" "std::vector< double >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< double >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter3);
    if (iter_t) {
      arg3 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_erase" "', argument " "3"" of type '" "std::vector< double >::iterator""'");
    }
  }
  result = std_vector_Sl_double_Sg__erase__SWIG_1(arg1,arg2,arg3);
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< double >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_erase(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[4] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorDouble_erase", 0, 3, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      swig::SwigPyIterator *iter = 0;
      int res = SWIG_ConvertPtr(argv[1], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
      _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter) != 0));
      if (_v) {
        return _wrap_VectorDouble_erase__SWIG_0(self, argc, argv);
      }
    }
  }
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      swig::SwigPyIterator *iter = 0;
      int res = SWIG_ConvertPtr(argv[1], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
      _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter) != 0));
      if (_v) {
        swig::SwigPyIterator *iter = 0;
        int res = SWIG_ConvertPtr(argv[2], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
        _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter) != 0));
        if (_v) {
          return _wrap_VectorDouble_erase__SWIG_1(self, argc, argv);
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorDouble_erase'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< double >::erase(std::vector< double >::iterator)\n"
    "    std::vector< double >::erase(std::vector< double >::iterator,std::vector< double >::iterator)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_new_VectorDouble__SWIG_3(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double >::size_type arg1 ;
  std::vector< double >::value_type *arg2 = 0 ;
  size_t val1 ;
  int ecode1 = 0 ;
  std::vector< double >::value_type temp2 ;
  double val2 ;
  int ecode2 = 0 ;
  std::vector< double > *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  ecode1 = SWIG_AsVal_size_t(swig_obj[0], &val1);
  if (!SWIG_IsOK(ecode1)) {
    SWIG_exception_fail(SWIG_ArgError(ecode1), "in method '" "new_VectorDouble" "', argument " "1"" of type '" "std::vector< double >::size_type""'");
  } 
  arg1 = static_cast< std::vector< double >::size_type >(val1);
  ecode2 = SWIG_AsVal_double(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "new_VectorDouble" "', argument " "2"" of type '" "std::vector< double >::value_type""'");
  } 
  temp2 = static_cast< std::vector< double >::value_type >(val2);
  arg2 = &temp2;
  result = (std::vector< double > *)new std::vector< double >(arg1,(std::vector< double >::value_type const &)*arg2);
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, SWIG_POINTER_NEW |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_new_VectorDouble(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "new_VectorDouble", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 0) {
    return _wrap_new_VectorDouble__SWIG_0(self, argc, argv);
  }
  if (argc == 1) {
    int _v;
    {
      int res = SWIG_AsVal_size_t(argv[0], NULL);
      _v = SWIG_CheckState(res);
    }
    if (_v) {
      return _wrap_new_VectorDouble__SWIG_2(self, argc, argv);
    }
  }
  if (argc == 1) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      return _wrap_new_VectorDouble__SWIG_1(self, argc, argv);
    }
  }
  if (argc == 2) {
    int _v;
    {
      int res = SWIG_AsVal_size_t(argv[0], NULL);
      _v = SWIG_CheckState(res);
    }
    if (_v) {
      {
        int res = SWIG_AsVal_double(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_new_VectorDouble__SWIG_3(self, argc, argv);
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'new_VectorDouble'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< double >::vector()\n"
    "    std::vector< double >::vector(std::vector< double > const &)\n"
    "    std::vector< double >::vector(std::vector< double >::size_type)\n"
    "    std::vector< double >::vector(std::vector< double >::size_type,std::vector< double >::value_type const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorDouble_push_back(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::value_type *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  std::vector< double >::value_type temp2 ;
  double val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorDouble_push_back", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_push_back" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_double(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble_push_back" "', argument " "2"" of type '" "std::vector< double >::value_type""'");
  } 
  temp2 = static_cast< std::vector< double >::value_type >(val2);
  arg2 = &temp2;
  (arg1)->push_back((std::vector< double >::value_type const &)*arg2);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_front(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::value_type *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_front" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (std::vector< double >::value_type *) &((std::vector< double > const *)arg1)->front();
  resultobj = SWIG_From_double(static_cast< double >(*result));
  (void)swig::container_owner<swig::traits<std::vector< double >::value_type>::category>::back_reference(resultobj, swig_obj[0]);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_back(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::value_type *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_back" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = (std::vector< double >::value_type *) &((std::vector< double > const *)arg1)->back();
  resultobj = SWIG_From_double(static_cast< double >(*result));
  (void)swig::container_owner<swig::traits<std::vector< double >::value_type>::category>::back_reference(resultobj, swig_obj[0]);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_assign(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::size_type arg2 ;
  std::vector< double >::value_type *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  std::vector< double >::value_type temp3 ;
  double val3 ;
  int ecode3 = 0 ;
  PyObject *swig_obj[3] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorDouble_assign", 3, 3, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_assign" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble_assign" "', argument " "2"" of type '" "std::vector< double >::size_type""'");
  } 
  arg2 = static_cast< std::vector< double >::size_type >(val2);
  ecode3 = SWIG_AsVal_double(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble_assign" "', argument " "3"" of type '" "std::vector< double >::value_type""'");
  } 
  temp3 = static_cast< std::vector< double >::value_type >(val3);
  arg3 = &temp3;
  (arg1)->assign(arg2,(std::vector< double >::value_type const &)*arg3);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_resize__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::size_type arg2 ;
  std::vector< double >::value_type *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  std::vector< double >::value_type temp3 ;
  double val3 ;
  int ecode3 = 0 ;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_resize" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble_resize" "', argument " "2"" of type '" "std::vector< double >::size_type""'");
  } 
  arg2 = static_cast< std::vector< double >::size_type >(val2);
  ecode3 = SWIG_AsVal_double(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble_resize" "', argument " "3"" of type '" "std::vector< double >::value_type""'");
  } 
  temp3 = static_cast< std::vector< double >::value_type >(val3);
  arg3 = &temp3;
  (arg1)->resize(arg2,(std::vector< double >::value_type const &)*arg3);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_resize(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[4] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorDouble_resize", 0, 3, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_size_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_VectorDouble_resize__SWIG_0(self, argc, argv);
      }
    }
  }
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_size_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        {
          int res = SWIG_AsVal_double(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          return _wrap_VectorDouble_resize__SWIG_1(self, argc, argv);
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorDouble_resize'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< double >::resize(std::vector< double >::size_type)\n"
    "    std::vector< double >::resize(std::vector< double >::size_type,std::vector< double >::value_type const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorDouble_insert__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::iterator arg2 ;
  std::vector< double >::value_type *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *iter2 = 0 ;
  int res2 ;
  std::vector< double >::value_type temp3 ;
  double val3 ;
  int ecode3 = 0 ;
  std::vector< double >::iterator result;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_insert" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], SWIG_as_voidptrptr(&iter2), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res2) || !iter2) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_insert" "', argument " "2"" of type '" "std::vector< double >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< double >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter2);
    if (iter_t) {
      arg2 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_insert" "', argument " "2"" of type '" "std::vector< double >::iterator""'");
    }
  }
  ecode3 = SWIG_AsVal_double(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble_insert" "', argument " "3"" of type '" "std::vector< double >::value_type""'");
  } 
  temp3 = static_cast< std::vector< double >::value_type >(val3);
  arg3 = &temp3;
  result = std_vector_Sl_double_Sg__insert__SWIG_0(arg1,arg2,(double const &)*arg3);
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< double >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_insert__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::iterator arg2 ;
  std::vector< double >::size_type arg3 ;
  std::vector< double >::value_type *arg4 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *iter2 = 0 ;
  int res2 ;
  size_t val3 ;
  int ecode3 = 0 ;
  std::vector< double >::value_type temp4 ;
  double val4 ;
  int ecode4 = 0 ;
  
  if ((nobjs < 4) || (nobjs > 4)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_insert" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], SWIG_as_voidptrptr(&iter2), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res2) || !iter2) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_insert" "', argument " "2"" of type '" "std::vector< double >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< double >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter2);
    if (iter_t) {
      arg2 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorDouble_insert" "', argument " "2"" of type '" "std::vector< double >::iterator""'");
    }
  }
  ecode3 = SWIG_AsVal_size_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorDouble_insert" "', argument " "3"" of type '" "std::vector< double >::size_type""'");
  } 
  arg3 = static_cast< std::vector< double >::size_type >(val3);
  ecode4 = SWIG_AsVal_double(swig_obj[3], &val4);
  if (!SWIG_IsOK(ecode4)) {
    SWIG_exception_fail(SWIG_ArgError(ecode4), "in method '" "VectorDouble_insert" "', argument " "4"" of type '" "std::vector< double >::value_type""'");
  } 
  temp4 = static_cast< std::vector< double >::value_type >(val4);
  arg4 = &temp4;
  std_vector_Sl_double_Sg__insert__SWIG_1(arg1,arg2,arg3,(double const &)*arg4);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_insert(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[5] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorDouble_insert", 0, 4, argv))) SWIG_fail;
  --argc;
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      swig::SwigPyIterator *iter = 0;
      int res = SWIG_ConvertPtr(argv[1], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
      _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter) != 0));
      if (_v) {
        {
          int res = SWIG_AsVal_double(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          return _wrap_VectorDouble_insert__SWIG_0(self, argc, argv);
        }
      }
    }
  }
  if (argc == 4) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< double,std::allocator< double > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      swig::SwigPyIterator *iter = 0;
      int res = SWIG_ConvertPtr(argv[1], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
      _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< double >::iterator > *>(iter) != 0));
      if (_v) {
        {
          int res = SWIG_AsVal_size_t(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          {
            int res = SWIG_AsVal_double(argv[3], NULL);
            _v = SWIG_CheckState(res);
          }
          if (_v) {
            return _wrap_VectorDouble_insert__SWIG_1(self, argc, argv);
          }
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorDouble_insert'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< double >::insert(std::vector< double >::iterator,std::vector< double >::value_type const &)\n"
    "    std::vector< double >::insert(std::vector< double >::iterator,std::vector< double >::size_type,std::vector< double >::value_type const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorDouble_reserve(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  std::vector< double >::size_type arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorDouble_reserve", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_reserve" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorDouble_reserve" "', argument " "2"" of type '" "std::vector< double >::size_type""'");
  } 
  arg2 = static_cast< std::vector< double >::size_type >(val2);
  (arg1)->reserve(arg2);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorDouble_capacity(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< double >::size_type result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorDouble_capacity" "', argument " "1"" of type '" "std::vector< double > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  result = ((std::vector< double > const *)arg1)->capacity();
  resultobj = SWIG_From_size_t(static_cast< size_t >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_delete_VectorDouble(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double > *arg1 = (std::vector< double > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, SWIG_POINTER_DISOWN |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "delete_VectorDouble" "', argument " "1"" of type '" "std::vector< double > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< double > * >(argp1);
  delete arg1;
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *VectorDouble_swigregister(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *obj;
  if (!SWIG_Python_UnpackTuple(args, "swigregister", 1, 1, &obj)) return NULL;
  SWIG_TypeNewClientData(SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t, SWIG_NewClientData(obj));
  return SWIG_Py_Void();
}

SWIGINTERN PyObject *VectorDouble_swiginit(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  return SWIG_Python_InitShadowInstance(args);
}

SWIGINTERN PyObject *_wrap_VectorString_iterator(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  PyObject **arg2 = (PyObject **) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  swig::SwigPyIterator *result = 0 ;
  
  arg2 = &swig_obj[0];
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_iterator" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (swig::SwigPyIterator *)std_vector_Sl_std_string_Sg__iterator(arg1,arg2);
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_swig__SwigPyIterator, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___nonzero__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  bool result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___nonzero__" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (bool)std_vector_Sl_std_string_Sg____nonzero__((std::vector< std::string > const *)arg1);
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___bool__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  bool result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___bool__" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (bool)std_vector_Sl_std_string_Sg____bool__((std::vector< std::string > const *)arg1);
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___len__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::size_type result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___len__" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = std_vector_Sl_std_string_Sg____len__((std::vector< std::string > const *)arg1);
  resultobj = SWIG_From_size_t(static_cast< size_t >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___getslice__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::difference_type arg2 ;
  std::vector< std::string >::difference_type arg3 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  ptrdiff_t val3 ;
  int ecode3 = 0 ;
  PyObject *swig_obj[3] ;
  std::vector< std::string,std::allocator< std::string > > *result = 0 ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorString___getslice__", 3, 3, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___getslice__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString___getslice__" "', argument " "2"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::difference_type >(val2);
  ecode3 = SWIG_AsVal_ptrdiff_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorString___getslice__" "', argument " "3"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg3 = static_cast< std::vector< std::string >::difference_type >(val3);
  try {
    result = (std::vector< std::string,std::allocator< std::string > > *)std_vector_Sl_std_string_Sg____getslice__(arg1,arg2,arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___setslice____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::difference_type arg2 ;
  std::vector< std::string >::difference_type arg3 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  ptrdiff_t val3 ;
  int ecode3 = 0 ;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___setslice__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString___setslice__" "', argument " "2"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::difference_type >(val2);
  ecode3 = SWIG_AsVal_ptrdiff_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorString___setslice__" "', argument " "3"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg3 = static_cast< std::vector< std::string >::difference_type >(val3);
  try {
    std_vector_Sl_std_string_Sg____setslice____SWIG_0(arg1,arg2,arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___setslice____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::difference_type arg2 ;
  std::vector< std::string >::difference_type arg3 ;
  std::vector< std::string,std::allocator< std::string > > *arg4 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  ptrdiff_t val3 ;
  int ecode3 = 0 ;
  int res4 = SWIG_OLDOBJ ;
  
  if ((nobjs < 4) || (nobjs > 4)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___setslice__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString___setslice__" "', argument " "2"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::difference_type >(val2);
  ecode3 = SWIG_AsVal_ptrdiff_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorString___setslice__" "', argument " "3"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg3 = static_cast< std::vector< std::string >::difference_type >(val3);
  {
    std::vector< std::string,std::allocator< std::string > > *ptr = (std::vector< std::string,std::allocator< std::string > > *)0;
    res4 = swig::asptr(swig_obj[3], &ptr);
    if (!SWIG_IsOK(res4)) {
      SWIG_exception_fail(SWIG_ArgError(res4), "in method '" "VectorString___setslice__" "', argument " "4"" of type '" "std::vector< std::string,std::allocator< std::string > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString___setslice__" "', argument " "4"" of type '" "std::vector< std::string,std::allocator< std::string > > const &""'"); 
    }
    arg4 = ptr;
  }
  try {
    std_vector_Sl_std_string_Sg____setslice____SWIG_1(arg1,arg2,arg3,(std::vector< std::string,std::allocator< std::string > > const &)*arg4);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res4)) delete arg4;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res4)) delete arg4;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___setslice__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[5] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorString___setslice__", 0, 4, argv))) SWIG_fail;
  --argc;
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        {
          int res = SWIG_AsVal_ptrdiff_t(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          return _wrap_VectorString___setslice____SWIG_0(self, argc, argv);
        }
      }
    }
  }
  if (argc == 4) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        {
          int res = SWIG_AsVal_ptrdiff_t(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          int res = swig::asptr(argv[3], (std::vector< std::string,std::allocator< std::string > >**)(0));
          _v = SWIG_CheckState(res);
          if (_v) {
            return _wrap_VectorString___setslice____SWIG_1(self, argc, argv);
          }
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorString___setslice__'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< std::string >::__setslice__(std::vector< std::string >::difference_type,std::vector< std::string >::difference_type)\n"
    "    std::vector< std::string >::__setslice__(std::vector< std::string >::difference_type,std::vector< std::string >::difference_type,std::vector< std::string,std::allocator< std::string > > const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorString___delslice__(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::difference_type arg2 ;
  std::vector< std::string >::difference_type arg3 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  ptrdiff_t val3 ;
  int ecode3 = 0 ;
  PyObject *swig_obj[3] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorString___delslice__", 3, 3, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___delslice__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString___delslice__" "', argument " "2"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::difference_type >(val2);
  ecode3 = SWIG_AsVal_ptrdiff_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorString___delslice__" "', argument " "3"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg3 = static_cast< std::vector< std::string >::difference_type >(val3);
  try {
    std_vector_Sl_std_string_Sg____delslice__(arg1,arg2,arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___delitem____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::difference_type arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___delitem__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString___delitem__" "', argument " "2"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::difference_type >(val2);
  try {
    std_vector_Sl_std_string_Sg____delitem____SWIG_0(arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___getitem____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  PySliceObject *arg2 = (PySliceObject *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  std::vector< std::string,std::allocator< std::string > > *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___getitem__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  {
    if (!PySlice_Check(swig_obj[1])) {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString___getitem__" "', argument " "2"" of type '" "PySliceObject *""'");
    }
    arg2 = (PySliceObject *) swig_obj[1];
  }
  try {
    result = (std::vector< std::string,std::allocator< std::string > > *)std_vector_Sl_std_string_Sg____getitem____SWIG_0(arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___setitem____SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  PySliceObject *arg2 = (PySliceObject *) 0 ;
  std::vector< std::string,std::allocator< std::string > > *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  int res3 = SWIG_OLDOBJ ;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___setitem__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  {
    if (!PySlice_Check(swig_obj[1])) {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString___setitem__" "', argument " "2"" of type '" "PySliceObject *""'");
    }
    arg2 = (PySliceObject *) swig_obj[1];
  }
  {
    std::vector< std::string,std::allocator< std::string > > *ptr = (std::vector< std::string,std::allocator< std::string > > *)0;
    res3 = swig::asptr(swig_obj[2], &ptr);
    if (!SWIG_IsOK(res3)) {
      SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "VectorString___setitem__" "', argument " "3"" of type '" "std::vector< std::string,std::allocator< std::string > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString___setitem__" "', argument " "3"" of type '" "std::vector< std::string,std::allocator< std::string > > const &""'"); 
    }
    arg3 = ptr;
  }
  try {
    std_vector_Sl_std_string_Sg____setitem____SWIG_0(arg1,arg2,(std::vector< std::string,std::allocator< std::string > > const &)*arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res3)) delete arg3;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res3)) delete arg3;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___setitem____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  PySliceObject *arg2 = (PySliceObject *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___setitem__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  {
    if (!PySlice_Check(swig_obj[1])) {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString___setitem__" "', argument " "2"" of type '" "PySliceObject *""'");
    }
    arg2 = (PySliceObject *) swig_obj[1];
  }
  try {
    std_vector_Sl_std_string_Sg____setitem____SWIG_1(arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___delitem____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  PySliceObject *arg2 = (PySliceObject *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___delitem__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  {
    if (!PySlice_Check(swig_obj[1])) {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString___delitem__" "', argument " "2"" of type '" "PySliceObject *""'");
    }
    arg2 = (PySliceObject *) swig_obj[1];
  }
  try {
    std_vector_Sl_std_string_Sg____delitem____SWIG_1(arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  } catch(std::invalid_argument &_e) {
    SWIG_exception_fail(SWIG_ValueError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___delitem__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorString___delitem__", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        _v = PySlice_Check(argv[1]);
      }
      if (_v) {
        return _wrap_VectorString___delitem____SWIG_1(self, argc, argv);
      }
    }
  }
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_VectorString___delitem____SWIG_0(self, argc, argv);
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorString___delitem__'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< std::string >::__delitem__(std::vector< std::string >::difference_type)\n"
    "    std::vector< std::string >::__delitem__(PySliceObject *)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorString___getitem____SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::difference_type arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  std::vector< std::string >::value_type *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___getitem__" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString___getitem__" "', argument " "2"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::difference_type >(val2);
  try {
    result = (std::vector< std::string >::value_type *) &std_vector_Sl_std_string_Sg____getitem____SWIG_1((std::vector< std::string > const *)arg1,arg2);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  }
  resultobj = SWIG_From_std_string(static_cast< std::string >(*result));
  (void)swig::container_owner<swig::traits<std::vector< std::string >::value_type>::category>::back_reference(resultobj, swig_obj[0]);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___getitem__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorString___getitem__", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        _v = PySlice_Check(argv[1]);
      }
      if (_v) {
        return _wrap_VectorString___getitem____SWIG_0(self, argc, argv);
      }
    }
  }
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_VectorString___getitem____SWIG_1(self, argc, argv);
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorString___getitem__'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< std::string >::__getitem__(PySliceObject *)\n"
    "    std::vector< std::string >::__getitem__(std::vector< std::string >::difference_type) const\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorString___setitem____SWIG_2(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::difference_type arg2 ;
  std::vector< std::string >::value_type *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  ptrdiff_t val2 ;
  int ecode2 = 0 ;
  int res3 = SWIG_OLDOBJ ;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString___setitem__" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_ptrdiff_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString___setitem__" "', argument " "2"" of type '" "std::vector< std::string >::difference_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::difference_type >(val2);
  {
    std::string *ptr = (std::string *)0;
    res3 = SWIG_AsPtr_std_string(swig_obj[2], &ptr);
    if (!SWIG_IsOK(res3)) {
      SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "VectorString___setitem__" "', argument " "3"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString___setitem__" "', argument " "3"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    arg3 = ptr;
  }
  try {
    std_vector_Sl_std_string_Sg____setitem____SWIG_2(arg1,arg2,(std::string const &)*arg3);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  }
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res3)) delete arg3;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res3)) delete arg3;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString___setitem__(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[4] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorString___setitem__", 0, 3, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        _v = PySlice_Check(argv[1]);
      }
      if (_v) {
        return _wrap_VectorString___setitem____SWIG_1(self, argc, argv);
      }
    }
  }
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        _v = PySlice_Check(argv[1]);
      }
      if (_v) {
        int res = swig::asptr(argv[2], (std::vector< std::string,std::allocator< std::string > >**)(0));
        _v = SWIG_CheckState(res);
        if (_v) {
          return _wrap_VectorString___setitem____SWIG_0(self, argc, argv);
        }
      }
    }
  }
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_ptrdiff_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        int res = SWIG_AsPtr_std_string(argv[2], (std::string**)(0));
        _v = SWIG_CheckState(res);
        if (_v) {
          return _wrap_VectorString___setitem____SWIG_2(self, argc, argv);
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorString___setitem__'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< std::string >::__setitem__(PySliceObject *,std::vector< std::string,std::allocator< std::string > > const &)\n"
    "    std::vector< std::string >::__setitem__(PySliceObject *)\n"
    "    std::vector< std::string >::__setitem__(std::vector< std::string >::difference_type,std::vector< std::string >::value_type const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorString_pop(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::value_type result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_pop" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  try {
    result = std_vector_Sl_std_string_Sg__pop(arg1);
  } catch(std::out_of_range &_e) {
    SWIG_exception_fail(SWIG_IndexError, (&_e)->what());
  }
  resultobj = SWIG_From_std_string(static_cast< std::string >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_append(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::value_type *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  int res2 = SWIG_OLDOBJ ;
  PyObject *swig_obj[2] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorString_append", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_append" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  {
    std::string *ptr = (std::string *)0;
    res2 = SWIG_AsPtr_std_string(swig_obj[1], &ptr);
    if (!SWIG_IsOK(res2)) {
      SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "VectorString_append" "', argument " "2"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString_append" "', argument " "2"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    arg2 = ptr;
  }
  std_vector_Sl_std_string_Sg__append(arg1,(std::string const &)*arg2);
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res2)) delete arg2;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res2)) delete arg2;
  return NULL;
}


SWIGINTERN PyObject *_wrap_new_VectorString__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **SWIGUNUSEDPARM(swig_obj)) {
  PyObject *resultobj = 0;
  std::vector< std::string > *result = 0 ;
  
  if ((nobjs < 0) || (nobjs > 0)) SWIG_fail;
  result = (std::vector< std::string > *)new std::vector< std::string >();
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, SWIG_POINTER_NEW |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_new_VectorString__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = 0 ;
  int res1 = SWIG_OLDOBJ ;
  std::vector< std::string > *result = 0 ;
  
  if ((nobjs < 1) || (nobjs > 1)) SWIG_fail;
  {
    std::vector< std::string,std::allocator< std::string > > *ptr = (std::vector< std::string,std::allocator< std::string > > *)0;
    res1 = swig::asptr(swig_obj[0], &ptr);
    if (!SWIG_IsOK(res1)) {
      SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "new_VectorString" "', argument " "1"" of type '" "std::vector< std::string > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "new_VectorString" "', argument " "1"" of type '" "std::vector< std::string > const &""'"); 
    }
    arg1 = ptr;
  }
  result = (std::vector< std::string > *)new std::vector< std::string >((std::vector< std::string > const &)*arg1);
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, SWIG_POINTER_NEW |  0 );
  if (SWIG_IsNewObj(res1)) delete arg1;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res1)) delete arg1;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_empty(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  bool result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_empty" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (bool)((std::vector< std::string > const *)arg1)->empty();
  resultobj = SWIG_From_bool(static_cast< bool >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_size(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::size_type result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_size" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = ((std::vector< std::string > const *)arg1)->size();
  resultobj = SWIG_From_size_t(static_cast< size_t >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_swap(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string > *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  void *argp2 = 0 ;
  int res2 = 0 ;
  PyObject *swig_obj[2] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorString_swap", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_swap" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], &argp2, SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t,  0 );
  if (!SWIG_IsOK(res2)) {
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "VectorString_swap" "', argument " "2"" of type '" "std::vector< std::string > &""'"); 
  }
  if (!argp2) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString_swap" "', argument " "2"" of type '" "std::vector< std::string > &""'"); 
  }
  arg2 = reinterpret_cast< std::vector< std::string > * >(argp2);
  (arg1)->swap(*arg2);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_begin(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::iterator result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_begin" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (arg1)->begin();
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< std::string >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_end(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::iterator result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_end" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (arg1)->end();
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< std::string >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_rbegin(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::reverse_iterator result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_rbegin" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (arg1)->rbegin();
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< std::string >::reverse_iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_rend(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::reverse_iterator result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_rend" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (arg1)->rend();
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< std::string >::reverse_iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_clear(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_clear" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  (arg1)->clear();
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_get_allocator(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  SwigValueWrapper< std::allocator< std::string > > result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_get_allocator" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = ((std::vector< std::string > const *)arg1)->get_allocator();
  resultobj = SWIG_NewPointerObj((new std::vector< std::string >::allocator_type(static_cast< const std::vector< std::string >::allocator_type& >(result))), SWIGTYPE_p_std__allocatorT_std__string_t, SWIG_POINTER_OWN |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_new_VectorString__SWIG_2(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string >::size_type arg1 ;
  size_t val1 ;
  int ecode1 = 0 ;
  std::vector< std::string > *result = 0 ;
  
  if ((nobjs < 1) || (nobjs > 1)) SWIG_fail;
  ecode1 = SWIG_AsVal_size_t(swig_obj[0], &val1);
  if (!SWIG_IsOK(ecode1)) {
    SWIG_exception_fail(SWIG_ArgError(ecode1), "in method '" "new_VectorString" "', argument " "1"" of type '" "std::vector< std::string >::size_type""'");
  } 
  arg1 = static_cast< std::vector< std::string >::size_type >(val1);
  result = (std::vector< std::string > *)new std::vector< std::string >(arg1);
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, SWIG_POINTER_NEW |  0 );
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_pop_back(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_pop_back" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  (arg1)->pop_back();
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_resize__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::size_type arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_resize" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString_resize" "', argument " "2"" of type '" "std::vector< std::string >::size_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::size_type >(val2);
  (arg1)->resize(arg2);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_erase__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::iterator arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *iter2 = 0 ;
  int res2 ;
  std::vector< std::string >::iterator result;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_erase" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], SWIG_as_voidptrptr(&iter2), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res2) || !iter2) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_erase" "', argument " "2"" of type '" "std::vector< std::string >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< std::string >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter2);
    if (iter_t) {
      arg2 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_erase" "', argument " "2"" of type '" "std::vector< std::string >::iterator""'");
    }
  }
  result = std_vector_Sl_std_string_Sg__erase__SWIG_0(arg1,arg2);
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< std::string >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_erase__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::iterator arg2 ;
  std::vector< std::string >::iterator arg3 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *iter2 = 0 ;
  int res2 ;
  swig::SwigPyIterator *iter3 = 0 ;
  int res3 ;
  std::vector< std::string >::iterator result;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_erase" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], SWIG_as_voidptrptr(&iter2), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res2) || !iter2) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_erase" "', argument " "2"" of type '" "std::vector< std::string >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< std::string >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter2);
    if (iter_t) {
      arg2 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_erase" "', argument " "2"" of type '" "std::vector< std::string >::iterator""'");
    }
  }
  res3 = SWIG_ConvertPtr(swig_obj[2], SWIG_as_voidptrptr(&iter3), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res3) || !iter3) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_erase" "', argument " "3"" of type '" "std::vector< std::string >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< std::string >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter3);
    if (iter_t) {
      arg3 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_erase" "', argument " "3"" of type '" "std::vector< std::string >::iterator""'");
    }
  }
  result = std_vector_Sl_std_string_Sg__erase__SWIG_1(arg1,arg2,arg3);
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< std::string >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_erase(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[4] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorString_erase", 0, 3, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      swig::SwigPyIterator *iter = 0;
      int res = SWIG_ConvertPtr(argv[1], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
      _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter) != 0));
      if (_v) {
        return _wrap_VectorString_erase__SWIG_0(self, argc, argv);
      }
    }
  }
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      swig::SwigPyIterator *iter = 0;
      int res = SWIG_ConvertPtr(argv[1], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
      _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter) != 0));
      if (_v) {
        swig::SwigPyIterator *iter = 0;
        int res = SWIG_ConvertPtr(argv[2], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
        _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter) != 0));
        if (_v) {
          return _wrap_VectorString_erase__SWIG_1(self, argc, argv);
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorString_erase'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< std::string >::erase(std::vector< std::string >::iterator)\n"
    "    std::vector< std::string >::erase(std::vector< std::string >::iterator,std::vector< std::string >::iterator)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_new_VectorString__SWIG_3(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string >::size_type arg1 ;
  std::vector< std::string >::value_type *arg2 = 0 ;
  size_t val1 ;
  int ecode1 = 0 ;
  int res2 = SWIG_OLDOBJ ;
  std::vector< std::string > *result = 0 ;
  
  if ((nobjs < 2) || (nobjs > 2)) SWIG_fail;
  ecode1 = SWIG_AsVal_size_t(swig_obj[0], &val1);
  if (!SWIG_IsOK(ecode1)) {
    SWIG_exception_fail(SWIG_ArgError(ecode1), "in method '" "new_VectorString" "', argument " "1"" of type '" "std::vector< std::string >::size_type""'");
  } 
  arg1 = static_cast< std::vector< std::string >::size_type >(val1);
  {
    std::string *ptr = (std::string *)0;
    res2 = SWIG_AsPtr_std_string(swig_obj[1], &ptr);
    if (!SWIG_IsOK(res2)) {
      SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "new_VectorString" "', argument " "2"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "new_VectorString" "', argument " "2"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    arg2 = ptr;
  }
  result = (std::vector< std::string > *)new std::vector< std::string >(arg1,(std::vector< std::string >::value_type const &)*arg2);
  resultobj = SWIG_NewPointerObj(SWIG_as_voidptr(result), SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, SWIG_POINTER_NEW |  0 );
  if (SWIG_IsNewObj(res2)) delete arg2;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res2)) delete arg2;
  return NULL;
}


SWIGINTERN PyObject *_wrap_new_VectorString(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[3] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "new_VectorString", 0, 2, argv))) SWIG_fail;
  --argc;
  if (argc == 0) {
    return _wrap_new_VectorString__SWIG_0(self, argc, argv);
  }
  if (argc == 1) {
    int _v;
    {
      int res = SWIG_AsVal_size_t(argv[0], NULL);
      _v = SWIG_CheckState(res);
    }
    if (_v) {
      return _wrap_new_VectorString__SWIG_2(self, argc, argv);
    }
  }
  if (argc == 1) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      return _wrap_new_VectorString__SWIG_1(self, argc, argv);
    }
  }
  if (argc == 2) {
    int _v;
    {
      int res = SWIG_AsVal_size_t(argv[0], NULL);
      _v = SWIG_CheckState(res);
    }
    if (_v) {
      int res = SWIG_AsPtr_std_string(argv[1], (std::string**)(0));
      _v = SWIG_CheckState(res);
      if (_v) {
        return _wrap_new_VectorString__SWIG_3(self, argc, argv);
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'new_VectorString'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< std::string >::vector()\n"
    "    std::vector< std::string >::vector(std::vector< std::string > const &)\n"
    "    std::vector< std::string >::vector(std::vector< std::string >::size_type)\n"
    "    std::vector< std::string >::vector(std::vector< std::string >::size_type,std::vector< std::string >::value_type const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorString_push_back(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::value_type *arg2 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  int res2 = SWIG_OLDOBJ ;
  PyObject *swig_obj[2] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorString_push_back", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_push_back" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  {
    std::string *ptr = (std::string *)0;
    res2 = SWIG_AsPtr_std_string(swig_obj[1], &ptr);
    if (!SWIG_IsOK(res2)) {
      SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "VectorString_push_back" "', argument " "2"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString_push_back" "', argument " "2"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    arg2 = ptr;
  }
  (arg1)->push_back((std::vector< std::string >::value_type const &)*arg2);
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res2)) delete arg2;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res2)) delete arg2;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_front(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::value_type *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_front" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (std::vector< std::string >::value_type *) &((std::vector< std::string > const *)arg1)->front();
  resultobj = SWIG_From_std_string(static_cast< std::string >(*result));
  (void)swig::container_owner<swig::traits<std::vector< std::string >::value_type>::category>::back_reference(resultobj, swig_obj[0]);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_back(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::value_type *result = 0 ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_back" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = (std::vector< std::string >::value_type *) &((std::vector< std::string > const *)arg1)->back();
  resultobj = SWIG_From_std_string(static_cast< std::string >(*result));
  (void)swig::container_owner<swig::traits<std::vector< std::string >::value_type>::category>::back_reference(resultobj, swig_obj[0]);
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_assign(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::size_type arg2 ;
  std::vector< std::string >::value_type *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  int res3 = SWIG_OLDOBJ ;
  PyObject *swig_obj[3] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorString_assign", 3, 3, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_assign" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString_assign" "', argument " "2"" of type '" "std::vector< std::string >::size_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::size_type >(val2);
  {
    std::string *ptr = (std::string *)0;
    res3 = SWIG_AsPtr_std_string(swig_obj[2], &ptr);
    if (!SWIG_IsOK(res3)) {
      SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "VectorString_assign" "', argument " "3"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString_assign" "', argument " "3"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    arg3 = ptr;
  }
  (arg1)->assign(arg2,(std::vector< std::string >::value_type const &)*arg3);
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res3)) delete arg3;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res3)) delete arg3;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_resize__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::size_type arg2 ;
  std::vector< std::string >::value_type *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  int res3 = SWIG_OLDOBJ ;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_resize" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString_resize" "', argument " "2"" of type '" "std::vector< std::string >::size_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::size_type >(val2);
  {
    std::string *ptr = (std::string *)0;
    res3 = SWIG_AsPtr_std_string(swig_obj[2], &ptr);
    if (!SWIG_IsOK(res3)) {
      SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "VectorString_resize" "', argument " "3"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString_resize" "', argument " "3"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    arg3 = ptr;
  }
  (arg1)->resize(arg2,(std::vector< std::string >::value_type const &)*arg3);
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res3)) delete arg3;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res3)) delete arg3;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_resize(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[4] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorString_resize", 0, 3, argv))) SWIG_fail;
  --argc;
  if (argc == 2) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_size_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        return _wrap_VectorString_resize__SWIG_0(self, argc, argv);
      }
    }
  }
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      {
        int res = SWIG_AsVal_size_t(argv[1], NULL);
        _v = SWIG_CheckState(res);
      }
      if (_v) {
        int res = SWIG_AsPtr_std_string(argv[2], (std::string**)(0));
        _v = SWIG_CheckState(res);
        if (_v) {
          return _wrap_VectorString_resize__SWIG_1(self, argc, argv);
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorString_resize'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< std::string >::resize(std::vector< std::string >::size_type)\n"
    "    std::vector< std::string >::resize(std::vector< std::string >::size_type,std::vector< std::string >::value_type const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorString_insert__SWIG_0(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::iterator arg2 ;
  std::vector< std::string >::value_type *arg3 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *iter2 = 0 ;
  int res2 ;
  int res3 = SWIG_OLDOBJ ;
  std::vector< std::string >::iterator result;
  
  if ((nobjs < 3) || (nobjs > 3)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_insert" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], SWIG_as_voidptrptr(&iter2), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res2) || !iter2) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_insert" "', argument " "2"" of type '" "std::vector< std::string >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< std::string >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter2);
    if (iter_t) {
      arg2 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_insert" "', argument " "2"" of type '" "std::vector< std::string >::iterator""'");
    }
  }
  {
    std::string *ptr = (std::string *)0;
    res3 = SWIG_AsPtr_std_string(swig_obj[2], &ptr);
    if (!SWIG_IsOK(res3)) {
      SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "VectorString_insert" "', argument " "3"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString_insert" "', argument " "3"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    arg3 = ptr;
  }
  result = std_vector_Sl_std_string_Sg__insert__SWIG_0(arg1,arg2,(std::string const &)*arg3);
  resultobj = SWIG_NewPointerObj(swig::make_output_iterator(static_cast< const std::vector< std::string >::iterator & >(result)),
    swig::SwigPyIterator::descriptor(),SWIG_POINTER_OWN);
  if (SWIG_IsNewObj(res3)) delete arg3;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res3)) delete arg3;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_insert__SWIG_1(PyObject *SWIGUNUSEDPARM(self), Py_ssize_t nobjs, PyObject **swig_obj) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::iterator arg2 ;
  std::vector< std::string >::size_type arg3 ;
  std::vector< std::string >::value_type *arg4 = 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  swig::SwigPyIterator *iter2 = 0 ;
  int res2 ;
  size_t val3 ;
  int ecode3 = 0 ;
  int res4 = SWIG_OLDOBJ ;
  
  if ((nobjs < 4) || (nobjs > 4)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_insert" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  res2 = SWIG_ConvertPtr(swig_obj[1], SWIG_as_voidptrptr(&iter2), swig::SwigPyIterator::descriptor(), 0);
  if (!SWIG_IsOK(res2) || !iter2) {
    SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_insert" "', argument " "2"" of type '" "std::vector< std::string >::iterator""'");
  } else {
    swig::SwigPyIterator_T<std::vector< std::string >::iterator > *iter_t = dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter2);
    if (iter_t) {
      arg2 = iter_t->get_current();
    } else {
      SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '" "VectorString_insert" "', argument " "2"" of type '" "std::vector< std::string >::iterator""'");
    }
  }
  ecode3 = SWIG_AsVal_size_t(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "VectorString_insert" "', argument " "3"" of type '" "std::vector< std::string >::size_type""'");
  } 
  arg3 = static_cast< std::vector< std::string >::size_type >(val3);
  {
    std::string *ptr = (std::string *)0;
    res4 = SWIG_AsPtr_std_string(swig_obj[3], &ptr);
    if (!SWIG_IsOK(res4)) {
      SWIG_exception_fail(SWIG_ArgError(res4), "in method '" "VectorString_insert" "', argument " "4"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "VectorString_insert" "', argument " "4"" of type '" "std::vector< std::string >::value_type const &""'"); 
    }
    arg4 = ptr;
  }
  std_vector_Sl_std_string_Sg__insert__SWIG_1(arg1,arg2,arg3,(std::string const &)*arg4);
  resultobj = SWIG_Py_Void();
  if (SWIG_IsNewObj(res4)) delete arg4;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res4)) delete arg4;
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_insert(PyObject *self, PyObject *args) {
  Py_ssize_t argc;
  PyObject *argv[5] = {
    0
  };
  
  if (!(argc = SWIG_Python_UnpackTuple(args, "VectorString_insert", 0, 4, argv))) SWIG_fail;
  --argc;
  if (argc == 3) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      swig::SwigPyIterator *iter = 0;
      int res = SWIG_ConvertPtr(argv[1], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
      _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter) != 0));
      if (_v) {
        int res = SWIG_AsPtr_std_string(argv[2], (std::string**)(0));
        _v = SWIG_CheckState(res);
        if (_v) {
          return _wrap_VectorString_insert__SWIG_0(self, argc, argv);
        }
      }
    }
  }
  if (argc == 4) {
    int _v;
    int res = swig::asptr(argv[0], (std::vector< std::string,std::allocator< std::string > >**)(0));
    _v = SWIG_CheckState(res);
    if (_v) {
      swig::SwigPyIterator *iter = 0;
      int res = SWIG_ConvertPtr(argv[1], SWIG_as_voidptrptr(&iter), swig::SwigPyIterator::descriptor(), 0);
      _v = (SWIG_IsOK(res) && iter && (dynamic_cast<swig::SwigPyIterator_T<std::vector< std::string >::iterator > *>(iter) != 0));
      if (_v) {
        {
          int res = SWIG_AsVal_size_t(argv[2], NULL);
          _v = SWIG_CheckState(res);
        }
        if (_v) {
          int res = SWIG_AsPtr_std_string(argv[3], (std::string**)(0));
          _v = SWIG_CheckState(res);
          if (_v) {
            return _wrap_VectorString_insert__SWIG_1(self, argc, argv);
          }
        }
      }
    }
  }
  
fail:
  SWIG_Python_RaiseOrModifyTypeError("Wrong number or type of arguments for overloaded function 'VectorString_insert'.\n"
    "  Possible C/C++ prototypes are:\n"
    "    std::vector< std::string >::insert(std::vector< std::string >::iterator,std::vector< std::string >::value_type const &)\n"
    "    std::vector< std::string >::insert(std::vector< std::string >::iterator,std::vector< std::string >::size_type,std::vector< std::string >::value_type const &)\n");
  return 0;
}


SWIGINTERN PyObject *_wrap_VectorString_reserve(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  std::vector< std::string >::size_type arg2 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  size_t val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  
  if (!SWIG_Python_UnpackTuple(args, "VectorString_reserve", 2, 2, swig_obj)) SWIG_fail;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_reserve" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  ecode2 = SWIG_AsVal_size_t(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "VectorString_reserve" "', argument " "2"" of type '" "std::vector< std::string >::size_type""'");
  } 
  arg2 = static_cast< std::vector< std::string >::size_type >(val2);
  (arg1)->reserve(arg2);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_VectorString_capacity(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  std::vector< std::string >::size_type result;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0 |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "VectorString_capacity" "', argument " "1"" of type '" "std::vector< std::string > const *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  result = ((std::vector< std::string > const *)arg1)->capacity();
  resultobj = SWIG_From_size_t(static_cast< size_t >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_delete_VectorString(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< std::string > *arg1 = (std::vector< std::string > *) 0 ;
  void *argp1 = 0 ;
  int res1 = 0 ;
  PyObject *swig_obj[1] ;
  
  if (!args) SWIG_fail;
  swig_obj[0] = args;
  res1 = SWIG_ConvertPtr(swig_obj[0], &argp1,SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, SWIG_POINTER_DISOWN |  0 );
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "delete_VectorString" "', argument " "1"" of type '" "std::vector< std::string > *""'"); 
  }
  arg1 = reinterpret_cast< std::vector< std::string > * >(argp1);
  delete arg1;
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *VectorString_swigregister(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *obj;
  if (!SWIG_Python_UnpackTuple(args, "swigregister", 1, 1, &obj)) return NULL;
  SWIG_TypeNewClientData(SWIGTYPE_p_std__vectorT_std__string_std__allocatorT_std__string_t_t, SWIG_NewClientData(obj));
  return SWIG_Py_Void();
}

SWIGINTERN PyObject *VectorString_swiginit(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  return SWIG_Python_InitShadowInstance(args);
}

SWIGINTERN PyObject *_wrap_initial_condition_call(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double,std::allocator< double > > *arg1 = 0 ;
  double arg2 ;
  int res1 = SWIG_OLDOBJ ;
  double val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  std::vector< double,std::allocator< double > > result;
  
  if (!SWIG_Python_UnpackTuple(args, "initial_condition_call", 2, 2, swig_obj)) SWIG_fail;
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res1 = swig::asptr(swig_obj[0], &ptr);
    if (!SWIG_IsOK(res1)) {
      SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "initial_condition_call" "', argument " "1"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "initial_condition_call" "', argument " "1"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg1 = ptr;
  }
  ecode2 = SWIG_AsVal_double(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "initial_condition_call" "', argument " "2"" of type '" "double""'");
  } 
  arg2 = static_cast< double >(val2);
  result = initial_condition_call((std::vector< double,std::allocator< double > > const &)*arg1,arg2);
  resultobj = swig::from(static_cast< std::vector< double,std::allocator< double > > >(result));
  if (SWIG_IsNewObj(res1)) delete arg1;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res1)) delete arg1;
  return NULL;
}


SWIGINTERN PyObject *_wrap_initial_condition_put(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double,std::allocator< double > > *arg1 = 0 ;
  double arg2 ;
  int res1 = SWIG_OLDOBJ ;
  double val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  std::vector< double,std::allocator< double > > result;
  
  if (!SWIG_Python_UnpackTuple(args, "initial_condition_put", 2, 2, swig_obj)) SWIG_fail;
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res1 = swig::asptr(swig_obj[0], &ptr);
    if (!SWIG_IsOK(res1)) {
      SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "initial_condition_put" "', argument " "1"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "initial_condition_put" "', argument " "1"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg1 = ptr;
  }
  ecode2 = SWIG_AsVal_double(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "initial_condition_put" "', argument " "2"" of type '" "double""'");
  } 
  arg2 = static_cast< double >(val2);
  result = initial_condition_put((std::vector< double,std::allocator< double > > const &)*arg1,arg2);
  resultobj = swig::from(static_cast< std::vector< double,std::allocator< double > > >(result));
  if (SWIG_IsNewObj(res1)) delete arg1;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res1)) delete arg1;
  return NULL;
}


SWIGINTERN PyObject *_wrap_beta_coefficients(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  double arg1 ;
  int arg2 ;
  double val1 ;
  int ecode1 = 0 ;
  int val2 ;
  int ecode2 = 0 ;
  PyObject *swig_obj[2] ;
  std::vector< double,std::allocator< double > > result;
  
  if (!SWIG_Python_UnpackTuple(args, "beta_coefficients", 2, 2, swig_obj)) SWIG_fail;
  ecode1 = SWIG_AsVal_double(swig_obj[0], &val1);
  if (!SWIG_IsOK(ecode1)) {
    SWIG_exception_fail(SWIG_ArgError(ecode1), "in method '" "beta_coefficients" "', argument " "1"" of type '" "double""'");
  } 
  arg1 = static_cast< double >(val1);
  ecode2 = SWIG_AsVal_int(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "beta_coefficients" "', argument " "2"" of type '" "int""'");
  } 
  arg2 = static_cast< int >(val2);
  result = beta_coefficients(arg1,arg2);
  resultobj = swig::from(static_cast< std::vector< double,std::allocator< double > > >(result));
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_construct_tridiagonal_coeffs(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  int arg1 ;
  double arg2 ;
  double arg3 ;
  double arg4 ;
  double arg5 ;
  std::vector< double,std::allocator< double > > *arg6 = 0 ;
  std::vector< double,std::allocator< double > > *arg7 = 0 ;
  std::vector< double,std::allocator< double > > *arg8 = 0 ;
  int val1 ;
  int ecode1 = 0 ;
  double val2 ;
  int ecode2 = 0 ;
  double val3 ;
  int ecode3 = 0 ;
  double val4 ;
  int ecode4 = 0 ;
  double val5 ;
  int ecode5 = 0 ;
  void *argp6 = 0 ;
  int res6 = 0 ;
  void *argp7 = 0 ;
  int res7 = 0 ;
  void *argp8 = 0 ;
  int res8 = 0 ;
  PyObject *swig_obj[8] ;
  
  if (!SWIG_Python_UnpackTuple(args, "construct_tridiagonal_coeffs", 8, 8, swig_obj)) SWIG_fail;
  ecode1 = SWIG_AsVal_int(swig_obj[0], &val1);
  if (!SWIG_IsOK(ecode1)) {
    SWIG_exception_fail(SWIG_ArgError(ecode1), "in method '" "construct_tridiagonal_coeffs" "', argument " "1"" of type '" "int""'");
  } 
  arg1 = static_cast< int >(val1);
  ecode2 = SWIG_AsVal_double(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "construct_tridiagonal_coeffs" "', argument " "2"" of type '" "double""'");
  } 
  arg2 = static_cast< double >(val2);
  ecode3 = SWIG_AsVal_double(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "construct_tridiagonal_coeffs" "', argument " "3"" of type '" "double""'");
  } 
  arg3 = static_cast< double >(val3);
  ecode4 = SWIG_AsVal_double(swig_obj[3], &val4);
  if (!SWIG_IsOK(ecode4)) {
    SWIG_exception_fail(SWIG_ArgError(ecode4), "in method '" "construct_tridiagonal_coeffs" "', argument " "4"" of type '" "double""'");
  } 
  arg4 = static_cast< double >(val4);
  ecode5 = SWIG_AsVal_double(swig_obj[4], &val5);
  if (!SWIG_IsOK(ecode5)) {
    SWIG_exception_fail(SWIG_ArgError(ecode5), "in method '" "construct_tridiagonal_coeffs" "', argument " "5"" of type '" "double""'");
  } 
  arg5 = static_cast< double >(val5);
  res6 = SWIG_ConvertPtr(swig_obj[5], &argp6, SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t,  0 );
  if (!SWIG_IsOK(res6)) {
    SWIG_exception_fail(SWIG_ArgError(res6), "in method '" "construct_tridiagonal_coeffs" "', argument " "6"" of type '" "std::vector< double,std::allocator< double > > &""'"); 
  }
  if (!argp6) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "construct_tridiagonal_coeffs" "', argument " "6"" of type '" "std::vector< double,std::allocator< double > > &""'"); 
  }
  arg6 = reinterpret_cast< std::vector< double,std::allocator< double > > * >(argp6);
  res7 = SWIG_ConvertPtr(swig_obj[6], &argp7, SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t,  0 );
  if (!SWIG_IsOK(res7)) {
    SWIG_exception_fail(SWIG_ArgError(res7), "in method '" "construct_tridiagonal_coeffs" "', argument " "7"" of type '" "std::vector< double,std::allocator< double > > &""'"); 
  }
  if (!argp7) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "construct_tridiagonal_coeffs" "', argument " "7"" of type '" "std::vector< double,std::allocator< double > > &""'"); 
  }
  arg7 = reinterpret_cast< std::vector< double,std::allocator< double > > * >(argp7);
  res8 = SWIG_ConvertPtr(swig_obj[7], &argp8, SWIGTYPE_p_std__vectorT_double_std__allocatorT_double_t_t,  0 );
  if (!SWIG_IsOK(res8)) {
    SWIG_exception_fail(SWIG_ArgError(res8), "in method '" "construct_tridiagonal_coeffs" "', argument " "8"" of type '" "std::vector< double,std::allocator< double > > &""'"); 
  }
  if (!argp8) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "construct_tridiagonal_coeffs" "', argument " "8"" of type '" "std::vector< double,std::allocator< double > > &""'"); 
  }
  arg8 = reinterpret_cast< std::vector< double,std::allocator< double > > * >(argp8);
  construct_tridiagonal_coeffs(arg1,arg2,arg3,arg4,arg5,*arg6,*arg7,*arg8);
  resultobj = SWIG_Py_Void();
  return resultobj;
fail:
  return NULL;
}


SWIGINTERN PyObject *_wrap_solve_tridiagonal(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double,std::allocator< double > > *arg1 = 0 ;
  std::vector< double,std::allocator< double > > *arg2 = 0 ;
  std::vector< double,std::allocator< double > > *arg3 = 0 ;
  std::vector< double,std::allocator< double > > *arg4 = 0 ;
  int res1 = SWIG_OLDOBJ ;
  int res2 = SWIG_OLDOBJ ;
  int res3 = SWIG_OLDOBJ ;
  int res4 = SWIG_OLDOBJ ;
  PyObject *swig_obj[4] ;
  std::vector< double,std::allocator< double > > result;
  
  if (!SWIG_Python_UnpackTuple(args, "solve_tridiagonal", 4, 4, swig_obj)) SWIG_fail;
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res1 = swig::asptr(swig_obj[0], &ptr);
    if (!SWIG_IsOK(res1)) {
      SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "solve_tridiagonal" "', argument " "1"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "solve_tridiagonal" "', argument " "1"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg1 = ptr;
  }
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res2 = swig::asptr(swig_obj[1], &ptr);
    if (!SWIG_IsOK(res2)) {
      SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "solve_tridiagonal" "', argument " "2"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "solve_tridiagonal" "', argument " "2"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg2 = ptr;
  }
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res3 = swig::asptr(swig_obj[2], &ptr);
    if (!SWIG_IsOK(res3)) {
      SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "solve_tridiagonal" "', argument " "3"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "solve_tridiagonal" "', argument " "3"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg3 = ptr;
  }
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res4 = swig::asptr(swig_obj[3], &ptr);
    if (!SWIG_IsOK(res4)) {
      SWIG_exception_fail(SWIG_ArgError(res4), "in method '" "solve_tridiagonal" "', argument " "4"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "solve_tridiagonal" "', argument " "4"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg4 = ptr;
  }
  result = solve_tridiagonal((std::vector< double,std::allocator< double > > const &)*arg1,(std::vector< double,std::allocator< double > > const &)*arg2,(std::vector< double,std::allocator< double > > const &)*arg3,(std::vector< double,std::allocator< double > > const &)*arg4);
  resultobj = swig::from(static_cast< std::vector< double,std::allocator< double > > >(result));
  if (SWIG_IsNewObj(res1)) delete arg1;
  if (SWIG_IsNewObj(res2)) delete arg2;
  if (SWIG_IsNewObj(res3)) delete arg3;
  if (SWIG_IsNewObj(res4)) delete arg4;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res1)) delete arg1;
  if (SWIG_IsNewObj(res2)) delete arg2;
  if (SWIG_IsNewObj(res3)) delete arg3;
  if (SWIG_IsNewObj(res4)) delete arg4;
  return NULL;
}


SWIGINTERN PyObject *_wrap_fractional_black_scholes_alpha(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  std::vector< double,std::allocator< double > > *arg1 = 0 ;
  double arg2 ;
  double arg3 ;
  double arg4 ;
  double arg5 ;
  int arg6 ;
  int arg7 ;
  double arg8 ;
  double arg9 ;
  std::string *arg10 = 0 ;
  int res1 = SWIG_OLDOBJ ;
  double val2 ;
  int ecode2 = 0 ;
  double val3 ;
  int ecode3 = 0 ;
  double val4 ;
  int ecode4 = 0 ;
  double val5 ;
  int ecode5 = 0 ;
  int val6 ;
  int ecode6 = 0 ;
  int val7 ;
  int ecode7 = 0 ;
  double val8 ;
  int ecode8 = 0 ;
  double val9 ;
  int ecode9 = 0 ;
  int res10 = SWIG_OLDOBJ ;
  PyObject *swig_obj[10] ;
  std::vector< double,std::allocator< double > > result;
  
  if (!SWIG_Python_UnpackTuple(args, "fractional_black_scholes_alpha", 10, 10, swig_obj)) SWIG_fail;
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res1 = swig::asptr(swig_obj[0], &ptr);
    if (!SWIG_IsOK(res1)) {
      SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "fractional_black_scholes_alpha" "', argument " "1"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "fractional_black_scholes_alpha" "', argument " "1"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg1 = ptr;
  }
  ecode2 = SWIG_AsVal_double(swig_obj[1], &val2);
  if (!SWIG_IsOK(ecode2)) {
    SWIG_exception_fail(SWIG_ArgError(ecode2), "in method '" "fractional_black_scholes_alpha" "', argument " "2"" of type '" "double""'");
  } 
  arg2 = static_cast< double >(val2);
  ecode3 = SWIG_AsVal_double(swig_obj[2], &val3);
  if (!SWIG_IsOK(ecode3)) {
    SWIG_exception_fail(SWIG_ArgError(ecode3), "in method '" "fractional_black_scholes_alpha" "', argument " "3"" of type '" "double""'");
  } 
  arg3 = static_cast< double >(val3);
  ecode4 = SWIG_AsVal_double(swig_obj[3], &val4);
  if (!SWIG_IsOK(ecode4)) {
    SWIG_exception_fail(SWIG_ArgError(ecode4), "in method '" "fractional_black_scholes_alpha" "', argument " "4"" of type '" "double""'");
  } 
  arg4 = static_cast< double >(val4);
  ecode5 = SWIG_AsVal_double(swig_obj[4], &val5);
  if (!SWIG_IsOK(ecode5)) {
    SWIG_exception_fail(SWIG_ArgError(ecode5), "in method '" "fractional_black_scholes_alpha" "', argument " "5"" of type '" "double""'");
  } 
  arg5 = static_cast< double >(val5);
  ecode6 = SWIG_AsVal_int(swig_obj[5], &val6);
  if (!SWIG_IsOK(ecode6)) {
    SWIG_exception_fail(SWIG_ArgError(ecode6), "in method '" "fractional_black_scholes_alpha" "', argument " "6"" of type '" "int""'");
  } 
  arg6 = static_cast< int >(val6);
  ecode7 = SWIG_AsVal_int(swig_obj[6], &val7);
  if (!SWIG_IsOK(ecode7)) {
    SWIG_exception_fail(SWIG_ArgError(ecode7), "in method '" "fractional_black_scholes_alpha" "', argument " "7"" of type '" "int""'");
  } 
  arg7 = static_cast< int >(val7);
  ecode8 = SWIG_AsVal_double(swig_obj[7], &val8);
  if (!SWIG_IsOK(ecode8)) {
    SWIG_exception_fail(SWIG_ArgError(ecode8), "in method '" "fractional_black_scholes_alpha" "', argument " "8"" of type '" "double""'");
  } 
  arg8 = static_cast< double >(val8);
  ecode9 = SWIG_AsVal_double(swig_obj[8], &val9);
  if (!SWIG_IsOK(ecode9)) {
    SWIG_exception_fail(SWIG_ArgError(ecode9), "in method '" "fractional_black_scholes_alpha" "', argument " "9"" of type '" "double""'");
  } 
  arg9 = static_cast< double >(val9);
  {
    std::string *ptr = (std::string *)0;
    res10 = SWIG_AsPtr_std_string(swig_obj[9], &ptr);
    if (!SWIG_IsOK(res10)) {
      SWIG_exception_fail(SWIG_ArgError(res10), "in method '" "fractional_black_scholes_alpha" "', argument " "10"" of type '" "std::string const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "fractional_black_scholes_alpha" "', argument " "10"" of type '" "std::string const &""'"); 
    }
    arg10 = ptr;
  }
  result = fractional_black_scholes_alpha((std::vector< double,std::allocator< double > > const &)*arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,(std::string const &)*arg10);
  resultobj = swig::from(static_cast< std::vector< double,std::allocator< double > > >(result));
  if (SWIG_IsNewObj(res1)) delete arg1;
  if (SWIG_IsNewObj(res10)) delete arg10;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res1)) delete arg1;
  if (SWIG_IsNewObj(res10)) delete arg10;
  return NULL;
}


SWIGINTERN PyObject *_wrap_calibration_objective(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  double arg1 ;
  std::vector< std::vector< double,std::allocator< double > >,std::allocator< std::vector< double,std::allocator< double > > > > *arg2 = 0 ;
  std::string *arg3 = 0 ;
  std::vector< double,std::allocator< double > > *arg4 = 0 ;
  double arg5 ;
  double arg6 ;
  double arg7 ;
  int arg8 ;
  int arg9 ;
  double arg10 ;
  double val1 ;
  int ecode1 = 0 ;
  void *argp2 = 0 ;
  int res2 = 0 ;
  int res3 = SWIG_OLDOBJ ;
  int res4 = SWIG_OLDOBJ ;
  double val5 ;
  int ecode5 = 0 ;
  double val6 ;
  int ecode6 = 0 ;
  double val7 ;
  int ecode7 = 0 ;
  int val8 ;
  int ecode8 = 0 ;
  int val9 ;
  int ecode9 = 0 ;
  double val10 ;
  int ecode10 = 0 ;
  PyObject *swig_obj[10] ;
  double result;
  
  if (!SWIG_Python_UnpackTuple(args, "calibration_objective", 10, 10, swig_obj)) SWIG_fail;
  ecode1 = SWIG_AsVal_double(swig_obj[0], &val1);
  if (!SWIG_IsOK(ecode1)) {
    SWIG_exception_fail(SWIG_ArgError(ecode1), "in method '" "calibration_objective" "', argument " "1"" of type '" "double""'");
  } 
  arg1 = static_cast< double >(val1);
  res2 = SWIG_ConvertPtr(swig_obj[1], &argp2, SWIGTYPE_p_std__vectorT_std__vectorT_double_std__allocatorT_double_t_t_std__allocatorT_std__vectorT_double_std__allocatorT_double_t_t_t_t,  0  | 0);
  if (!SWIG_IsOK(res2)) {
    SWIG_exception_fail(SWIG_ArgError(res2), "in method '" "calibration_objective" "', argument " "2"" of type '" "std::vector< std::vector< double,std::allocator< double > >,std::allocator< std::vector< double,std::allocator< double > > > > const &""'"); 
  }
  if (!argp2) {
    SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "calibration_objective" "', argument " "2"" of type '" "std::vector< std::vector< double,std::allocator< double > >,std::allocator< std::vector< double,std::allocator< double > > > > const &""'"); 
  }
  arg2 = reinterpret_cast< std::vector< std::vector< double,std::allocator< double > >,std::allocator< std::vector< double,std::allocator< double > > > > * >(argp2);
  {
    std::string *ptr = (std::string *)0;
    res3 = SWIG_AsPtr_std_string(swig_obj[2], &ptr);
    if (!SWIG_IsOK(res3)) {
      SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "calibration_objective" "', argument " "3"" of type '" "std::string const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "calibration_objective" "', argument " "3"" of type '" "std::string const &""'"); 
    }
    arg3 = ptr;
  }
  {
    std::vector< double,std::allocator< double > > *ptr = (std::vector< double,std::allocator< double > > *)0;
    res4 = swig::asptr(swig_obj[3], &ptr);
    if (!SWIG_IsOK(res4)) {
      SWIG_exception_fail(SWIG_ArgError(res4), "in method '" "calibration_objective" "', argument " "4"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    if (!ptr) {
      SWIG_exception_fail(SWIG_ValueError, "invalid null reference " "in method '" "calibration_objective" "', argument " "4"" of type '" "std::vector< double,std::allocator< double > > const &""'"); 
    }
    arg4 = ptr;
  }
  ecode5 = SWIG_AsVal_double(swig_obj[4], &val5);
  if (!SWIG_IsOK(ecode5)) {
    SWIG_exception_fail(SWIG_ArgError(ecode5), "in method '" "calibration_objective" "', argument " "5"" of type '" "double""'");
  } 
  arg5 = static_cast< double >(val5);
  ecode6 = SWIG_AsVal_double(swig_obj[5], &val6);
  if (!SWIG_IsOK(ecode6)) {
    SWIG_exception_fail(SWIG_ArgError(ecode6), "in method '" "calibration_objective" "', argument " "6"" of type '" "double""'");
  } 
  arg6 = static_cast< double >(val6);
  ecode7 = SWIG_AsVal_double(swig_obj[6], &val7);
  if (!SWIG_IsOK(ecode7)) {
    SWIG_exception_fail(SWIG_ArgError(ecode7), "in method '" "calibration_objective" "', argument " "7"" of type '" "double""'");
  } 
  arg7 = static_cast< double >(val7);
  ecode8 = SWIG_AsVal_int(swig_obj[7], &val8);
  if (!SWIG_IsOK(ecode8)) {
    SWIG_exception_fail(SWIG_ArgError(ecode8), "in method '" "calibration_objective" "', argument " "8"" of type '" "int""'");
  } 
  arg8 = static_cast< int >(val8);
  ecode9 = SWIG_AsVal_int(swig_obj[8], &val9);
  if (!SWIG_IsOK(ecode9)) {
    SWIG_exception_fail(SWIG_ArgError(ecode9), "in method '" "calibration_objective" "', argument " "9"" of type '" "int""'");
  } 
  arg9 = static_cast< int >(val9);
  ecode10 = SWIG_AsVal_double(swig_obj[9], &val10);
  if (!SWIG_IsOK(ecode10)) {
    SWIG_exception_fail(SWIG_ArgError(ecode10), "in method '" "calibration_objective" "', argument " "10"" of type '" "double""'");
  } 
  arg10 = static_cast< double >(val10);
  result = (double)calibration_objective(arg1,(std::vector< std::vector< double,std::allocator< double > >,std::allocator< std::vector< double,std::allocator< double > > > > const &)*arg2,(std::string const &)*arg3,(std::vector< double,std::allocator< double > > const &)*arg4,arg5,arg6,arg7,arg8,arg9,arg10);
  resultobj = SWIG_From_double(static_cast< double >(result));
  if (SWIG_IsNewObj(res3)) delete arg3;
  if (SWIG_IsNewObj(res4)) delete arg4;
  return resultobj;
fail:
  if (SWIG_IsNewObj(res3)) delete arg3;
  if (SWIG_IsNewObj(res4)) delete arg4;
  return NULL;
}


static PyMethodDef SwigMethods[] = {
	 { "SWIG_PyInstanceMethod_New", SWIG_PyInstanceMethod_New, METH_O, NULL},
	 { "delete_SwigPyIterator", _wrap_delete_SwigPyIterator, METH_O, NULL},
	 { "SwigPyIterator_value", _wrap_SwigPyIterator_value, METH_O, NULL},
	 { "SwigPyIterator_incr", _wrap_SwigPyIterator_incr, METH_VARARGS, NULL},
	 { "SwigPyIterator_decr", _wrap_SwigPyIterator_decr, METH_VARARGS, NULL},
	 { "SwigPyIterator_distance", _wrap_SwigPyIterator_distance, METH_VARARGS, NULL},
	 { "SwigPyIterator_equal", _wrap_SwigPyIterator_equal, METH_VARARGS, NULL},
	 { "SwigPyIterator_copy", _wrap_SwigPyIterator_copy, METH_O, NULL},
	 { "SwigPyIterator_next", _wrap_SwigPyIterator_next, METH_O, NULL},
	 { "SwigPyIterator___next__", _wrap_SwigPyIterator___next__, METH_O, NULL},
	 { "SwigPyIterator_previous", _wrap_SwigPyIterator_previous, METH_O, NULL},
	 { "SwigPyIterator_advance", _wrap_SwigPyIterator_advance, METH_VARARGS, NULL},
	 { "SwigPyIterator___eq__", _wrap_SwigPyIterator___eq__, METH_VARARGS, NULL},
	 { "SwigPyIterator___ne__", _wrap_SwigPyIterator___ne__, METH_VARARGS, NULL},
	 { "SwigPyIterator___iadd__", _wrap_SwigPyIterator___iadd__, METH_VARARGS, NULL},
	 { "SwigPyIterator___isub__", _wrap_SwigPyIterator___isub__, METH_VARARGS, NULL},
	 { "SwigPyIterator___add__", _wrap_SwigPyIterator___add__, METH_VARARGS, NULL},
	 { "SwigPyIterator___sub__", _wrap_SwigPyIterator___sub__, METH_VARARGS, NULL},
	 { "SwigPyIterator_swigregister", SwigPyIterator_swigregister, METH_O, NULL},
	 { "VectorDouble_iterator", _wrap_VectorDouble_iterator, METH_O, NULL},
	 { "VectorDouble___nonzero__", _wrap_VectorDouble___nonzero__, METH_O, NULL},
	 { "VectorDouble___bool__", _wrap_VectorDouble___bool__, METH_O, NULL},
	 { "VectorDouble___len__", _wrap_VectorDouble___len__, METH_O, NULL},
	 { "VectorDouble___getslice__", _wrap_VectorDouble___getslice__, METH_VARARGS, NULL},
	 { "VectorDouble___setslice__", _wrap_VectorDouble___setslice__, METH_VARARGS, NULL},
	 { "VectorDouble___delslice__", _wrap_VectorDouble___delslice__, METH_VARARGS, NULL},
	 { "VectorDouble___delitem__", _wrap_VectorDouble___delitem__, METH_VARARGS, NULL},
	 { "VectorDouble___getitem__", _wrap_VectorDouble___getitem__, METH_VARARGS, NULL},
	 { "VectorDouble___setitem__", _wrap_VectorDouble___setitem__, METH_VARARGS, NULL},
	 { "VectorDouble_pop", _wrap_VectorDouble_pop, METH_O, NULL},
	 { "VectorDouble_append", _wrap_VectorDouble_append, METH_VARARGS, NULL},
	 { "VectorDouble_empty", _wrap_VectorDouble_empty, METH_O, NULL},
	 { "VectorDouble_size", _wrap_VectorDouble_size, METH_O, NULL},
	 { "VectorDouble_swap", _wrap_VectorDouble_swap, METH_VARARGS, NULL},
	 { "VectorDouble_begin", _wrap_VectorDouble_begin, METH_O, NULL},
	 { "VectorDouble_end", _wrap_VectorDouble_end, METH_O, NULL},
	 { "VectorDouble_rbegin", _wrap_VectorDouble_rbegin, METH_O, NULL},
	 { "VectorDouble_rend", _wrap_VectorDouble_rend, METH_O, NULL},
	 { "VectorDouble_clear", _wrap_VectorDouble_clear, METH_O, NULL},
	 { "VectorDouble_get_allocator", _wrap_VectorDouble_get_allocator, METH_O, NULL},
	 { "VectorDouble_pop_back", _wrap_VectorDouble_pop_back, METH_O, NULL},
	 { "VectorDouble_erase", _wrap_VectorDouble_erase, METH_VARARGS, NULL},
	 { "new_VectorDouble", _wrap_new_VectorDouble, METH_VARARGS, NULL},
	 { "VectorDouble_push_back", _wrap_VectorDouble_push_back, METH_VARARGS, NULL},
	 { "VectorDouble_front", _wrap_VectorDouble_front, METH_O, NULL},
	 { "VectorDouble_back", _wrap_VectorDouble_back, METH_O, NULL},
	 { "VectorDouble_assign", _wrap_VectorDouble_assign, METH_VARARGS, NULL},
	 { "VectorDouble_resize", _wrap_VectorDouble_resize, METH_VARARGS, NULL},
	 { "VectorDouble_insert", _wrap_VectorDouble_insert, METH_VARARGS, NULL},
	 { "VectorDouble_reserve", _wrap_VectorDouble_reserve, METH_VARARGS, NULL},
	 { "VectorDouble_capacity", _wrap_VectorDouble_capacity, METH_O, NULL},
	 { "delete_VectorDouble", _wrap_delete_VectorDouble, METH_O, NULL},
	 { "VectorDouble_swigregister", VectorDouble_swigregister, METH_O, NULL},
	 { "VectorDouble_swiginit", VectorDouble_swiginit, METH_VARARGS, NULL},
	 { "VectorString_iterator", _wrap_VectorString_iterator, METH_O, NULL},
	 { "VectorString___nonzero__", _wrap_VectorString___nonzero__, METH_O, NULL},
	 { "VectorString___bool__", _wrap_VectorString___bool__, METH_O, NULL},
	 { "VectorString___len__", _wrap_VectorString___len__, METH_O, NULL},
	 { "VectorString___getslice__", _wrap_VectorString___getslice__, METH_VARARGS, NULL},
	 { "VectorString___setslice__", _wrap_VectorString___setslice__, METH_VARARGS, NULL},
	 { "VectorString___delslice__", _wrap_VectorString___delslice__, METH_VARARGS, NULL},
	 { "VectorString___delitem__", _wrap_VectorString___delitem__, METH_VARARGS, NULL},
	 { "VectorString___getitem__", _wrap_VectorString___getitem__, METH_VARARGS, NULL},
	 { "VectorString___setitem__", _wrap_VectorString___setitem__, METH_VARARGS, NULL},
	 { "VectorString_pop", _wrap_VectorString_pop, METH_O, NULL},
	 { "VectorString_append", _wrap_VectorString_append, METH_VARARGS, NULL},
	 { "VectorString_empty", _wrap_VectorString_empty, METH_O, NULL},
	 { "VectorString_size", _wrap_VectorString_size, METH_O, NULL},
	 { "VectorString_swap", _wrap_VectorString_swap, METH_VARARGS, NULL},
	 { "VectorString_begin", _wrap_VectorString_begin, METH_O, NULL},
	 { "VectorString_end", _wrap_VectorString_end, METH_O, NULL},
	 { "VectorString_rbegin", _wrap_VectorString_rbegin, METH_O, NULL},
	 { "VectorString_rend", _wrap_VectorString_rend, METH_O, NULL},
	 { "VectorString_clear", _wrap_VectorString_clear, METH_O, NULL},
	 { "VectorString_get_allocator", _wrap_VectorString_get_allocator, METH_O, NULL},
	 { "VectorString_pop_back", _wrap_VectorString_pop_back, METH_O, NULL},
	 { "VectorString_erase", _wrap_VectorString_erase, METH_VARARGS, NULL},
	 { "new_VectorString", _wrap_new_VectorString, METH_VARARGS, NULL},
	 { "VectorString_push_back", _wrap_VectorString_push_back, METH_VARARGS, NULL},
	 { "VectorString_front", _wrap_VectorString_front, METH_O, NULL},
	 { "VectorString_back", _wrap_VectorString_back, METH_O, NULL},
	 { "VectorString_assign", _wrap_VectorString_assign, METH_VARARGS, NULL},
	 { "VectorString_resize", _wrap_VectorString_resize, METH_VARARGS, NULL},
	 { "VectorString_insert", _wrap_VectorString_insert, METH_VARARGS, NULL},
	 { "VectorString_reserve", _wrap_VectorString_reserve, METH_VARARGS, NULL},
	 { "VectorString_capacity", _wrap_VectorString_capacity, METH_O, NULL},
	 { "delete_VectorString", _wrap_delete_VectorString, METH_O, NULL},
	 { "VectorString_swigregister", VectorString_swigregister, METH_O, NULL},
	 { "VectorString_swiginit", VectorString_swiginit, METH_VARARGS, NULL},
	 { "initial_condition_call", _wrap_initial_condition_call, METH_VARARGS, NULL},
	 { "initial_condition_put", _wrap_initial_condition_put, METH_VARARGS, NULL},
	 { "beta_coefficients", _wrap_beta_coefficients, METH_VARARGS, NULL},
	 { "construct_tridiagonal_coeffs", _wrap_construct_tridiagonal_coeffs, METH_VARARGS, NULL},
	 { "solve_tridiagonal", _wrap_solve_tridiagonal, METH_VARARGS, NULL},
	 { "fractional_black_scholes_alpha", _wrap_fractional_black_scholes_alpha, METH_VARARGS, NULL},
	 { "calibration_objective", _wrap_calibration_objective, METH_VARARGS, NULL},
	 { NULL, NULL, 0, NULL }
};

static PyMethodDef SwigMethods_proxydocs[] = {
	 { NULL, NULL, 0, NULL }
};


/* -------- TYPE CONVERSION AND EQUIVALENCE RULES (BEGIN) -------- */

static swig_type_info _swigt__p_allocator_type = {"_p_allocator_type", "allocator_type *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_char = {"_p_char", "char *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_difference_type = {"_p_difference_type", "difference_type *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_p_PyObject = {"_p_p_PyObject", "PyObject **", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_size_type = {"_p_size_type", "size_type *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_std__allocatorT_double_t = {"_p_std__allocatorT_double_t", "std::vector< double >::allocator_type *|std::allocator< double > *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_std__allocatorT_std__string_t = {"_p_std__allocatorT_std__string_t", "std::vector< std::string >::allocator_type *|std::allocator< std::string > *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_std__invalid_argument = {"_p_std__invalid_argument", "std::invalid_argument *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_std__vectorT_double_std__allocatorT_double_t_t = {"_p_std__vectorT_double_std__allocatorT_double_t_t", "std::vector< double,std::allocator< double > > *|std::vector< double > *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_std__vectorT_std__string_std__allocatorT_std__string_t_t = {"_p_std__vectorT_std__string_std__allocatorT_std__string_t_t", "std::vector< std::string,std::allocator< std::string > > *|std::vector< std::string > *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_std__vectorT_std__vectorT_double_std__allocatorT_double_t_t_std__allocatorT_std__vectorT_double_std__allocatorT_double_t_t_t_t = {"_p_std__vectorT_std__vectorT_double_std__allocatorT_double_t_t_std__allocatorT_std__vectorT_double_std__allocatorT_double_t_t_t_t", "std::vector< std::vector< double,std::allocator< double > >,std::allocator< std::vector< double,std::allocator< double > > > > *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_swig__SwigPyIterator = {"_p_swig__SwigPyIterator", "swig::SwigPyIterator *", 0, 0, (void*)0, 0};
static swig_type_info _swigt__p_value_type = {"_p_value_type", "value_type *", 0, 0, (void*)0, 0};

static swig_type_info *swig_type_initial[] = {
  &_swigt__p_allocator_type,
  &_swigt__p_char,
  &_swigt__p_difference_type,
  &_swigt__p_p_PyObject,
  &_swigt__p_size_type,
  &_swigt__p_std__allocatorT_double_t,
  &_swigt__p_std__allocatorT_std__string_t,
  &_swigt__p_std__invalid_argument,
  &_swigt__p_std__vectorT_double_std__allocatorT_double_t_t,
  &_swigt__p_std__vectorT_std__string_std__allocatorT_std__string_t_t,
  &_swigt__p_std__vectorT_std__vectorT_double_std__allocatorT_double_t_t_std__allocatorT_std__vectorT_double_std__allocatorT_double_t_t_t_t,
  &_swigt__p_swig__SwigPyIterator,
  &_swigt__p_value_type,
};

static swig_cast_info _swigc__p_allocator_type[] = {  {&_swigt__p_allocator_type, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_char[] = {  {&_swigt__p_char, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_difference_type[] = {  {&_swigt__p_difference_type, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_p_PyObject[] = {  {&_swigt__p_p_PyObject, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_size_type[] = {  {&_swigt__p_size_type, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_std__allocatorT_double_t[] = {  {&_swigt__p_std__allocatorT_double_t, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_std__allocatorT_std__string_t[] = {  {&_swigt__p_std__allocatorT_std__string_t, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_std__invalid_argument[] = {  {&_swigt__p_std__invalid_argument, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_std__vectorT_double_std__allocatorT_double_t_t[] = {  {&_swigt__p_std__vectorT_double_std__allocatorT_double_t_t, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_std__vectorT_std__string_std__allocatorT_std__string_t_t[] = {  {&_swigt__p_std__vectorT_std__string_std__allocatorT_std__string_t_t, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_std__vectorT_std__vectorT_double_std__allocatorT_double_t_t_std__allocatorT_std__vectorT_double_std__allocatorT_double_t_t_t_t[] = {  {&_swigt__p_std__vectorT_std__vectorT_double_std__allocatorT_double_t_t_std__allocatorT_std__vectorT_double_std__allocatorT_double_t_t_t_t, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_swig__SwigPyIterator[] = {  {&_swigt__p_swig__SwigPyIterator, 0, 0, 0},{0, 0, 0, 0}};
static swig_cast_info _swigc__p_value_type[] = {  {&_swigt__p_value_type, 0, 0, 0},{0, 0, 0, 0}};

static swig_cast_info *swig_cast_initial[] = {
  _swigc__p_allocator_type,
  _swigc__p_char,
  _swigc__p_difference_type,
  _swigc__p_p_PyObject,
  _swigc__p_size_type,
  _swigc__p_std__allocatorT_double_t,
  _swigc__p_std__allocatorT_std__string_t,
  _swigc__p_std__invalid_argument,
  _swigc__p_std__vectorT_double_std__allocatorT_double_t_t,
  _swigc__p_std__vectorT_std__string_std__allocatorT_std__string_t_t,
  _swigc__p_std__vectorT_std__vectorT_double_std__allocatorT_double_t_t_std__allocatorT_std__vectorT_double_std__allocatorT_double_t_t_t_t,
  _swigc__p_swig__SwigPyIterator,
  _swigc__p_value_type,
};


/* -------- TYPE CONVERSION AND EQUIVALENCE RULES (END) -------- */

static swig_const_info swig_const_table[] = {
{0, 0, 0, 0.0, 0, 0}};

#ifdef __cplusplus
}
#endif
/* -----------------------------------------------------------------------------
 * Type initialization:
 * This problem is tough by the requirement that no dynamic
 * memory is used. Also, since swig_type_info structures store pointers to
 * swig_cast_info structures and swig_cast_info structures store pointers back
 * to swig_type_info structures, we need some lookup code at initialization.
 * The idea is that swig generates all the structures that are needed.
 * The runtime then collects these partially filled structures.
 * The SWIG_InitializeModule function takes these initial arrays out of
 * swig_module, and does all the lookup, filling in the swig_module.types
 * array with the correct data and linking the correct swig_cast_info
 * structures together.
 *
 * The generated swig_type_info structures are assigned statically to an initial
 * array. We just loop through that array, and handle each type individually.
 * First we lookup if this type has been already loaded, and if so, use the
 * loaded structure instead of the generated one. Then we have to fill in the
 * cast linked list. The cast data is initially stored in something like a
 * two-dimensional array. Each row corresponds to a type (there are the same
 * number of rows as there are in the swig_type_initial array). Each entry in
 * a column is one of the swig_cast_info structures for that type.
 * The cast_initial array is actually an array of arrays, because each row has
 * a variable number of columns. So to actually build the cast linked list,
 * we find the array of casts associated with the type, and loop through it
 * adding the casts to the list. The one last trick we need to do is making
 * sure the type pointer in the swig_cast_info struct is correct.
 *
 * First off, we lookup the cast->type name to see if it is already loaded.
 * There are three cases to handle:
 *  1) If the cast->type has already been loaded AND the type we are adding
 *     casting info to has not been loaded (it is in this module), THEN we
 *     replace the cast->type pointer with the type pointer that has already
 *     been loaded.
 *  2) If BOTH types (the one we are adding casting info to, and the
 *     cast->type) are loaded, THEN the cast info has already been loaded by
 *     the previous module so we just ignore it.
 *  3) Finally, if cast->type has not already been loaded, then we add that
 *     swig_cast_info to the linked list (because the cast->type) pointer will
 *     be correct.
 * ----------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#if 0
} /* c-mode */
#endif
#endif

#if 0
#define SWIGRUNTIME_DEBUG
#endif


SWIGRUNTIME void
SWIG_InitializeModule(void *clientdata) {
  size_t i;
  swig_module_info *module_head, *iter;
  int init;
  
  /* check to see if the circular list has been setup, if not, set it up */
  if (swig_module.next==0) {
    /* Initialize the swig_module */
    swig_module.type_initial = swig_type_initial;
    swig_module.cast_initial = swig_cast_initial;
    swig_module.next = &swig_module;
    init = 1;
  } else {
    init = 0;
  }
  
  /* Try and load any already created modules */
  module_head = SWIG_GetModule(clientdata);
  if (!module_head) {
    /* This is the first module loaded for this interpreter */
    /* so set the swig module into the interpreter */
    SWIG_SetModule(clientdata, &swig_module);
  } else {
    /* the interpreter has loaded a SWIG module, but has it loaded this one? */
    iter=module_head;
    do {
      if (iter==&swig_module) {
        /* Our module is already in the list, so there's nothing more to do. */
        return;
      }
      iter=iter->next;
    } while (iter!= module_head);
    
    /* otherwise we must add our module into the list */
    swig_module.next = module_head->next;
    module_head->next = &swig_module;
  }
  
  /* When multiple interpreters are used, a module could have already been initialized in
       a different interpreter, but not yet have a pointer in this interpreter.
       In this case, we do not want to continue adding types... everything should be
       set up already */
  if (init == 0) return;
  
  /* Now work on filling in swig_module.types */
#ifdef SWIGRUNTIME_DEBUG
  printf("SWIG_InitializeModule: size %lu\n", (unsigned long)swig_module.size);
#endif
  for (i = 0; i < swig_module.size; ++i) {
    swig_type_info *type = 0;
    swig_type_info *ret;
    swig_cast_info *cast;
    
#ifdef SWIGRUNTIME_DEBUG
    printf("SWIG_InitializeModule: type %lu %s\n", (unsigned long)i, swig_module.type_initial[i]->name);
#endif
    
    /* if there is another module already loaded */
    if (swig_module.next != &swig_module) {
      type = SWIG_MangledTypeQueryModule(swig_module.next, &swig_module, swig_module.type_initial[i]->name);
    }
    if (type) {
      /* Overwrite clientdata field */
#ifdef SWIGRUNTIME_DEBUG
      printf("SWIG_InitializeModule: found type %s\n", type->name);
#endif
      if (swig_module.type_initial[i]->clientdata) {
        type->clientdata = swig_module.type_initial[i]->clientdata;
#ifdef SWIGRUNTIME_DEBUG
        printf("SWIG_InitializeModule: found and overwrite type %s \n", type->name);
#endif
      }
    } else {
      type = swig_module.type_initial[i];
    }
    
    /* Insert casting types */
    cast = swig_module.cast_initial[i];
    while (cast->type) {
      /* Don't need to add information already in the list */
      ret = 0;
#ifdef SWIGRUNTIME_DEBUG
      printf("SWIG_InitializeModule: look cast %s\n", cast->type->name);
#endif
      if (swig_module.next != &swig_module) {
        ret = SWIG_MangledTypeQueryModule(swig_module.next, &swig_module, cast->type->name);
#ifdef SWIGRUNTIME_DEBUG
        if (ret) printf("SWIG_InitializeModule: found cast %s\n", ret->name);
#endif
      }
      if (ret) {
        if (type == swig_module.type_initial[i]) {
#ifdef SWIGRUNTIME_DEBUG
          printf("SWIG_InitializeModule: skip old type %s\n", ret->name);
#endif
          cast->type = ret;
          ret = 0;
        } else {
          /* Check for casting already in the list */
          swig_cast_info *ocast = SWIG_TypeCheck(ret->name, type);
#ifdef SWIGRUNTIME_DEBUG
          if (ocast) printf("SWIG_InitializeModule: skip old cast %s\n", ret->name);
#endif
          if (!ocast) ret = 0;
        }
      }
      
      if (!ret) {
#ifdef SWIGRUNTIME_DEBUG
        printf("SWIG_InitializeModule: adding cast %s\n", cast->type->name);
#endif
        if (type->cast) {
          type->cast->prev = cast;
          cast->next = type->cast;
        }
        type->cast = cast;
      }
      cast++;
    }
    /* Set entry in modules->types array equal to the type */
    swig_module.types[i] = type;
  }
  swig_module.types[i] = 0;
  
#ifdef SWIGRUNTIME_DEBUG
  printf("**** SWIG_InitializeModule: Cast List ******\n");
  for (i = 0; i < swig_module.size; ++i) {
    int j = 0;
    swig_cast_info *cast = swig_module.cast_initial[i];
    printf("SWIG_InitializeModule: type %lu %s\n", (unsigned long)i, swig_module.type_initial[i]->name);
    while (cast->type) {
      printf("SWIG_InitializeModule: cast type %s\n", cast->type->name);
      cast++;
      ++j;
    }
    printf("---- Total casts: %d\n",j);
  }
  printf("**** SWIG_InitializeModule: Cast List ******\n");
#endif
}

/* This function will propagate the clientdata field of type to
* any new swig_type_info structures that have been added into the list
* of equivalent types.  It is like calling
* SWIG_TypeClientData(type, clientdata) a second time.
*/
SWIGRUNTIME void
SWIG_PropagateClientData(void) {
  size_t i;
  swig_cast_info *equiv;
  static int init_run = 0;
  
  if (init_run) return;
  init_run = 1;
  
  for (i = 0; i < swig_module.size; i++) {
    if (swig_module.types[i]->clientdata) {
      equiv = swig_module.types[i]->cast;
      while (equiv) {
        if (!equiv->converter) {
          if (equiv->type && !equiv->type->clientdata)
          SWIG_TypeClientData(equiv->type, swig_module.types[i]->clientdata);
        }
        equiv = equiv->next;
      }
    }
  }
}

#ifdef __cplusplus
#if 0
{
  /* c-mode */
#endif
}
#endif



#ifdef __cplusplus
extern "C" {
#endif
  
  /* Python-specific SWIG API */
#define SWIG_newvarlink()                             SWIG_Python_newvarlink()
#define SWIG_addvarlink(p, name, get_attr, set_attr)  SWIG_Python_addvarlink(p, name, get_attr, set_attr)
#define SWIG_InstallConstants(d, constants)           SWIG_Python_InstallConstants(d, constants)
  
  /* -----------------------------------------------------------------------------
   * global variable support code.
   * ----------------------------------------------------------------------------- */
  
  typedef struct swig_globalvar {
    char       *name;                  /* Name of global variable */
    PyObject *(*get_attr)(void);       /* Return the current value */
    int       (*set_attr)(PyObject *); /* Set the value */
    struct swig_globalvar *next;
  } swig_globalvar;
  
  typedef struct swig_varlinkobject {
    PyObject_HEAD
    swig_globalvar *vars;
  } swig_varlinkobject;
  
  SWIGINTERN PyObject *
  swig_varlink_repr(swig_varlinkobject *SWIGUNUSEDPARM(v)) {
#if PY_VERSION_HEX >= 0x03000000
    return PyUnicode_InternFromString("<Swig global variables>");
#else
    return PyString_FromString("<Swig global variables>");
#endif
  }
  
  SWIGINTERN PyObject *
  swig_varlink_str(swig_varlinkobject *v) {
#if PY_VERSION_HEX >= 0x03000000
    PyObject *str = PyUnicode_InternFromString("(");
    PyObject *tail;
    PyObject *joined;
    swig_globalvar *var;
    for (var = v->vars; var; var=var->next) {
      tail = PyUnicode_FromString(var->name);
      joined = PyUnicode_Concat(str, tail);
      Py_DecRef(str);
      Py_DecRef(tail);
      str = joined;
      if (var->next) {
        tail = PyUnicode_InternFromString(", ");
        joined = PyUnicode_Concat(str, tail);
        Py_DecRef(str);
        Py_DecRef(tail);
        str = joined;
      }
    }
    tail = PyUnicode_InternFromString(")");
    joined = PyUnicode_Concat(str, tail);
    Py_DecRef(str);
    Py_DecRef(tail);
    str = joined;
#else
    PyObject *str = PyString_FromString("(");
    swig_globalvar *var;
    for (var = v->vars; var; var=var->next) {
      PyString_ConcatAndDel(&str,PyString_FromString(var->name));
      if (var->next) PyString_ConcatAndDel(&str,PyString_FromString(", "));
    }
    PyString_ConcatAndDel(&str,PyString_FromString(")"));
#endif
    return str;
  }
  
  SWIGINTERN void
  swig_varlink_dealloc(swig_varlinkobject *v) {
    swig_globalvar *var = v->vars;
    while (var) {
      swig_globalvar *n = var->next;
      free(var->name);
      free(var);
      var = n;
    }
  }
  
  SWIGINTERN PyObject *
  swig_varlink_getattr(swig_varlinkobject *v, char *n) {
    PyObject *res = NULL;
    swig_globalvar *var = v->vars;
    while (var) {
      if (strcmp(var->name,n) == 0) {
        res = (*var->get_attr)();
        break;
      }
      var = var->next;
    }
    if (res == NULL && !PyErr_Occurred()) {
      PyErr_Format(PyExc_AttributeError, "Unknown C global variable '%s'", n);
    }
    return res;
  }
  
  SWIGINTERN int
  swig_varlink_setattr(swig_varlinkobject *v, char *n, PyObject *p) {
    int res = 1;
    swig_globalvar *var = v->vars;
    while (var) {
      if (strcmp(var->name,n) == 0) {
        res = (*var->set_attr)(p);
        break;
      }
      var = var->next;
    }
    if (res == 1 && !PyErr_Occurred()) {
      PyErr_Format(PyExc_AttributeError, "Unknown C global variable '%s'", n);
    }
    return res;
  }
  
  SWIGINTERN PyTypeObject*
  swig_varlink_type(void) {
    static char varlink__doc__[] = "Swig var link object";
    static PyTypeObject varlink_type;
    static int type_init = 0;
    if (!type_init) {
      const PyTypeObject tmp = {
#if PY_VERSION_HEX >= 0x03000000
        PyVarObject_HEAD_INIT(NULL, 0)
#else
        PyObject_HEAD_INIT(NULL)
        0,                                  /* ob_size */
#endif
        "swigvarlink",                      /* tp_name */
        sizeof(swig_varlinkobject),         /* tp_basicsize */
        0,                                  /* tp_itemsize */
        (destructor) swig_varlink_dealloc,  /* tp_dealloc */
        0,                                  /* tp_print */
        (getattrfunc) swig_varlink_getattr, /* tp_getattr */
        (setattrfunc) swig_varlink_setattr, /* tp_setattr */
        0,                                  /* tp_compare */
        (reprfunc) swig_varlink_repr,       /* tp_repr */
        0,                                  /* tp_as_number */
        0,                                  /* tp_as_sequence */
        0,                                  /* tp_as_mapping */
        0,                                  /* tp_hash */
        0,                                  /* tp_call */
        (reprfunc) swig_varlink_str,        /* tp_str */
        0,                                  /* tp_getattro */
        0,                                  /* tp_setattro */
        0,                                  /* tp_as_buffer */
        0,                                  /* tp_flags */
        varlink__doc__,                     /* tp_doc */
        0,                                  /* tp_traverse */
        0,                                  /* tp_clear */
        0,                                  /* tp_richcompare */
        0,                                  /* tp_weaklistoffset */
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, /* tp_iter -> tp_weaklist */
        0,                                  /* tp_del */
        0,                                  /* tp_version_tag */
#if PY_VERSION_HEX >= 0x03040000
        0,                                  /* tp_finalize */
#endif
#if PY_VERSION_HEX >= 0x03080000
        0,                                  /* tp_vectorcall */
#endif
#if (PY_VERSION_HEX >= 0x03080000) && (PY_VERSION_HEX < 0x03090000)
        0,                                  /* tp_print */
#endif
#ifdef COUNT_ALLOCS
        0,                                  /* tp_allocs */
        0,                                  /* tp_frees */
        0,                                  /* tp_maxalloc */
        0,                                  /* tp_prev */
        0                                   /* tp_next */
#endif
      };
      varlink_type = tmp;
      type_init = 1;
      if (PyType_Ready(&varlink_type) < 0)
      return NULL;
    }
    return &varlink_type;
  }
  
  /* Create a variable linking object for use later */
  SWIGINTERN PyObject *
  SWIG_Python_newvarlink(void) {
    swig_varlinkobject *result = PyObject_NEW(swig_varlinkobject, swig_varlink_type());
    if (result) {
      result->vars = 0;
    }
    return ((PyObject*) result);
  }
  
  SWIGINTERN void 
  SWIG_Python_addvarlink(PyObject *p, const char *name, PyObject *(*get_attr)(void), int (*set_attr)(PyObject *p)) {
    swig_varlinkobject *v = (swig_varlinkobject *) p;
    swig_globalvar *gv = (swig_globalvar *) malloc(sizeof(swig_globalvar));
    if (gv) {
      size_t size = strlen(name)+1;
      gv->name = (char *)malloc(size);
      if (gv->name) {
        memcpy(gv->name, name, size);
        gv->get_attr = get_attr;
        gv->set_attr = set_attr;
        gv->next = v->vars;
      }
    }
    v->vars = gv;
  }
  
  SWIGINTERN PyObject *
  SWIG_globals(void) {
    static PyObject *globals = 0;
    if (!globals) {
      globals = SWIG_newvarlink();
    }
    return globals;
  }
  
  /* -----------------------------------------------------------------------------
   * constants/methods manipulation
   * ----------------------------------------------------------------------------- */
  
  /* Install Constants */
  SWIGINTERN void
  SWIG_Python_InstallConstants(PyObject *d, swig_const_info constants[]) {
    PyObject *obj = 0;
    size_t i;
    for (i = 0; constants[i].type; ++i) {
      switch(constants[i].type) {
      case SWIG_PY_POINTER:
        obj = SWIG_InternalNewPointerObj(constants[i].pvalue, *(constants[i]).ptype,0);
        break;
      case SWIG_PY_BINARY:
        obj = SWIG_NewPackedObj(constants[i].pvalue, constants[i].lvalue, *(constants[i].ptype));
        break;
      default:
        obj = 0;
        break;
      }
      if (obj) {
        PyDict_SetItemString(d, constants[i].name, obj);
        Py_DECREF(obj);
      }
    }
  }
  
  /* -----------------------------------------------------------------------------*/
  /* Fix SwigMethods to carry the callback ptrs when needed */
  /* -----------------------------------------------------------------------------*/
  
  SWIGINTERN void
  SWIG_Python_FixMethods(PyMethodDef *methods,
    swig_const_info *const_table,
    swig_type_info **types,
    swig_type_info **types_initial) {
    size_t i;
    for (i = 0; methods[i].ml_name; ++i) {
      const char *c = methods[i].ml_doc;
      if (!c) continue;
      c = strstr(c, "swig_ptr: ");
      if (c) {
        int j;
        swig_const_info *ci = 0;
        const char *name = c + 10;
        for (j = 0; const_table[j].type; ++j) {
          if (strncmp(const_table[j].name, name, 
              strlen(const_table[j].name)) == 0) {
            ci = &(const_table[j]);
            break;
          }
        }
        if (ci) {
          void *ptr = (ci->type == SWIG_PY_POINTER) ? ci->pvalue : 0;
          if (ptr) {
            size_t shift = (ci->ptype) - types;
            swig_type_info *ty = types_initial[shift];
            size_t ldoc = (c - methods[i].ml_doc);
            size_t lptr = strlen(ty->name)+2*sizeof(void*)+2;
            char *ndoc = (char*)malloc(ldoc + lptr + 10);
            if (ndoc) {
              char *buff = ndoc;
              memcpy(buff, methods[i].ml_doc, ldoc);
              buff += ldoc;
              memcpy(buff, "swig_ptr: ", 10);
              buff += 10;
              SWIG_PackVoidPtr(buff, ptr, ty->name, lptr);
              methods[i].ml_doc = ndoc;
            }
          }
        }
      }
    }
  } 
  
  /* -----------------------------------------------------------------------------
   * Method creation and docstring support functions
   * ----------------------------------------------------------------------------- */
  
  /* -----------------------------------------------------------------------------
   * Function to find the method definition with the correct docstring for the
   * proxy module as opposed to the low-level API
   * ----------------------------------------------------------------------------- */
  
  SWIGINTERN PyMethodDef *SWIG_PythonGetProxyDoc(const char *name) {
    /* Find the function in the modified method table */
    size_t offset = 0;
    int found = 0;
    while (SwigMethods_proxydocs[offset].ml_meth != NULL) {
      if (strcmp(SwigMethods_proxydocs[offset].ml_name, name) == 0) {
        found = 1;
        break;
      }
      offset++;
    }
    /* Use the copy with the modified docstring if available */
    return found ? &SwigMethods_proxydocs[offset] : NULL;
  }
  
  /* -----------------------------------------------------------------------------
   * Wrapper of PyInstanceMethod_New() used in Python 3
   * It is exported to the generated module, used for -fastproxy
   * ----------------------------------------------------------------------------- */
  
  SWIGINTERN PyObject *SWIG_PyInstanceMethod_New(PyObject *SWIGUNUSEDPARM(self), PyObject *func) {
    if (PyCFunction_Check(func)) {
      PyCFunctionObject *funcobj = (PyCFunctionObject *)func;
      PyMethodDef *ml = SWIG_PythonGetProxyDoc(funcobj->m_ml->ml_name);
      if (ml)
      func = PyCFunction_NewEx(ml, funcobj->m_self, funcobj->m_module);
    }
#if PY_VERSION_HEX >= 0x03000000
    return PyInstanceMethod_New(func);
#else
    return PyMethod_New(func, NULL, NULL);
#endif
  }
  
  /* -----------------------------------------------------------------------------
   * Wrapper of PyStaticMethod_New()
   * It is exported to the generated module, used for -fastproxy
   * ----------------------------------------------------------------------------- */
  
  SWIGINTERN PyObject *SWIG_PyStaticMethod_New(PyObject *SWIGUNUSEDPARM(self), PyObject *func) {
    if (PyCFunction_Check(func)) {
      PyCFunctionObject *funcobj = (PyCFunctionObject *)func;
      PyMethodDef *ml = SWIG_PythonGetProxyDoc(funcobj->m_ml->ml_name);
      if (ml)
      func = PyCFunction_NewEx(ml, funcobj->m_self, funcobj->m_module);
    }
    return PyStaticMethod_New(func);
  }
  
#ifdef __cplusplus
}
#endif

/* -----------------------------------------------------------------------------*
 *  Partial Init method
 * -----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
#endif

SWIGEXPORT 
#if PY_VERSION_HEX >= 0x03000000
PyObject*
#else
void
#endif
SWIG_init(void) {
  PyObject *m, *d, *md, *globals;
  
#if PY_VERSION_HEX >= 0x03000000
  static struct PyModuleDef SWIG_module = {
    PyModuleDef_HEAD_INIT,
    SWIG_name,
    NULL,
    -1,
    SwigMethods,
    NULL,
    NULL,
    NULL,
    NULL
  };
#endif
  
#if defined(SWIGPYTHON_BUILTIN)
  static SwigPyClientData SwigPyObject_clientdata = {
    0, 0, 0, 0, 0, 0, 0
  };
  static PyGetSetDef this_getset_def = {
    (char *)"this", &SwigPyBuiltin_ThisClosure, NULL, NULL, NULL
  };
  static SwigPyGetSet thisown_getset_closure = {
    SwigPyObject_own,
    SwigPyObject_own
  };
  static PyGetSetDef thisown_getset_def = {
    (char *)"thisown", SwigPyBuiltin_GetterClosure, SwigPyBuiltin_SetterClosure, NULL, &thisown_getset_closure
  };
  PyTypeObject *builtin_pytype;
  int builtin_base_count;
  swig_type_info *builtin_basetype;
  PyObject *tuple;
  PyGetSetDescrObject *static_getset;
  PyTypeObject *metatype;
  PyTypeObject *swigpyobject;
  SwigPyClientData *cd;
  PyObject *public_interface, *public_symbol;
  PyObject *this_descr;
  PyObject *thisown_descr;
  PyObject *self = 0;
  int i;
  
  (void)builtin_pytype;
  (void)builtin_base_count;
  (void)builtin_basetype;
  (void)tuple;
  (void)static_getset;
  (void)self;
  
  /* Metaclass is used to implement static member variables */
  metatype = SwigPyObjectType();
  assert(metatype);
#endif
  
  (void)globals;
  
  /* Create singletons now to avoid potential deadlocks with multi-threaded usage after module initialization */
  SWIG_This();
  SWIG_Python_TypeCache();
  SwigPyPacked_type();
#ifndef SWIGPYTHON_BUILTIN
  SwigPyObject_type();
#endif
  
  /* Fix SwigMethods to carry the callback ptrs when needed */
  SWIG_Python_FixMethods(SwigMethods, swig_const_table, swig_types, swig_type_initial);
  
#if PY_VERSION_HEX >= 0x03000000
  m = PyModule_Create(&SWIG_module);
#else
  m = Py_InitModule(SWIG_name, SwigMethods);
#endif
  
  md = d = PyModule_GetDict(m);
  (void)md;
  
  SWIG_InitializeModule(0);
  
#ifdef SWIGPYTHON_BUILTIN
  swigpyobject = SwigPyObject_TypeOnce();
  
  SwigPyObject_stype = SWIG_MangledTypeQuery("_p_SwigPyObject");
  assert(SwigPyObject_stype);
  cd = (SwigPyClientData*) SwigPyObject_stype->clientdata;
  if (!cd) {
    SwigPyObject_stype->clientdata = &SwigPyObject_clientdata;
    SwigPyObject_clientdata.pytype = swigpyobject;
  } else if (swigpyobject->tp_basicsize != cd->pytype->tp_basicsize) {
    PyErr_SetString(PyExc_RuntimeError, "Import error: attempted to load two incompatible swig-generated modules.");
# if PY_VERSION_HEX >= 0x03000000
    return NULL;
# else
    return;
# endif
  }
  
  /* All objects have a 'this' attribute */
  this_descr = PyDescr_NewGetSet(SwigPyObject_type(), &this_getset_def);
  (void)this_descr;
  
  /* All objects have a 'thisown' attribute */
  thisown_descr = PyDescr_NewGetSet(SwigPyObject_type(), &thisown_getset_def);
  (void)thisown_descr;
  
  public_interface = PyList_New(0);
  public_symbol = 0;
  (void)public_symbol;
  
  PyDict_SetItemString(md, "__all__", public_interface);
  Py_DECREF(public_interface);
  for (i = 0; SwigMethods[i].ml_name != NULL; ++i)
  SwigPyBuiltin_AddPublicSymbol(public_interface, SwigMethods[i].ml_name);
  for (i = 0; swig_const_table[i].name != 0; ++i)
  SwigPyBuiltin_AddPublicSymbol(public_interface, swig_const_table[i].name);
#endif
  
  SWIG_InstallConstants(d,swig_const_table);
  
  
  // thread safe initialization
  swig::container_owner_attribute();
  
#if PY_VERSION_HEX >= 0x03000000
  return m;
#else
  return;
#endif
}

