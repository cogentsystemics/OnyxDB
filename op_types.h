


#ifndef __OPTYPES_H_INCLUDED__
#define __OPTYPES_H_INCLUDED__   13

#ifdef __cplusplus
extern "C" {
#endif

#define deftype(X, Y)    typedef typeof(X) Y

#ifndef _unixcube_openware_op_typestruct_def_
    #define _unixcube_openware_op_typestruct_def_
    #define typestruct  typedef struct
#endif

#define A(x)(y)  A[x][y]

#define _a(A)       [A]
#define m2(A,B)     [A][B]
#define m3(A,B,C)   [A][B][C]
  
#define block_   {
#define _block   }
  
#define loop     for
//#define loop(a){b;...}     for(a;b;c)b,...;


#define ret(X)    return(X)
#define rt(X)    ret(X)
  
#define _u unsigned
#define _s static
#define _v void
#define _sh short
#define _sn signed

#define _c char
#define _i int
#define _l long
#define _f float
#define _d double
#define _ll _l _l

#define _uc  _u _c
#define _ui  _u _i
#define _ul  _u _l
#define _uf  _u _f
#define _ud  _u _d
#define _ull _u _ll


#define _co const
#define _cc const char
#define _ci const int
#define _cf const float
#define _cd const _d
#define _cl _co _l

#define re(r)         return(r)

#define outc_(m)      printf("debug msg:%s\n",m)
#define outi_(n)      printf("debug num:%i\n",n)
#define outs_(s)      printf(s)

#define outc outc_
#define outi outi_
#define outs outs_

#define cnew(a,x)       (typeof(a)*)malloc(sizeof(typeof(a))*x)
#define new_c(x) (_c*)malloc((sizeof(_c)*x)
#define new_i(x) (_i*)malloc((sizeof(_i)*x)
#define new_l(x) (_l*)malloc((sizeof(_l)*x)
#define new_ul(x) (_ul*)malloc((sizeof(_ul)*x)
#define new_ull(x) (_ull*)malloc((sizeof(_ull)*x)

#define atype(s)         (typeof(s)*)malloc(sizeof(typeof(s)))
#define atypex(s,n)      (typeof(s)*)malloc( ( sizeof(typeof(s) ) *n)
#define atypesz(s,n)      (typeof(s)*)malloc( (( sizeof(typeof(s) ) *n)+1)
#define o(x) printf(x)

#ifdef __cplusplus
}
#endif

#endif
