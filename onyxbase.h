

/*   onyxbase.h :: suppliment to onyxcore.h, provides types and function primitives for
 * for parsing, organizing and querying databases with integrated inline scripting language
 * based on string, data and number processing functionality ....
 */




#include <unixcube/op_types.h>
#include <unixcube/op_fio.h>

#include <gmp.h>
#include <unixcube/op_basex3.h>
#include <unixcube/onyxcore.h>


typedef struct onyxword_t onyxword;
typedef struct onyxnum_t onyxnum;
typedef struct onyxbase_t onyxbase;
typedef struct onyxquery_t onyxquery;


static _c* dlm_sz[]={
    "</!>",
    "{{",
    "}}",
    "str[",
    "]str",
    "{/}"
};

static _c* dlm2_sz[]={
    ";",
    ":",
    "[",
    "]",
    "((",
    "))"
};

static _c* dlm3_sz[]={
    ":--",
    "--:",
    "::-",
    "-::",
    ":-",
    "-:"
};

static _c* onyx_keywords[]={
  "fieldnames",
    "num_base",
    "baseX",
    "base10",
    "base2",
    "base16",
    "base13",
    "base18",
    "basew",
    "base32",
    "base44",
    "base64"
};
static _c* excl_sz0 = "n";

static _c* rep_sz0 = "s";

_c* ch_indx_lwr[] = {"abcdefghijklmnopqrstuvwxyz","ABCDEFGHIJKLMNOPQRSTUVWXYZ"};
_c* ch_indx_upr[] = {"ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz"};

_c n_sz[] = {"0123456789abcdefghijklmnopqrstuvwxyz"};

_c* case_lwr = "abcdefghijklmnopqrstuvwxyz";
_c* case_upr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

_ul chrcmp_ci(_c c0, _c c1)
{
    _ul i0;
    _ul len0 = szlen2(case_lwr);
    _ul match0, match1;

    match0 = match1 = 0;
    for(i0 = 0 ; i0 < len0; ++i0){
        if(c0 == case_lwr[i0] || c0 == case_upr[i0] ){
            match0 = 1;
            break;
        }
    }
    if(match0){
        if( c1 == case_lwr[i0] || c1 == case_upr[i0] )
            match1 = 1;

    }else{
        if( c0 == c1)
            match1 = 1;

    }
    return(match1);
}

_ul strcmp_ci(_c* s0, _c* s1 )
{
    _ul i0;
    _c* sp0;
    _ul len_sz, len_abc;

    len_sz = szlen2(s0);

    _ul match0 = chrcmp_ci(s0[0], s1[0]);

    if(match0){
        for( i0=0; i0 < len_sz && match0; ++i0){
            match0 = chrcmp_ci(s0[i0], s1[i0] );

        }
    }

    return(match0);
}

_ul slen(_c* s)
{
    _ul i0=0;
    _c* sp0=0;
    for( sp0 = s; *sp0; ++sp0, ++i0);
    return(i0);
}


_i dlmchk(_c* s0, _c* dlm)
{
    _ul i0, i1;
    _c* sp0, *sp1;
    _i match=0;

    for(sp0 = s0, sp1 = dlm ; *sp0 && *sp1; ++sp0, ++sp1){
        if(*sp0 == *sp1){
            match=1;
        }else{
            match=0;
            break;
        }
    }
    return(match);
}

typedef struct onyxword_t {
    _c* word;
    _ul wlen;

    mpz_t* mp_n;
    onyxnum* num;
    _ul n, nBase, is_n, bitlen;

    _v** vpx;
    _v* rstack, *_vp0, *_vp1, *_vp2;
    _v* (*func)(_v);
    _v* (*func_)(void*, ...);
    _ul vpx_c, vpx_typ, vpx_stat;

    struct onyxword_t* master, *stem;
    struct onyxword_t* head, *prev, *next;
    _ul ref[8];
    _ul index, stem_index, archetype, subtype, state;

    _c* szptr_orig, *szptr_next;
}onyxword;

onyxword* onyxword_new(){ return( (onyxword*)malloc(sizeof(onyxword))  ); }
_v onyxword_setptrs(onyxword* wnod, onyxword* hed, onyxword* pre, onyxword* nex)
{ wnod->head = hed;wnod->prev=pre;wnod->next=nex;  }
_v onyxword_setmaster( onyxword* wnod, onyxword *mas, onyxword* stm)
   { wnod->master=mas; wnod->stem=stm; }

typedef struct onyxnum_t {
    _ul n32;
    _c* n_sz;
    _ul nsz_len;

    mpz_t* n;
    _ul n_base, bitwidth;

    struct onyxnum_t* head, *prev, *next;

    _ul index, archetype, subtype, state;

}onyxnum;

typedef struct onyxbase_t {
    _c* nam, *dat;
    _ul namlen, datlen, n;

    onyxword* wordset;
    _ul words;
    
    struct onyxbase_t* master, *stem;
    struct onyxbase_t* head, *prev, *next;
    struct onyxbase_t* db_x, *db_y, *db_z;

    _ul index, stem_index, archetype, subtype, state;
}onyxbase;

onyxbase* onyxbase_new(){ return( (onyxbase*)malloc(sizeof(onyxbase))  ); }
_v onyxbase_setptrs(onyxbase* nod, onyxbase* hed, onyxbase* pre, onyxbase* nex )
{  nod->head = hed; nod->prev=pre; nod->next=nex;  }
_v onyxbase_setmaster( onyxbase* nod, onyxbase* mas, onyxbase* stm )
{nod->master=mas; nod->stem=stm; }

typedef struct onyxquery_t {
    onyxword* dbQuery;

    onyxbase* dbObj;
    onyxword* dbEntry;
    onyxword* dbWord;

    struct onyxquery_t* master, *stem;
    struct onyxquery_t* head, *prev, *next;

    _ul query_hits, query_index_current;

    _ul index, stem_index, archetype, subtype, state;

}onyxquery;

onyxquery* onyxquery_new(){ return( (onyxquery*)malloc(sizeof(onyxquery)) ); }
_v onyxquery_setptrs( onyxquery* nod, onyxquery* hed, onyxquery* pre, onyxquery* nex )
{ nod->head=hed; nod->prev=pre; nod->next=nex; }

typedef struct onyx_obj_t {

    onyxword* wrdObj;
    onyxbase* dbsObj;
    onyxquery* qryObj;

    onyxword* szRaw, *szExcludes, *szDelims, *szFiles;

    struct onyx_obj_t* master, *stem;
    struct onyx_obj_t* head, *prev, *next;

    _ul index, stem_index, archetype, subtype, state;

}onyx_obj;

onyx_obj* onyx_obj_new(){ return( (onyx_obj*)malloc( sizeof(onyx_obj) )  ); }
_v onyx_obj_setptrs( onyx_obj* obj, onyx_obj* hed, onyx_obj* pre, onyx_obj* nex )
{ obj->head=hed; obj->prev=pre; obj->next=nex; }


typedef struct onyxcode_object_t {
    _c* nam, dat;
    _ul namlen, datlen;

    onyxbase* funcObj, *varObj, *shellObj, *fileObj;
    onyxquery* queryObj;
    _ul nFunc, nVar, nShell, nFile, nQuery;

    struct onyxcode_object_t* head, *prev, *next;

    _ul index, archetype, subtype, obj_state, exec_state;

}onyxcode_object;
onyxcode_object* onyxcode_object_new(){ return( (onyxcode_object*)malloc(sizeof(onyxcode_object))); }

typedef struct onyx_functor_t{
    _c* func_name;

    _v* (*func)(_v);
    _v* (*func_)(void*, ...);

    _ul param_count;
    _ul param_indx;
    onyxword* param_set;

    struct onyx_functor_t* master, *head, *prev, *next, *tail;
}onyx_functor;

typedef struct onyxfunc_t {
    onyxword* sz_blocs; // <func> NAME { stackframe_return } { stackframe_pass} { instruction_set }
    onyxword* ret_stack, *live_stack, *inp_stack;
    onyxword* lineset;

    _ul is_cyclic, cycle_state, cycle_indx;

    struct onyxfunc_t* master, *stem;
    struct onyxfunc_t* prev, *spawn;

    _ul index, archetype, subtype, state;

}onyxfunc;

typedef struct onyxcode_expr_t {
    onyxword* words;
    onyxword* stack;
    onyxword* retvals;

    _ul archetype, subtype, superstate, substate;

}onyxcode_expr;

onyxcode_expr* onyx_expr_new(){ return( (onyxcode_expr*)malloc(sizeof(onyxcode_expr) ) ); }
_v onyxcode_expr_setptrs( onyxcode_expr* nod, onyxcode_expr* wrd, onyxcode_expr* stk, onyxcode_expr* ret )
{ nod->words=wrd;nod->stack=stk;nod->retvals=ret; }

typedef struct onyxcode_line_t {
    onyxcode_expr* expr;
    onyxbase* dbObj;
    struct onyxcode_line_t* head, *prev, *next;

    _ul index, archetype, subtype, superstate, substate;
    _i is_final;
}onyxcode_line;

onyxcode_line* onyxcode_line_new()
{
    return( (onyxcode_line*) malloc( sizeof(onyxcode_line )));
}

onyxcode_line* onyxcode_line_setptrs( onyxcode_line* nod, onyxbase* db, onyxcode_expr* expr, onyxcode_line* hed,
                                      onyxcode_line* pre, onyxcode_line* nex)
{   nod->dbObj=db; nod->expr=expr;nod->head=hed;nod->prev=pre;nod->next=nex; }


onyxfunc* onyxfunc_new(){    return(  (onyxfunc*) malloc( sizeof(onyxfunc)));}

onyxfunc* onyx_function_obj_constuct( _c* code, onyxword* funcFrame_dlm )
{

    _ul i0, i1, i2;

    onyxfunc* fObj = onyxfunc_new();
    onyxword* w0, *w1, *w2, *w3;
    onyxbase* db0, *db1, *db2;

    // parse blocks of func parameters ( name, ret_stack, in_stack, codelines)  into onyxfunc structure ....


}
/*
_ul onyxword_destruct(onyxword* wObj)
{
    _ul len = slen( wObj->word);
    len = ( len != wObj->wlen ) ?
                ( ( len > wObj->wlen ) ? :  (wObj->wlen) ) :
                    ( len ); 

    free(wObj->word );
    free(wObj );

    return(len);
}*/

typedef struct onyxtype_t {
    _c* type_name;
    _ul nlen;
    _c* type_data;
    _ul dlen;

    onyxword* t;
    _ul word_count;
    _v* (*func)(_v*, ...);
    _ul inp_params;

    struct onyxtype_t* head, *prev, *next;
}onyxtype;

onyxtype* onyxtype_new(){ return( (onyxtype*)malloc(sizeof(onyxtype)) ); }

onyxtype* onyxtype_initnew(){
    onyxtype* typ = onyxtype_new();

    return(typ);
}


onyxnum* onyxnum_spawn(_c* szn, _ul base)
{
    _ul i0, i1;
    _c* sp0, *sp1;

    onyxnum* n = (onyxnum*)malloc(sizeof(onyxnum));
    n->n = (mpz_t*)malloc(sizeof(mpz_t));

    n->nsz_len = szlen2(szn);
    n->n_sz = szalloc( n->nsz_len);
    for( i0=0, sp0 = szn; *sp0; ++sp0, ++i0)
        n->n_sz[i0] = *sp0;
    n->n_sz[i0++] = 0;

    n->n_base = base;

    mpz_init(n->n);
    n->n32 = baseX_sz(n->n_sz, n->n_base );
    mpz_set_ui( n->n, n->n32);

    return(n);
}

_ul szlen_to(_c* sz, _c to)
{
    _ul i0;
    _c* sp0, *sp1;
    for( sp0 = sz; *sp0 != to && *sp0; ++sp0)i0++;
    return(i0);
}

onyxword* onyxword_destruct2( onyxword* w)
{
    onyxword* nex = w->next ? : 0;
    free(w->word);
    free(w);
    return(nex);
}

_v onyxword_destruct_set( onyxword* wObjs )
{
    onyxword* w0=0;
    onyxword* n0 = 0;
    for( w0 = wObjs; w0->next;  ){
        n0 = onyxword_destruct2( w0 );
        w0 = n0;
    }
}

onyxbase* onyxbase_destruct( onyxbase* db )
{
    onyxbase* nex = db->next ? : 0;
    onyxword_destruct_set(db->wordset);
    free(db);
    return(nex);
}

_v onyxbase_destruct_objects(onyxbase* dbObj)
{
    onyxbase* db0, *db1;
    db0 = db1 = 0;
    
    for( db0 = dbObj; db0->next; ){
        db1 = onyxbase_destruct(db0);
        db0 = db1;
    }
}

_ul onyxbase_destruct_set(onyxbase* dbObj)
{
    onyxbase* dbp0, *dbp1, *dbp2;
    onyxword* wp0, *wp1, *wp2;
    for( dbp0 = dbObj; dbp0->next; ){

        dbp1 = dbp0;
        dbp0 = dbp0->next;

        for( wp0 = dbp1->wordset; wp0->next; ){
            wp1 = wp0->next;

            if(wp0->word)
                free(wp0->word);
            free( wp0 );

            wp0 = wp1;

        }

        if(dbp1->nam)
            free(dbp1->nam);
        if(dbp1->dat)
            free(dbp1->dat);

        free(dbp1);
    }
    return(0);
}

onyxword* onyxword_copy( onyxword* wObj_a )
{
    _c* sp0;
    _ul i0, i1;
    onyxword* wObj_b = onyxword_new();

    wObj_b->wlen = wObj_a->wlen ? : slen( wObj_a->word );
    wObj_b->word = szalloc( wObj_b->wlen );

    for( sp0 = wObj_a->word, i0 = 0; *sp0; ++sp0, ++i0)
        wObj_b->word[ i0 ] = *sp0;

    wObj_b->word[ wObj_b->wlen ] = 0;

    return(wObj_b);

}


//create duplicate wordset and copy all data from original and return the head pointer to new object...
onyxword* onyxword_copy_wordset(onyxword* wordset)
{
    _ul i0, i1;
    _c* sp0, *sp1;
    onyxword* wObj, *wp0, *wp1;


    wObj = onyxword_new();

    wObj->index = wordset->index;

    for(wp0 = wordset, wp1 = wObj; wp0->next; wp0 = wp0->next ){

        wp1->wlen = wp0->wlen;
        wp1->word = szalloc( wp1->wlen );
        for( i0=0, sp0 = wp0->word; *sp0; ++sp0, ++i0)
            wp1->word[ i0] = *sp0;
        wp1->word[i0++]=0;


        wp1->archetype = wp0->archetype;
        wp1->subtype = wp0->subtype;

        if(wp0->func)
            wp1->func = wp0->func;
        if(wp0->func_)
            wp1->func_ = wp0->func_;

        wp1->next = onyxword_new();


        wp1->next->index = wp1->index + 1;
        wp1->next->master = wObj;
        wp1->next->head = wObj;
        wp1->next->prev = wp1;
        wp1->next->next = 0;

        wp1 = wp1->next;


    }
    return(wObj);
}

onyxword* onyxword_copy_rawsz(_c* rsz, _c to )
{
    _c* sp0, *sp1;
    _ul i0, i1;
    onyxword* wObj_b = onyxword_new();
    onyxword* wp0 = wObj_b;

    wObj_b->wlen = szlen_to(rsz, to);
    wObj_b->word = szalloc( wObj_b->wlen );

    for( sp0 = rsz ; *sp0 != to && *sp0; ){

        for( ; *sp0 != ' '; ++sp0);
        sp0++;

        for( sp1 = sp0, i0=0; *sp0 != ' ' && *sp0 != to; ++sp0)++i0;

        wp0->wlen = i0 + 1;
        wp0->word = szalloc(wp0->wlen);

        for( i0=0; sp1 != sp0; ++sp1, ++i0)
            wp0->word[i0] = *sp1;
        wp0->word[ wp0->wlen ] = 0;

        wp0->next = onyxword_new();
        wp0->next->index = wp0->index + 1;
        onyxword_setptrs(wp0->next, wObj_b, wp0, 0);
        wp0 = wp0->next;

    }


    return(wObj_b);
}

onyxword* onyxword_copy_rawsz_set(_i cnt, _c** rsz )
{
    _ul i0, i1;
    _c* t0, *t1;
    onyxword* wObj;
    onyxword* wp0, *wp1;

    wObj = onyxword_new();
    wp0 = wObj;

    wp0->index = 0;
    wp0->archetype = 1;
    wp0->subtype = 1;

    for( i0=0; i0 < cnt; ++i0){

        wp0->wlen = szlen2( rsz[i0] );
        wp0->word = szalloc(wp0->wlen);
        for( i1=0, t0 = rsz[i0]; *t0; ++t0, ++i1)
            wp0->word[i1] = *t0;
        wp0->word[i1++] = 0;

        wp0->next = onyxword_new();
        onyxword_setmaster(wp0->next, wObj, 0);
        onyxword_setptrs(wp0->next, wObj, wp0, 0);

        wp0->next->archetype = wp0->archetype + 1;
        wp0->next->subtype = 1;

        wp0 = wp0->next;
    }
    return(wObj);
}

_ul onyxword_getcase( _c ch )
{
    _ul ch_case = 0;
    _c* sp0_L = 0;
    _c* sp0_U = 0;

    for( sp0_L = ch_indx_lwr[0], sp0_U = ch_indx_upr[1] ; *sp0_L && *sp0_U; ++sp0_L, ++sp0_U ){

        if( ch == *sp0_L ){
            ch_case = 1;
            break;
        } else if( ch == *sp0_U ){
            ch_case = 2;
            break;
        }
    }
    return(ch_case);
}

onyxword* onyxword_copy_ci( onyxword* wObj_a )
{
    _c* sp0;
    _ul i0, i1;
    onyxword* wObj_b = onyxword_new();

    wObj_b->wlen = wObj_a->wlen ? : slen( wObj_a->word );
    wObj_b->word = szalloc( wObj_b->wlen );

    for( sp0 = wObj_a->word, i0 = 0; *sp0; ++sp0, ++i0){
        wObj_b->word[ i0 ] = *sp0;
    }

    wObj_b->word[ wObj_b->wlen ] = 0;

    return(wObj_b);

}


_ul isnum_deci(_c ch)
{
    _ul i0;
    _ul isnum=0;
    for( i0=0; i0 < 10; ++i0){
        if( ch == n_sz[i0]){
            isnum=1;
            break;
        }else{
            isnum=0;
        }
    }
    return(isnum);
}

_ul isnum_baseN(_c ch, _ul base)
{
    _ul i0;
    _ul isnum=0;
    for( i0=0; i0 < base; ++i0){
        if( ch == n_sz[i0]){
            isnum=1;
            break;
        }else{
            isnum=0;
        }
    }
    return(isnum);
}

_ul onyxword_isnum_deci( onyxword* wObj )
{
    _ul is_num=0;
    _c* sp0 = 0;
    for( sp0 = wObj->word; *sp0; ++sp0 ){
        is_num = isnum_deci( *sp0 );
        if(!is_num)
            break;
    }
    return(is_num);
}

_ul num_getbase(_c* num_sz )
{

    return(0);
}

onyxword* onyxword_parse_at(_c* at, _c dlm)
{
    _ul i0, i1;
    _c* sp0, *sp1;
    onyxword* w = onyxword_new();

    for(sp0 = at, i0=0; *sp0 != dlm && *sp0; ++sp0)++i0;

    w->wlen = i0+1;
    w->word = szalloc( w->wlen );

    for( sp1 = at, i1=0; sp1 != sp0 && i1 <= i0; ++sp1, ++i1)
        w->word[i1] = *sp1;

    w->word[i1++]=0;

    w->szptr_orig = at;
    w->szptr_next = sp0++;

    return(w);
}

onyxword* onyxword_parse_sz(_c* at, _c dlm)
{
    _ul i0, i1;
    _c* spA, *sp0, *sp1;

    onyxword* w0bj = onyxword_new();
    onyxword* w = w0bj;

    onyxword_setptrs(w, w0bj, 0, 0);

    for ( spA = at; *spA; ) {

        for(sp0 = spA, i0=0; *sp0 != dlm && *sp0; ++sp0)++i0;

        w->wlen = i0+1;
        w->word = szalloc( w->wlen );

        for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1, ++i1)
            w->word[i1] = *sp1;

        w->word[i1++]=0;

        w->next = onyxword_new();
        onyxword_setptrs( w->next, w0bj, w, 0 );
        w->next->word = 0; w->next->wlen = 0;
        w->next->index = w->index + 1;
        w = w->next;

        spA = sp1;
        spA++;

    }

    return(w0bj);
}

// eliminates all newline occurances from within parsed blocks
onyxword* onyxword_parse_sz2(_c* at, _c dlm)
{
    _ul i0, i1;
    _c* spA, *sp0, *sp1;

    onyxword* w0bj = onyxword_new();
    onyxword* w = w0bj;

    onyxword_setptrs(w, w0bj, 0, 0);

    for ( spA = at; *spA; ) {

        for(sp0 = spA, i0=0; *sp0 != dlm && *sp0; ++sp0)++i0;

        w->wlen = i0+1;
        w->word = szalloc( w->wlen );

        for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){
            if(*sp1 != '\n'){
                w->word[i1] = *sp1;
                i1++;
            }
        }

        w->word[i1++]=0;

        w->next = onyxword_new();
        onyxword_setptrs( w->next, w0bj, w, 0 );
        w->next->word = 0; w->next->wlen = 0;
        w->next->index = w->index + 1;
        w = w->next;

        spA = sp1;
        spA++;

    }

    return(w0bj);
}


//excludes the copying of any given list (string) of
//characters to parsed text blocks.
onyxword* onyxword_parse_sz3(_c* at, _c dlm, _c* excl)
{
    _ul i0, i1;
    _c* spA, *sp0, *sp1;

    _c* exclP = 0;
    _i exclude = 0;

    onyxword* w0bj = onyxword_new();
    onyxword* w = w0bj;

    onyxword_setptrs(w, w0bj, 0, 0);

    for ( spA = at; *spA; ) {

        for(sp0 = spA, i0=0; *sp0 != dlm && *sp0; ++sp0)++i0;

        w->wlen = i0+1;
        w->word = szalloc( w->wlen );

        for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){

            for(exclude = 0, exclP = excl; *exclP; ++exclP){
                if(*sp1 == *exclP){
                    exclude = 1;
                    break;
                } else {
                    exclude = 0;
                }
            }
            if(!exclude){
                w->word[i1] = *sp1;
                i1++;
            }
        }
        w->word[i1++]=0;

        w->next = onyxword_new();
        onyxword_setptrs( w->next, w0bj, w, 0 );
        w->next->word = 0; w->next->wlen = 0;
        w->next->index = w->index + 1;
        w = w->next;

        spA = sp1;
        spA++;

    }
    return(w0bj);
}

/*    same thing as xx_sz3(..) but with exclude-on-copy symbolics..
 */
onyxword* onyxword_parse_sz3b(_c* at, _c dlm, _c* excl)
{
    _ul i0, i1;
    _c* spA, *sp0, *sp1;

    _c exch=0;
    _c* exclP = 0;
    _i exclude = 0;

    onyxword* w0bj = onyxword_new();
    onyxword* w = w0bj;

    onyxword_setptrs(w, w0bj, 0, 0);

    for ( spA = at; *spA; ) {

        for(sp0 = spA, i0=0; *sp0 != dlm && *sp0; ++sp0)++i0;

        w->wlen = i0+1;
        w->word = szalloc( w->wlen );

        for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){

            for(exclude = 0, exclP = excl; *exclP; ++exclP){

                exch = (*exclP == 'n') ? ('\n') :
                    ( (*exclP == 'r') ? ('\r') :
                        ( (*exclP == 's') ? ( ' ') :
                            ( *exclP) ) );

                if(*sp1 == exch){
                    exclude = 1;
                    break;
                } else {
                    exclude = 0;
                }
            }
            if(!exclude){
                w->word[i1] = *sp1;
                i1++;
            }
        }
        w->word[i1++]=0;

        w->next = onyxword_new();
        onyxword_setptrs( w->next, w0bj, w, 0 );
        w->next->word = 0; w->next->wlen = 0;
        w->next->index = w->index + 1;
        w = w->next;

        spA = sp1;
        spA++;

    }
    return(w0bj);
}

// parse text into blocs and create the linked list from
// thus, chaining one to the next in order of consequence
// as they appear in the string "at". exclude excl characters...
onyxword* onyxword_parse_sz4(_c* at, _c* dlm, _c* excl)
{
    _ul i0, i1;
    _c* spA, *sp0, *sp1;

    _ul dlm_len = slen(dlm);

    _c* exclP = 0;
    _i exclude = 0;

    onyxword* w0bj = onyxword_new();
    onyxword* w = w0bj;

    onyxword_setptrs(w, w0bj, 0, 0);

    for ( spA = at; *spA; ) {

        for(sp0 = spA, i0=0; !dlmchk(sp0, dlm) && *sp0; ++sp0)++i0;

        w->wlen = i0+1;
        w->word = szalloc( w->wlen );

        for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){

            for(exclude = 0, exclP = excl; *exclP; ++exclP){
                if(*sp1 == *exclP){
                    exclude = 1;
                    break;
                } else {
                    exclude = 0;
                }
            }
            if(!exclude){
                w->word[i1] = *sp1;
                i1++;
            }
        }
        w->word[i1++]=0;

        w->next = onyxword_new();
        onyxword_setptrs( w->next, w0bj, w, 0 );
        w->next->word = 0; w->next->wlen = 0;
        w->next->index = w->index + 1;
        w = w->next;

        sp1 += (dlm_len-1);
        spA = sp1;
        spA++;

    }
    return(w0bj);
}

/*   same thing as xx_sz4(..) except that this one parses all data as-is,
 *   excluding nothing. perfect for top-level parsing for things like config
 *   data or internal db-intrinsic symbols/delimeters/tokens to be later used
 *   on subsequent chunks of data/code.
 */
onyxword* onyxword_parse_sz4b(_c* at, _c* dlm)
{
    _ul i0, i1;
    _c* spA, *sp0, *sp1;

    _ul dlm_len = slen(dlm);

    _c* exclP = 0;
    _i exclude = 0;

    onyxword* w0bj = onyxword_new();
    onyxword* w = w0bj;

    onyxword_setptrs(w, w0bj, 0, 0);

    for ( spA = at; *spA; ) {

        for(sp0 = spA, i0=0; !dlmchk(sp0, dlm) && *sp0; ++sp0)++i0;

        w->wlen = i0+1;
        w->word = szalloc( w->wlen );

        for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){


                w->word[i1] = *sp1;
                i1++;

        }
        w->word[i1++]=0;

        w->next = onyxword_new();
        onyxword_setptrs( w->next, w0bj, w, 0 );
        w->next->word = 0; w->next->wlen = 0;
        w->next->index = w->index + 1;
        w = w->next;

        sp1 += (dlm_len-1);
        spA = sp1;
        spA++;

    }
    return(w0bj);
}

//parse using set of delimeters, storing within each node the index
//of the delimeter caught, for later use as archetype reference for
//that bloc type respective of that delimeter...
onyxword* onyxword_parse_sz5(_c* at, _c** dlm, _ul dlm_c, _c* excl)
{
    _ul i0, i1;
    _c* spA, *sp0, *sp1;

    _ul dlm_len = 0;
    _ul dlm_i0, dlm_mat;

    _c* exclP = 0;
    _i exclude = 0;

    onyxword* w0bj = onyxword_new();
    onyxword* w = w0bj;

    onyxword_setptrs(w, w0bj, 0, 0);

    for ( spA = at; *spA; ) {

        dlm_mat = dlm_i0 = dlm_len = 0;

        for(sp0 = spA, i0=0; !dlm_mat && *sp0; ++sp0){
            for( dlm_i0 = dlm_mat; (dlm_i0 < dlm_c) && !dlm_mat; ++dlm_i0){
                dlm_mat = dlmchk(sp0, dlm[dlm_i0] );
                if(dlm_mat)
                    break;
            }
            if(dlm_mat){
                dlm_len = slen( dlm[dlm_i0] );
                dlm_mat = dlm_i0;
                break;
            }
            ++i0;
        }

        w->wlen = i0+1;
        w->word = szalloc( w->wlen );

        for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){

            for(exclude = 0, exclP = excl; *exclP; ++exclP){
                if(*sp1 == *exclP){
                    exclude = 1;
                    break;
                } else {
                    exclude = 0;
                }
            }
            if(!exclude){
                w->word[i1] = *sp1;
                i1++;
            }
        }
        w->word[i1++]=0;

        w->next = onyxword_new();
        w->next->archetype = dlm_mat;
        onyxword_setptrs( w->next, w0bj, w, 0 );
        w->next->word = 0; w->next->wlen = 0;
        w->next->index = w->index + 1;
        w = w->next;

        sp1 += (dlm_len-1);
        spA = sp1;
        spA++;
    }
    return(w0bj);
}


// just the same as xx_sz5(..) except that our exclude-on-copy function-
// ality accepts symbolic token characters and interprets them follows:
// 'n' == '\n', 's' == ' ', 'r' == '\r',..
onyxword* onyxword_parse_sz5a(_c* at, _c** dlm, _ul dlm_c, _c* excl)
{
    _ul i0, i1;
    _c* spA, *sp0, *sp1;

    _ul dlm_len = 0;
    _ul dlm_i0, dlm_mat;

    _c exch;
    _c* exclP = 0;
    _i exclude = 0;

    onyxword* w0bj = onyxword_new();
    onyxword* w = w0bj;

    onyxword_setptrs(w, w0bj, 0, 0);

    for ( spA = at; *spA; ) {

        dlm_mat = dlm_i0 = dlm_len = 0;

        for(sp0 = spA, i0=0; !dlm_mat && *sp0; ++sp0){
            for( dlm_i0 = dlm_mat; (dlm_i0 < dlm_c) && !dlm_mat; ++dlm_i0){
                dlm_mat = dlmchk(sp0, dlm[dlm_i0] );
                if(dlm_mat)
                    break;
            }
            if(dlm_mat){
                dlm_len = slen( dlm[dlm_i0] );
                dlm_mat = dlm_i0;
                break;
            }
            ++i0;
        }

        w->wlen = i0+1;
        w->word = szalloc( w->wlen );

        for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){

            for(exclude = 0, exclP = excl, exch = 0; *exclP; ++exclP){

                exch = (*exclP == 'n') ? ('\n') : //grotesque if elegant overuse of
                    ( (*exclP == 'r') ? ('\r') :  //nested conditional expressions..
                        ( (*exclP == 's') ? ( ' ') :
                            ( *exclP) ) );


                if( *sp1 == exch ){
                    exclude = 1;
                    break;
                } else {
                    exclude = 0;
                }
                exch = 0;
            }
            if(!exclude){
                w->word[i1] = *sp1;
                i1++;
            }
        }
        w->word[i1++]=0;

        w->next = onyxword_new();
        w->next->archetype = dlm_mat;
        onyxword_setptrs( w->next, w0bj, w, 0 );
        w->next->word = 0; w->next->wlen = 0;
        w->next->index = w->index + 1;
        w = w->next;

        sp1 += (dlm_len-1);
        spA = sp1;
        spA++;
    }
    return(w0bj);
}


/* this variant does the same as xx_sz5a except that it offers control over character
*  repeats, also disallow copied strings from beginning with reap_chz
* */
onyxword* onyxword_parse_sz5b(_c* at, _c** dlm, _ul dlm_c, _c* excl, _c* reap_chz, _ul reap_lim)
{
   _ul i0, i1;
   _c* spA, *sp0, *sp1;

   _ul r0, r1;
   _c* spR, *spR2;

   _ul dlm_len = 0;
   _ul dlm_i0, dlm_mat;

   _c reap_ch, last_ch;
   _c* reapP=0;
   _ul reap_n, is_lim;

   _c exch;
   _c* exclP = 0;
   _i exclude = 0;

   onyxword* w0bj = onyxword_new();
   onyxword* w = w0bj;

   onyxword_setptrs(w, w0bj, 0, 0);

   for ( spA = at, last_ch = 0; *spA; ) {

       dlm_mat = dlm_i0 = dlm_len = 0;

       for(sp0 = spA, i0=0; !dlm_mat && *sp0; ++sp0){
           for( dlm_i0 = dlm_mat; (dlm_i0 < dlm_c) && !dlm_mat; ++dlm_i0){
               dlm_mat = dlmchk(sp0, dlm[dlm_i0] );
               if(dlm_mat)
                   break;
           }
           if(dlm_mat){
               dlm_len = slen( dlm[dlm_i0] );
               dlm_mat = dlm_i0;
               break;
           }
           ++i0;
       }

       w->wlen = i0+1;
       w->word = szalloc( w->wlen );

       for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){

           for(exclude = 0, exclP = excl, exch = 0; *exclP; ++exclP){

               exch = (*exclP == 'n') ? ('\n') :
                   ( (*exclP == 'r') ? ('\r') :
                       ( (*exclP == 's') ? ( ' ') :
                           ( *exclP) ) );


               if( *sp1 == exch ){
                   exclude = 1;
                   break;
               } else {
                   exclude = 0;
               }
               exch = 0;
           }

           for( reapP = reap_chz; *reapP; ++reapP){

               reap_ch = (*reapP == 'n') ? ('\n') :
                   ( (*reapP == 'r') ? ('\r') :
                       ( (*reapP == 's') ? ( ' ') :
                           ( *reapP) ) );

               if( *sp1 == reap_ch ){

                   if( i1 == 0 ){
                       for(; *sp1 == reap_ch; ++sp1);
                   }
                   if( *(sp1+1) == reap_ch ){
                       for( spR = sp1, r0=0 ; *spR == reap_ch; ++spR, ++r0);

                       if( r0 > reap_lim){

                           for( r1=0; r1 < reap_lim; ++r1){
                               w->word[i1] = reap_ch;
                               i1++;
                           }
                           sp1 += r0;
                       }
                   }
                   break;
               }
           }

           if( !exclude ){

               w->word[i1] = *sp1;
               i1++;
           }
       }
       last_ch = *sp1;

       w->word[i1++]=0;

       w->next = onyxword_new();
       w->next->archetype = dlm_mat;
       onyxword_setptrs( w->next, w0bj, w, 0 );
       w->next->word = 0; w->next->wlen = 0;
       w->next->index = w->index + 1;
       w = w->next;

       sp1 += (dlm_len-1);
       spA = sp1;
       spA++;
   }
   return(w0bj);
}

/*  same as xx_sz5b(..) except that it uses a pre-parsed onyxword list of delimeters rather than
 *  a two-dimensional char* array ..
 */
onyxword* onyxword_parse_sz5c(_c* at, onyxword* dlm, _ul dlm_c, _c* excl, _c* reap_chz, _ul reap_lim)
{
   _ul i0, i1;
   _c* spA, *sp0, *sp1;

   _ul r0, r1;
   _c* spR, *spR2;

   _ul dlm_len = 0;
   _ul dlm_i0, dlm_mat;
   onyxword* dWp0=0;

   _c reap_ch, last_ch;
   _c* reapP=0;
   _ul reap_n, is_lim;

   _c exch;
   _c* exclP = 0;
   _i exclude = 0;

   onyxword* w0bj = onyxword_new();
   onyxword* w = w0bj;

   onyxword_setptrs(w, w0bj, 0, 0);

   for ( spA = at, last_ch = 0; *spA; ) {

       dlm_mat = dlm_i0 = dlm_len = 0;

       for(sp0 = spA, i0=0; !dlm_mat && *sp0; ++sp0){

           for( dWp0 = dlm; dWp0->next; dWp0 = dWp0->next ){
               dlm_mat = dlmchk(sp0, dWp0->word );
               if(dlm_mat)
                   break;
           }

           if(dlm_mat){
               dlm_len = slen( dWp0->word );
               dlm_mat = dWp0->index;
               break;
           }
           ++i0;
       }

       w->wlen = i0+1;
       w->word = szalloc( w->wlen );

       for( sp1 = spA, i1=0; sp1 != sp0 && i1 <= i0; ++sp1){

           for(exclude = 0, exclP = excl, exch = 0; *exclP; ++exclP){

               exch = (*exclP == 'n') ? ('\n') :
                   ( (*exclP == 'r') ? ('\r') :
                       ( (*exclP == 's') ? ( ' ') :
                           ( *exclP) ) );


               if( *sp1 == exch ){
                   exclude = 1;
                   break;
               } else {
                   exclude = 0;
               }
               exch = 0;
           }

           for( reapP = reap_chz; *reapP; ++reapP){

               reap_ch = (*reapP == 'n') ? ('\n') :
                   ( (*reapP == 'r') ? ('\r') :
                       ( (*reapP == 's') ? ( ' ') :
                           ( *reapP) ) );

               if( *sp1 == reap_ch ){

                   if( i1 == 0 ){
                       for(; *sp1 == reap_ch; ++sp1);
                   }
                   if( *(sp1+1) == reap_ch ){
                       for( spR = sp1, r0=0 ; *spR == reap_ch; ++spR, ++r0);

                       if( r0 > reap_lim){

                           for( r1=0; r1 < reap_lim; ++r1){
                               w->word[i1] = reap_ch;
                               i1++;
                           }
                           sp1 += r0;
                       }
                   }
                   break;
               }
           }

           if( !exclude ){

               w->word[i1] = *sp1;
               i1++;
           }
       }
       last_ch = *sp1;

       w->word[i1++]=0;

       w->next = onyxword_new();
       w->next->archetype = dlm_mat;
       onyxword_setptrs( w->next, w0bj, w, 0 );
       w->next->word = 0; w->next->wlen = 0;
       w->next->index = w->index + 1;
       w = w->next;

       sp1 += (dlm_len-1);
       spA = sp1;
       spA++;
   }
   return(w0bj);
}


// break lines and construct database from each of thier wordsets from raw scratch data...
onyxbase* onyxbase_lineset( _c* dat_in )
{
    _ul i0, i1, i2;
    _c* sp0, *sp1, *sp2;
    onyxbase* dbObj;
    onyxbase* db0, *db1, *db2;
    onyxword* wp0, *wp1, *wp2;

    _i line_good = 0;

    db0 = onyxbase_new();
    db0->index = 0;

    dbObj = db0;
    dbObj->words = 0;

    for( sp0 = dat_in; *sp0; ){

        if( sp0[0] == ' ' ){

            line_good = 0;
            for( sp1 = sp0; *sp1 && *sp1 != '\n'; ++sp1){
                if(*sp1 != ' '){
                    line_good = 1;
                    break;
                }
            }

        } else if( sp0[0] == '\0' || sp0[0] == 0 || sp0[0] == '\n'){
            line_good = 0;
        } else {
            line_good = 1;
        }

        if( line_good) {
            db0->wordset = onyxword_new();
            onyxword_setptrs( db0->wordset, db0->wordset, 0, 0);
            onyxword_setmaster( db0->wordset, db0->wordset, 0);
            db0->wordset->index = 0;
            db0->words = 0;
            wp0 = db0->wordset;

            for(sp1 = sp0, i0=0 ; *sp0 && *sp0 != '\n'; ++sp0 )i0++;
            db0->datlen = i0;
            db0->dat = szalloc( i0+1 );

            for( i0=0; sp1 != sp0; ++sp1, ++i0)
                db0->dat[i0] = *sp1;

            db0->dat[i0]=0;

            if( db0->dat[0])

            for( sp1 = db0->dat; *sp1; ){

                if( sp1[0] == '\0' || sp1[0] == 0 )
                    sp1++;

                for( ; *sp1 == ' '; ++sp1);

                for( sp2 = sp1, i0=0; *sp1 && *sp1 != ' '; ++sp1)++i0;
                wp0->wlen = i0;
                wp0->word = szalloc( wp0->wlen );
                for( i0=0; sp2 != sp1; ++sp2, ++i0)
                    wp0->word[i0] = *sp2;
                wp0->word[i0]=0;

                wp0->next = onyxword_new();
                onyxword_setptrs( wp0->next, dbObj->wordset, wp0, 0 );
                onyxword_setmaster( wp0->next, dbObj->wordset, 0);
                wp0->next->index = wp0->index + 1;
                wp0 = wp0->next;
                db0->words++;
            }
            db0->next = onyxbase_new();
            onyxbase_setptrs( db0->next, dbObj, db0, 0);
            onyxbase_setmaster( db0->next, dbObj, 0);
            db0->next->index = db0->index + 1;
            db0 = db0->next;

            dbObj->words++;
        }
        sp0++;
    }
    return(dbObj);
}

// break lines and construct database from each of thier wordsets from raw scratch data...
onyxbase* onyxbase_lineset2( _c* dat_in, _c line_dlm )
{
    _ul i0, i1, i2;
    _c* sp0, *sp1, *sp2;
    onyxbase* dbObj;
    onyxbase* db0, *db1, *db2;
    onyxword* wp0, *wp1, *wp2;

    _i line_good = 0;

    db0 = onyxbase_new();
    db0->index = 0;

    dbObj = db0;
    dbObj->words = 0;

    for( sp0 = dat_in; *sp0; ){

        if( sp0[0] == ' ' ){

            line_good = 0;
            for( sp1 = sp0; *sp1 && *sp1 != line_dlm; ++sp1){
                if(*sp1 != ' '){
                    line_good = 1;
                    break;
                }
            }

        } else if( sp0[0] == '\0' || sp0[0] == 0 || sp0[0] == line_dlm){
            line_good = 0;
        } else {
            line_good = 1;
        }

        if( line_good) {
            db0->wordset = onyxword_new();
            onyxword_setptrs( db0->wordset, db0->wordset, 0, 0);
            onyxword_setmaster( db0->wordset, db0->wordset, 0);
            db0->wordset->index = 0;
            db0->words = 0;
            wp0 = db0->wordset;

            for(sp1 = sp0, i0=0 ; *sp0 && *sp0 != line_dlm; ++sp0 )i0++;
            db0->datlen = i0;
            db0->dat = szalloc( i0+1 );

            for( i0=0; sp1 != sp0; ++sp1, ++i0)
                db0->dat[i0] = *sp1;

            db0->dat[i0]=0;

            if( db0->dat[0])

            for( sp1 = db0->dat; *sp1; ){

                if( sp1[0] == '\0' || sp1[0] == 0 )
                    sp1++;

                for( ; *sp1 == ' '; ++sp1);

                for( sp2 = sp1, i0=0; *sp1 && *sp1 != ' '; ++sp1)++i0;
                wp0->wlen = i0;
                wp0->word = szalloc( wp0->wlen );
                for( i0=0; sp2 != sp1; ++sp2, ++i0)
                    wp0->word[i0] = *sp2;
                wp0->word[i0]=0;

                wp0->next = onyxword_new();
                onyxword_setptrs( wp0->next, dbObj->wordset, wp0, 0 );
                onyxword_setmaster( wp0->next, dbObj->wordset, 0);
                wp0->next->index = wp0->index + 1;
                wp0 = wp0->next;
                db0->words++;
            }
            db0->next = onyxbase_new();
            onyxbase_setptrs( db0->next, dbObj, db0, 0);
            onyxbase_setmaster( db0->next, dbObj, 0);
            db0->next->index = db0->index + 1;
            db0 = db0->next;

            dbObj->words++;
        }
        sp0++;
    }
    return(dbObj);
}

onyxbase* onyxbase_x(_c* line, onyxword* wp_x )
{
    _ul i0, i1;
    _ul len0;
    _ul match0;

    onyxbase* dbObj;
    onyxbase* db0, *db1, *db2;
    onyxword* wp0, *wp1, *wp2;

    _c* sp0, *sp1, *sp2;

    dbObj = onyxbase_new();
    db0 = dbObj;

    for( sp0 = line; *sp0; ){

        for( wp0 = wp_x; wp0->next; wp0 = wp0->next ){
            match0 = ucstrcmp3( sp0, wp0->word );
            if(match0)
                break;
        }
        if(match0){

            db0->next = onyxbase_new();
            if(wp0->archetype == 1){
                db0->archetype = 1;
            }else if(wp0->archetype == 2){

            }else if(wp0->archetype == 3){

            }else if(wp0->archetype == 4){

            }

        }

        sp0++;
    }


}

onyxbase* onyxbase_db( onyxword* wset, _c** dlm_set, _ul dlm_c, _c* excl_set)
{
    _ul i0, i1;

    _c* dlmSz;
    _c* sp0, *sp1, *sp2;

    onyxword* wObj = 0;
    onyxword* w0, *w1, *w2;
    onyxbase* db = onyxbase_new();
    onyxbase* dbObj = db;

    for( wObj = wset, i0 = 0; wObj->next; wObj = wObj->next, ++i0){

        dbObj->wordset = onyxword_parse_sz5b(wObj->word, dlm_set, dlm_c, excl_set, rep_sz0, 1 );

        dbObj->next = onyxbase_new();

        dbObj->next->index = dbObj->index + 1;

        onyxbase_setptrs( dbObj->next, db, dbObj, 0);

        dbObj = dbObj->next;
    }

    return( db );
}

onyxbase* onyxbase_db2( _c* db_fnam , _c* db_symtable_fnam, _c dlm_c )
{
    _c* fin_db = fget_data(db_fnam);
    _c* fin_dlm = fget_data(db_symtable_fnam);

    onyxbase* db_obj;
    onyxword* dlm_w;

    onyxbase* db0, *db1, *db2;
    onyxword* wp0, *wp1, *wp2;
    db0 = db1 = db2 = 0;
    wp0 = wp1 = wp2 = 0;

    dlm_w = onyxword_parse_sz2(fin_dlm, dlm_c);



}

onyxquery* onyxbase_query(_c* term, onyxbase* db )
{
    _ul i0, i1;

    _i match;

    onyxquery* qT0, *qT1;
    onyxbase* dP0, *dP1;
    onyxword* wP0, *wP1;

    qT0 = onyxquery_new();
    qT0->dbObj = onyxbase_new();
    qT0->dbObj->wordset = onyxword_new();

    qT0->index = qT0->dbObj->index = qT0->dbObj->wordset->index = 0;

    qT0->dbEntry = qT0->dbObj->wordset;

    qT0->dbObj->words = 0;

    match = 0;

    qT0->query_hits = 0;
    qT0->query_index_current = 0;
    for( dP0 = db; dP0->next; dP0 = dP0->next ){

        for( wP0 = wP1 = dP0->wordset; wP0->next; wP0 = wP0->next ){
            match = dlmchk( term, wP0->word );

            if(match){
                qT0->dbObj->wordset->next = wP0;

                qT0->dbObj->wordset->next->ref[0] = wP0->index;

                qT0->dbObj->wordset->next->index = qT0->dbObj->wordset->index + 1;

                qT0->dbObj->wordset = qT0->dbObj->wordset->next;

                qT0->dbObj->words++;

                qT0->query_hits++;
            }

            match = 0;
        }
    }
    qT0->dbObj->wordset = qT0->dbEntry;

    return(qT0);
}

onyxquery* onyxbase_query2(_c* term, onyxbase* db )
{
    _ul i0, i1;

    _i match;

    onyxquery* qT0, *qT1, *qT2;;
    onyxbase* dP0, *dP1;
    onyxword* wP0, *wP1;

    qT0 = onyxquery_new();
    qT0->query_hits = 0;
    qT0->query_index_current = 0;

    qT1 = qT2 = qT0;

    match = 0;
    for( dP0 = db; dP0->next; dP0 = dP0->next ){

        for( wP0 = wP1 = dP0->wordset; wP0->next; wP0 = wP0->next ){
            match = dlmchk( term, wP0->word );

            if(match){
                qT0->dbObj = dP0;
                qT0->dbEntry = onyxword_copy(wP1);

                qT0->next = onyxquery_new();
                qT0->next->index = qT0->index + 1;
                onyxquery_setptrs(qT0->next, qT1, qT0, 0);

                qT0 = qT0->next;

                qT1->query_hits++;
            }

            match = 0;
        }
    }

    return(qT1);
}

onyxbase* onyxbase_query3( _c* term, onyxbase* db)
{
    _ul i0;
    _c* sp0;
    onyxbase* qObj,  *qt0;
    onyxbase* db0, *db1;
    onyxword* wp0, *wp1;

    qObj = onyxbase_new();
    qt0 = qObj;

    _i match, matches_found;


    qObj->wordset = onyxword_copy_wordset(db->wordset);

    qObj->index = db->index;

    matches_found = 0;
    for( db0 = db ; db0->next; db0 = db0->next ){

        for( wp0 = wp1 = db0->wordset; wp0->next; wp0 = wp0->next){
            match = ucstrcmp3(term, wp0->word );
            if(match){

                if(!matches_found)
                    matches_found = 1;

                qt0->next = onyxbase_new();
                qt0->next->wordset = onyxword_copy_wordset(wp1);

                qt0->next->index = db0->index;

                qt0->next->master = qt0->next->head = qObj;
                qt0->next->prev = qt0;
                qt0->next->next = 0;

                qt0 = qt0->next;
                break;
            }
        }
    }
    if(!matches_found){

        onyxword_destruct2(qObj->wordset);
        free(qObj);

        qObj = 0;
    }
    return(qObj);
}

onyxbase* onyxshell_exec(_c* sh_code )
{
    _c* din=0;
    onyxbase* dbObj=0;
    fput_data("./onyx-tmp.sh", sh_code);
    system("chmod +rwx ./onyx-tmp.sh");
    system("./onyx-tmp.sh > ./onyx-tmp.out.oxdb");

    din = fget_data( "./onyx-tmp.out.oxdb");

    dbObj = onyxbase_lineset(din);

    free(din);
    return(dbObj);
}

onyxbase* onyxshell_execlines(_c* sh_code )
{
    _c* din=0;
    onyxbase* dbObj=0;
    onyxbase* dbLines=0;
    
    dbLines = onyxbase_lineset(sh_code);
    
    for( dbObj=dbLines; dbObj->next; dbObj=dbObj->next ){
    
        fput_data("./onyx-tmp.sh", dbObj->dat );
        system("chmod +rwx ./onyx-tmp.sh");
        system("./onyx-tmp.sh > ./onyx-tmp.out.oxdb");
    
        din = fget_data("./onyx-tmp.out.oxdb");
        
        dbObj->stem = onyxbase_lineset(din);
        
        free(din);
    
    }
    
    return(dbLines);
}

/* where sz_src is the "leave-off" point in the parsing chain. this same point will be stored
 * into the onyxword (ret) object so that a cycle of chain-parsing may continue from hence..
 */
onyxword* onyx_parse_link(_c* sz_src, onyxword* dlm_set )
{
    _ul i0, i1;
    _c* sp0, *sp1;
    onyxword* wObj_r, *wObj_p;
    onyxword* wUtl0;

    _ul match=0;

    wObj_r = 0; //onyxword_new();

    for( sp0 = sz_src; *sp0; ){

        for( wUtl0 = dlm_set; wUtl0->next; wUtl0 = wUtl0->next ){
            match = ucstrcmp3( sp0, wUtl0->word );
            if(match)
                break;
        }
        if(match)
            break;
        else
            sp0++;
    }
    if(match){
        wObj_r = onyxword_new();
        wObj_p = wObj_r;

        onyxword_setmaster(wObj_r, wObj_r, 0);
        onyxword_setptrs(wObj_r, wObj_r, 0, 0);

        for( ; *sp0 != ' '; ++sp0);
        sp0++;

        sp1 = sp0;

        for(i0=0; *sp0 != ' '; ++sp0)i0++;

        wObj_p->archetype = wUtl0->archetype;
        wObj_p->subtype = wUtl0->subtype;

        wObj_p->wlen = i0 + 1;
        wObj_p->word = szalloc( wObj_p->wlen );
        for( i0=0; sp1 != sp0; ++i0, ++sp1)
            wObj_p->word[ i0 ] = *sp1;
        wObj_p->word[i0++]=0;

        wObj_p->szptr_orig = sz_src;
        wObj_p->szptr_next = sp0++;

    }else{
        sp0++;
    }

    return(wObj_r);
}

onyxbase* onyx_parse_str(_c* str, onyxword* dlms )
{
    onyxbase* dbObj, *db0, *db1;
    onyxword* wObj, *wp0, *wp1;

    wObj = onyx_parse_link(str, dlms);

    for( wp0 = wObj; ; ){
        wp0->next = onyx_parse_link( wp0->szptr_next, dlms );
        if(!wp0->next){
            break;
        }
    }

    return(wObj);
}

_v sz_seek_spc(_c* sp0)
{
    for( ;*sp0 != ' ' && *sp0; ++sp0);
    sp0++;
}

_ul sz_cnt_to_spc(_c* sp0)
{
    _ul i0;
    for( i0=0;*sp0 != ' ' && *sp0; ++sp0)i0++;
    sp0++;
    return(i0);
}
_v sz_seek_spc_or(_c* sp0, _c cor)
{
    for( ;*sp0 != ' ' && *sp0 != cor && *sp0; ++sp0);
    sp0++;
}

/*
 *  onyxObj->funcObj = onyxbase_new();
    onyxObj->varObj = onyxbase_new();
    onyxObj->shellObj = onyxbase_new();
    onyxObj->fileObj = onyxbase_new();
    onyxObj->queryObj = onyxquery_new();
    */

onyxbase* onyxcode_translate_symbolix( onyxcode_object* toplev, _c* query_term )
{
    _ul mat;
    _ul i0, i1;
    _c* sp0, *sp1, *sp2;

    onyxbase* resObj=0;
    onyxbase* dbX[5], *dbY[5], *dbZ[5];
    onyxword* wp0, *wp1, *wp2, *wp3;

    for( mat=0, dbX[0] = toplev->funcObj,
         dbX[1]=toplev->varObj, dbX[2]=toplev->shellObj,
         dbX[3]=toplev->fileObj, dbX[4]=toplev->queryObj;
          !mat; ){

        for( i0 =0, mat=0; i0 < 5; ++i0){
            for( dbY[0] = dbX[i0]; dbY[0]->next; dbY[0] = dbY[0]->next ){
                mat  = ucstrcmp3( query_term, dbY[0]->nam );
                if(mat)
                    break;
            }
            if(mat)
                break;
        }
        if(mat){
            resObj = dbY[0];
            break;
        }

    }

    return(resObj);
}

onyxcode_object* onyxcode_funcparse( _c* code )
{

    _ul i0[3], i1[3];
    _c* sp0, *sp1, *sp2;
    _c* spA, *spB, *spC;
    
    onyxbase* dbObj;
    onyxbase* vObj;
    onyxbase* shObj;
    onyxbase* fileObj;

    onyxquery* qObj;

    onyxcode_object* onyxObj;


    onyxbase* db0, *db1, *db2, *dbX, *dbF;
    onyxword* wp0, *wp1, *wp2, *wpX;
    onyxword* wTmp0, *wTmp1;
    onyxquery* qp0, *qp1;

    _i brack_skip=0;
    _i varnam_match=0;
    _i eval=0;

    onyxObj = onyxcode_object_new();
    
    onyxObj->funcObj = onyxbase_new();
    onyxObj->varObj = onyxbase_new();
    onyxObj->shellObj = onyxbase_new();
    onyxObj->fileObj = onyxbase_new();
    onyxObj->queryObj = onyxquery_new();

    dbObj = onyxObj->funcObj;
    vObj = onyxObj->varObj;
    shObj = onyxObj->shellObj;
    fileObj = onyxObj->fileObj;
    qObj = onyxObj->queryObj;

    dbObj->index=0;
    vObj->index = 0;
    shObj->index = 0;
    fileObj->index = 0;
    qObj->index = 0;
    
    for( sp0 = code, db0 = dbObj, db1=vObj, db2=shObj, dbF=fileObj; *sp0; ){
    
        if( ucstrcmp3( sp0, "func" ) ){
            
            
            for( ; *sp0 != ' ' && *sp0; ++sp0);
            sp0++;                                  //make sure the *sp0 ptr is not on  the ' '
                                                   
            for( sp1 = sp0 , i0[0]=0;*sp0 != ' ' && *sp0 != '{' && *sp0 ; ++sp0) i0[0]++;
            
            db0->namlen = i0[0] + 1;
            db0->nam = szalloc( db0->namlen );
            
            for( i0[1]=0; sp1 != sp0; ++i0[1], ++sp1)
                db0->nam[ i0[1] ] = *sp1;
            db0->nam[i0[1]++] = 0;
            
            for( ; *sp0 != '{'; ++sp0);
            sp0++;
            
            for( i0[0] = 0, sp1 = sp0; *sp0 != '}' && *sp0; ++sp0) i0[0]++;
            
            db0->datlen = i0[0] + 1;
            db0->dat = szalloc( db0->datlen );
            
            for( i0[1]=0; sp1 != sp0; ++i0[1], ++sp1)
                db0->dat[ i0[1] ] = *sp1;
            db0->dat[i0[1]++] = 0;
            
            db0->stem = onyxbase_lineset2(db0->dat, ';');
            
            db0->next = onyxbase_new();
            db0->next->index = db0->index + 1;
            db0->next->next = 0;
            
            db0 = db0->next;
            
        
        }else if( ucstrcmp3( sp0, "var" ) ){//TOD0: fix bug in onyxword_copy_rawsz(.) algo that causes it
                                            //crash upon passing string with only one word element preced-
                                            //ing the semicolon. . .
            for( ;*sp0 != ' ' && *sp0; ++sp0);
            sp0++;for( ;*sp0 != ' ' && *sp0; ++sp0);
            sp0++;
            
            for( sp1 = sp0, i0[0]=0; *sp0 != ' '; ++sp0) i0[0]++;
        
            db1->namlen = i0[0] + 1;
            db1->nam = szalloc( i0[0] );
            
            for( i0[1]=0; sp1 != sp0; ++i0[1], ++sp1)
                db1->nam[i0[1]] = *sp1;
            db1->nam[i0[1]++] = 0;

            sp0++;

            if( ucstrcmp3( sp0, "set" ) ){

                db1->wordset = onyxword_copy_rawsz( sp0, ';');

            }else if( ucstrcmp3( sp0, "eval") ){

                wTmp0 = onyxword_copy_rawsz( sp0, ';' );
                for( dbX = vObj; dbX->next; dbX = dbX->next ){
                    varnam_match = ucstrcmp3( wTmp0->word, dbX->nam );
                    if(varnam_match)
                        break;
                }
                if( varnam_match )
                    eval = ucstrcmp3( wTmp0->next->word, dbX->wordset->word );
                else
                    eval = -1;

                db1->archetype = 2     ;
                db1->subtype = eval;

                varnam_match = 0;

                printf("--!eval:%s\n", eval? ("yes") : ("no") );

            }else if( ucstrcmp3( sp0, "copy") ){

                wTmp0 = onyxword_copy_rawsz( sp0, ';');

                for(dbX = vObj ; dbX->next; dbX = dbX->next){
                    varnam_match = ucstrcmp3( wTmp0->word, dbX->nam );

                    if(varnam_match)
                        break;
                }
                if(varnam_match){
                    db1->wordset = onyxword_copy(dbX->wordset );
                } else {
                    db1->wordset = onyxword_copy(wTmp0);
                    //destroy wTmp0 here!
                }

                varnam_match= 0;
            }else if( ucstrcmp3( sp0, "buf" ) ){

                wTmp0 = onyxword_copy_rawsz(sp0, ';');

                printf(" ..BUF: %s -", wTmp0->word );

                db1->n =  baseX_sz( wTmp0->word, 10 );

                printf("%i\n", db1->n );

                db1->datlen = db1->n + 1;
                db1->dat = szalloc( db1->datlen );

            }
        
            db1->next = onyxbase_new();
            db1->next->index = db1->index + 1;
            db1->next->next = 0;
            
            db1 = db1->next;
        
        }else if( ucstrcmp3( sp0, "shell") ){
        
            for( ;*sp0 != ' ' && *sp0; ++sp0);
            sp0++;
            
            for( sp1 = sp0, i0[0]=0; *sp0 != ' ' && *sp0 != '{'; ++sp0) i0[0]++;
        
            brack_skip = (*sp0 == ' ') ? 1 : 0; // todo later: write algo above to 
                                                // accnt for +1 space at end of nam ..
            db2->namlen = i0[0] + 1;
            db2->nam = szalloc( i0[0] );
            
            for( i0[1]=0; sp1 != sp0; ++i0[1], ++sp1)
                db2->nam[i0[1]] = *sp1;
            db2->nam[i0[1]++] = 0;
            
            
            for( ; *sp0 != '{' && *(sp0+1) != '{'; ++sp0);
            sp0 += brack_skip == 1 ? 3 : 2;
            
            for( i0[0] = 0, sp1 = sp0; *sp0 != '}' && *(sp0+1) != '}' && *sp0; ++sp0) i0[0]++;

            db2->datlen = i0[0] + 1;
            db2->dat = szalloc( db2->datlen );
            
            for( i0[1]=0; sp1 != sp0; ++i0[1], ++sp1)
                db2->dat[ i0[1] ] = *sp1;
            db2->dat[i0[1]++] = 0;

            db2->next = onyxbase_new();
            db2->next->index = db2->index + 1;
            db2->next->next = 0;
            
            db2 = db2->next;
        
        }
        /* Now its time to stop all this giggley little bullshit and get down
           * to the grit and grime of this fucking program. stop your goddamn
           * little laughedy, fiddle-fucking bullshit and try to act serious.
           */
        else if( ucstrcmp3(sp0, "file"  ) ){

            for( ;*sp0 != ' ' && *sp0; ++sp0);
            sp0++;

            for( sp1 = sp0, i0[0]=0; *sp0 != ' '; ++sp0) i0[0]++;

            dbF->namlen = i0[0] + 1;
            dbF->nam = szalloc( dbF->namlen );

            for( i0[1]=0; sp1 != sp0; ++i0[1], ++sp1)
                dbF->nam[i0[1]] = *sp1;
            dbF->nam[i0[1]++] = 0;

            sp0++;

            if( ucstrcmp3(sp0, "read" ) ){

                dbF->dat = fget_data( dbF->nam );
                dbF->datlen = szlen2( dbF->dat );

                printf("    FILE READ[%i]: %s\n",dbF->datlen, dbF->dat);

            }else if( ucstrcmp3( sp0, "write" ) ){

                //spA = sp0 + szlen_to(sp0, ' ');
                //dbF->wordset = onyxword_translate_fromsrc(onyxObj, spA );

            }else if( ucstrcmp3( sp0, "open" ) ){

            }else if( ucstrcmp3( sp0, "pipe" ) ){

            }

        }else if( ucstrcmp3(sp0, "loop" )){


        }else if( ucstrcmp3( sp0, "query" ) ){

        }
    
        sp0++;
    }
    
    
    return(onyxObj);
}

onyxbase* onyxbase_find_obj(onyxcode_object* oObj, _c* onam)
{
    _ul i0, i1;
    onyxbase* dbX[5], *dbP;
    onyxword* wp0, *wp1, *wp2;

    _ul match=0;

    dbX[0] = oObj->funcObj;
    dbX[1] = oObj->varObj;
    dbX[2] = oObj->shellObj;
    dbX[3] = oObj->fileObj;
    dbX[4] = oObj->queryObj->dbObj;

    match = 0;
    for(i0=0; i0 < 5; ++i0){

        for( dbP = dbX[i0]; dbP->next; dbP = dbP->next ){
            match = ucstrcmp3( dbP->nam, onam );
            if(match){
                break;
            }else{
                for(wp0 = dbX[i0]->wordset ; wp0->next ; wp0 = wp0->next){
                    match = ucstrcmp3( wp0->word, onam );
                    if(match)
                        break;
                }
            }
            if(match)
                break;
        }
        if(match)
            break;
    }
    if(!match)
        dbP=0;
    return(dbP);
}

onyxbase* onyxbase_eval_expression(onyxword* expr, onyxword* dlmset )
{
    _ul i0, i1, i2;
    _c* sp0, *sp1, *sp2;

    onyxbase* dbExpr;

    onyxword* wp0, *wp2;
    onyxword* dlm0, *dlm1, *dlm2;

    _ul match=0;

    dbExpr = onyxbase_new();

    for(wp0 = expr; wp0->next; wp0 = wp0->next){

        for( dlm0 = dlmset, match=0; dlm0->next && !match; dlm0 = dlm0->next){
            match = ucstrcmp3( wp0->word, dlm0->word );
            if(match)
                break;
        }
        if(match){
            wp0->next->rstack = dlm0->func_( wp0->next->word );
            wp0->state = 3;
            wp0->next->state = 3;

            wp0->archetype = 2;

            wp0->next->archetype = dlm0->index;



        } else {
            wp0->archetype = 0;
            wp0->next->archetype = 0;
        }

    }

    dbExpr->wordset = expr;

    return(dbExpr);
}

static _v onyxcode_execute(onyxcode_object* cObj )
{

    onyxbase* fObj, *vObj, *shObj;
    onyxword* wp0, *wp1, *wp2;


    for(shObj= cObj->shellObj ; shObj->next; shObj = shObj->next){
        shObj->stem = onyxshell_exec( shObj->dat );
    }

}

onyxword* onyxword_line_exec(onyxword* line)
{
    _ul i0, i1, i2;
    _c* sp0, *sp1, *sp2;
    onyxword* wp0, *wp1, *wp2, *wp3;

    onyxtype* typeObj = onyxtype_new();

    for(wp0 = line ; wp0->next ; wp0=wp0->next){
        sp0 = wp0->word;
        if( sp0[0] == '.' ){

        }else if( sp0[0] == '%' ){

        }

    }

}




/*                |-=\  |=  ||==\     /==\|-=\   |-=\ + \    /-\
                  |~    | \ ||  /    | [ /|^     |^   | |   |--/
                  |___/ |  \||_/      \_/ |      |    | |-_/ \_/
*/
