
/*
#ifndef _unixcube_openware_onyxcore_h_
  #define _unixcube_openware_onyxcore_h_ 1313

#define _CDD_IOBUFFER_BYTELENGTH_ 10000
*/


/* 
    type definitions

*/

typedef struct szseg_t
{
    _ul i0, i1, i2;
    _ul seg, off;
    _ul len;
    
    struct szseg_t* corel;
    
}szseg;

szseg* sseg_new(_ul beg_i, _ul end_i, szseg* with)
{
    szseg* szr = (szseg*)malloc(sizeof(szseg));
    szr->i0=0;szr->i1=0;szr->i2=0;
    szr->seg=beg_i;szr->off=end_i;
    szr->len=end_i-beg_i;
    szr->corel=(with != 0) ? (with) : (0);
    return(szr);
}

/* sz_t   and   cSz  types .. */


typedef struct _sz_t{
    _c* sz;
    _c* sznam;
    _ul szlen, sznamlen, szindex, szstate;

    struct _sz_t* head, *next;
}sz_t;


typedef struct cSz_t{
   _c* sz;
   _ul szlen;

   _c* ptr[4];
   _ul ptr_i[4];

   _ul  dat_typ;
   _c* dat_num;
   
   _ul obj_indx, obj_state;

   sz_t* szt;

   struct cSz_t* master;
   struct cSz_t* head, *prev, *next;
   struct cSz_t* tangent;   
}cSz;

cSz* csz_alloc_obj(){ return((cSz*)malloc(sizeof(cSz))); }


/*    p_ref  and   o_seg   types.... */


typedef struct p_ref_t{
    _i match_1, match_2;
    _l code;

    _ul beg_seg, beg_off;
    _ul end_seg, end_off;

    _ul seglen_whole, seglen_sub;
    struct p_ref_t* head, *prev, *next;

    _ul index, state, type;
}p_ref;

p_ref* pref_new(){ return( (p_ref*)malloc(sizeof(p_ref)) ); }

typedef struct o_seg_t{
    p_ref* parse;

    _c* runseg;
    _ul runcode;
    _ul runlen;


    struct o_seg_t* head, *prev, *next;

    _ul length, index, state, is_tail, type;
}o_seg;


/*
 *              linenode and wordnode (wordn) types....
 * */


typedef struct onyx_linenode_t{
    _c* line;
    _ul llen;

    struct onyx_linenode_t* head, *prev, *next;
    _ul index, state;
}linenode;

linenode* line_loke(){ return( (linenode*)malloc(sizeof(linenode)) ); }





typedef struct onyx_wordnode_t{
    _c* word;
    _ul wordlen;

    struct onyx_wordnode_t *head, *prev, *next;
    _ul index, state, type;
}wordnode;
#define wordn   wordnode

wordn* wordloke(){ re( (wordn*)malloc((sizeof(wordn))) ); }



/*  wnode stuff ....             */


typedef struct onyx_wnode_t {
    _c* w;
    _ul wlen;
    _ul wtype;
    
    struct onyx_wnode_t* parent, *head, *prev, *next, *tangent;
    _ul index, archetype, subtype, state;
}wnode;



_v wnode_ptrs(wnode* w, wnode* par, wnode* hed, wnode* pre, wnode* nex, wnode* tan)
{ w->parent=par;w->head=hed;w->next=nex;w->prev=pre;w->tangent=tan; }

_v wnode_attribs(wnode* w, _ul indx, _ul atyp, _ul subt, _ul stat)
{ w->index=indx;w->archetype=atyp;w->subtype=subt;w->state=stat; }


/*   lnode  and dnode  stuff */

typedef struct onyx_lnode_t {
    _c* l;
    _ul llen;
    _ul ltype;
    
    wnode* wObj, *wRef;
    _ul wCount;
    _ul wRefSymmetry; //which line contains our ROW?COL model, if any. else 0
    
    struct onyx_lnode_t* parent, *head, *prev, *next, *tangent;
    _ul index, archetype, subtype, state;
}lnode;

typedef struct onyx_dnode_t{
    _c* sz;
    _ul n_ul;
    mpz_t* n_ap;

    lnode* line;
    wnode* word, *wutl;
    _ul s_len, n_ap_len, w_len;
    _ul w_cnt;

    struct onyx_dnode_t* root, *head, *prev, *next, *stem;
    _ul is_root, is_stem;
    _ul index, archetype, subtype, state;
}dnode;

_v lnode_ptrs(lnode* w, lnode* par, lnode* hed, lnode* pre, lnode* nex, lnode* tan)
{ w->parent=par;w->head=hed;w->next=nex;w->prev=pre;w->tangent=tan; }

_v lnode_attribs(lnode* w, _ul indx, _ul atyp, _ul subt, _ul stat)
{ w->index=indx;w->archetype=atyp;w->subtype=subt;w->state=stat; }


/*               onyx executable script  types and objects .....
 *
 *
 * */

typedef struct onyx_funcline_t {
    _c* sz;
    _ul len;
    struct onyx_funcline_t* head, *prev, *next;
    _ul index, archetype, subtype, state;
}funcline;

typedef struct wnode_pair_t {
    wnode* alpha, *beta, *gamma, *delta, *omega;
    _ul szlen[5];

    struct wnodee_pair_t* head, *prev, *next;
    _ul index, archetype, subtype, state;
}wnode_x2;

typedef struct onyx_line_lvalue_t {

    wnode* ws0, *ws1;
    wnode* wp0, *wp1;
    _ul ws0_c, ws1_c;
    _ul wp0_i, wp1_i;

    _ul basic_type, archetype, subtype, state_n, state_t;
}onyx_line_lvalue;

typestruct onyx_logical_object_t onyx_logical_object;
typestruct onyx_logical_construc_t onyx_logical_construct;
typestruct onyx_solution_t onyx_solution;
typestruct onyx_equation_t onyx_equation;
typestruct onyx_formula_t onyx_formula;

typedef struct onyx_operand_t onyx_operand;
typestruct onyx_dFrame_t onyx_dFrame;
#define onyx_dataframe  onyx_dFrame
typestruct onyx_execution_frame_t onyx_execution_frame;
typestruct onyx_func_t onyx_func;
typestruct onyx_operation_t onyx_operation;
typestruct onyx_appstack_t onyx_appstack;



typedef struct onyx_logical_object_t{
    struct onyx_logical_construc_t* construct;

    struct onyx_logical_object_t* x, *y, *z;
    struct onyx_logical_object_t* pos, *neg;

    _ul x_v, y_v, z_v, pos_v, neg_v;

    _ul state, index, type;
}onyx_logical_object;

typedef struct onyx_logical_construc_t{
    onyx_dataframe* elem_in;
    onyx_dataframe* elem_out;
    _ul elem_in_i, elem_out_i;
    _ul elem_in_c, elem_out_c;

    struct onyx_logical_construc_t* x, *y, *z;
    _ul xi, yi, zi;
    _ul state, index, type;
}onyx_logical_construct;

typedef struct onyx_solution_t{
    _c* filename;
    _ul fnam_len, file_siz;

    onyx_logical_object* elements;
    onyx_dFrame* solution_data;

    _ul solution_state, solution_type;
}onyx_solution;

typedef struct onyx_equation_t{


}onyx_equation;

typedef struct onyx_formula_t{
    onyx_logical_object* logic_obj;

    _ul state, type;
}onyx_formula;

typedef struct onyx_operand_t{
    onyx_operand* op0;
    onyx_operand* op1;

    _c* code_sz;
    _ul csz_len;
    _i opcode;

    _ul is_keyword;
    _ul keyword;
    _ul keycode;

    onyx_operand* parent;
    onyx_operand* head, *prev, *next;
    onyx_operand* stem;

    _ul index, count, state, type, is_stem;
    _ul types_in_frame;
}onyx_operand;

typedef struct onyx_dFrame_t{
    _c* name;
    _i nam_len;

    _c* sz_sub0, *sz_sub1;
    _ul sub0_len, sub1_len;

    _c* sz;
    _ul n;
    _ul* n_arr;
    _d flt;

    _ul length;
    _c* len_sz;
    _ul lsz_len;
    _ul nbase;

    mpz_t* n_ap;
    _ul bitwidth_ap, bytewidth_ap;

    struct onyx_dFrame_t* dfram;
    wnode* lval_w, *rval_w;
    _ul archetype, subtype, etype;
    _ul lv_c, rv_c;

    struct onyx_dFrame_t* parent;
    struct onyx_dFrame_t* head, *prev, *next;
    struct onyx_dFrame_t* stem;
    _ul parent_iter, stem_iter;

    _ul framset_mem_tot, framset_cnt;
    _ul thisframe_mem_inuse;
    _ul index, state, type;
    _ul stem_index;
    _i op_type;
}onyx_dFrame;

typedef struct onyx_execution_frame_t{
    onyx_func* parent_fnct;

    _c* opkey_sz;
    _c* codeline_sz;
    _ul linesz_len;

    _ul opline, hasop, opcode;
    onyx_dFrame* line_resolv;
    onyx_dFrame* linepass_alpha;
    onyx_dFrame* dataObjects;
    _ul n_ops, soltn_stat;

    wnode_x2* LRpair;

    onyx_line_lvalue* line_lval;
    onyx_operand* line_ops;
    onyx_operand* this_op;
    dnode* dObj, dUtl;
    wnode* wObj, fig_set, *lval, *rval; //fig_set: onyx refers to an lexical element as a "figure", fuck it. why not?
    lnode* lObj;

    struct onyx_execution_frame_t* parent; //needed so we dont forget from which frame spawned us after forks...
    struct onyx_execution_frame_t* head, *prev, *next;
    struct onyx_execution_frame_t* fork;

    _ul index, state, type;
}onyx_execution_frame;

typedef struct onyx_func_t {
    _c* name;
    _c* param_str;
    _c* body;

    _ul name_len, param_len, body_len;

    _i ret_type;
    _v* ret_stack;

    onyx_dFrame* p_frames;
    onyx_dFrame* p_util;
     _ul framset_mem_tot, framset_count;

    dnode* dataset;
    lnode* lineset;
    wnode* wordset;

    onyx_dataframe* vars_local, *vars_body;
    onyx_dataframe* var_util, *bod_util;
    _ul locvars_mem_tot, locvars_count, bodvars_c;

    onyx_execution_frame* exec_stack;
    onyx_execution_frame* exec_frame;
    _ul stack_frame_tot, stack_frame_i;

    struct onyx_func_t* head, *prev, *next;

    _ul index, state, is_end;
}onyx_func;

onyx_func* onyx_func_tail(onyx_func* f)
{
    onyx_func* util = f;
    for( ; util->next; util = util->next);
    return(util);
}

typedef struct onyx_operation_t
{
    onyx_solution* solution;
    struct onyx_operation_t* head, *prev, *next;
    _ul index, state, type;
}onyx_operation;


typedef struct onyx_appstack_t{

    onyx_dataframe* argz_in;

    onyx_dataframe* global_variables;
    onyx_execution_frame* global_exec_frames;   //these are the "lines" of code
    onyx_func* global_functions;

}onyx_appstack;

/*
     end of type definitions
*/


_ul szlen2(_c* sz)
{
    _ul i0;
    _c* sp = sz;
    for(i0=0; *sp; ++sp, ++i0);;
    return(i0);
}

_ul countch(_c* s0, _c ch);
_c* str2splice(_c* s0, _c* s1, _c delim);
_c* cmdxsplice(_i vN, _c** szV, _c delim);
_c* strxsplice_endi(_i vN, _c** szV, _c delim);
_c* strchain(_i argc, _c** argv);
_c* scatlines(_c* s0, _c excl);
_c* scatlines2(_c* s0, _c* s1);
_c* srepllines(_c* s0, _c excl, _c repl);

_ul ucstrcmp2(_c* s0, _c* s1);
_ul ucstrcmp3(_c* s0, _c* s1);


static const char* apnd_cmd = " > ./cmd.out.put";
static unsigned int apnd_szlen;

_c* locsiz(_ul n) { return( (_c*)malloc(sizeof(_c)*n) ); }
_c* locsiz2(_ul n, _ul pad) { return( (_c*)malloc((sizeof(_c)*n)+pad) ); }

_c* szalloc( _ul len )
{
    _c* sz = (_c*)malloc( (sizeof(_c) * len) +1);
    sz[len]=0;
    return(sz);
}

static _c* case_low = {
'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
                               'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r',
                               's', 't', 'u', 'v', 'w', 'x', 'y', 'z'
};
static _c* case_up = {
                              'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
                               'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
                               'T', 'U', 'V', 'W', 'X', 'Y', 'Z'
};

//if argument string is found anywhere in the **array, return 1;
_i argscan(_c* term, _i argc, _c** argv)
{
    _ul i, mat;

    for(i=mat=0; i < argc; ++i){
    
        if( ucstrcmp3(term, argv[i] )){
            mat = 1;
            break;
        }
    }
    return(mat);
}

_i catchcase(_c c) // <- dry, quaint sort of geek humor. that's me, ..
{                  //    Nathan A. Wallace the "normal guy". :)
    _ul i0, i1;
    _i chr_case = 0;
    _ul szup_len = szlen2(case_low);
    _ul szlow_len = szlen2(case_up);
    _ul match_up, match_low;
    for( i0 = match_up = match_low = i1 = 0; i0 < szlow_len; ){
        if(c == case_low[ i0 ]){
            match_low = i0;
            i0++;
            continue;
        }else{
            match_low = 0;
            break;
        }
    }
    if(!match_low){
        for( ; i1 < szup_len; ){
            if( c == case_up[ i1 ]){
                match_up = i1;
                i1++;
                continue;
            }else{
                match_up = 0;
                break;
            }
        }
        if(!match_up){
            chr_case = 0; //NO MATCH in any case :(
        }else{
            chr_case = 2; // 2 == is_uppercase
        }
    }else{
        chr_case = 1; // 1 == is_lowercase
    }
    
    return(chr_case);
}

_i chrcmp_case(_c c0, _c c1)
{
    

}

//return first occuring position of chr in vect, if any. otherwise -1
_l chpos(_c* vect, _c chr )
{
    _ul i0, i1;
    _ul match0, match1;
    _ul vcnt = szlen2(vect);
    for( i0 = 0; i0 < vcnt; ){
        if( vect[i0] == chr )
        {
            match0 = 1;
            break;
        
        }else{
            match0 = 0;
            i0++;
        }
    }
    
    if(match0){
        match1 = i0;
    }else{
        match1 = -1;
    }    
        
    return(match1);
}

_ul szcpy(_c* s0, _c* s1)
{
    _ul i0;
    _c* p0, *p1;
    for(p0 = s1, i0=0; *p0; ++p0, ++i0)
        s0[i0] = *p0;
    return(i0);
}

_ul szcpy_to(_c* s0, _c* s1, _ul to)
{
    _ul i0;
    _c* p0, *p1;
    for(p0 = s1, i0=0; *p0 && *p0 != to; ++p0, ++i0){
        s0[i0] = *p0;
        if(i0 >= to)
            break;
    }
    return(i0);
}
/*
wrote this badass little son of a bitch while trying to figure out a relatively
efficient system for executing the lines of onyx script code.. dont ask how/why ..
*/
static _ul szcpy_exclch(_c* s0, _c* s1, _c excl_ch, _c spec_excp )
{
    _ul i0;
    _c* sp0, *sp1;
    for(i0=0, sp0 = s0, sp1 = s1;  *sp1; ++sp0, ++sp1){
        if( *sp1 == spec_excp ){
            s0[ i0 ] = excl_ch;
            i0++;

        }else{
            if( *sp1 != excl_ch ){
                s0[ i0 ] = *sp1;
                i0++;
            }
        }
    }
    s0[i0++]=0;
    return(i0);
}
/*
_ul szcount_extant_ch(_c* s0 )
{
    _ul i0, i1;
    for( i0=i1=0; s0[i0]; ++i0 )
        if( s0[i0] != ' ' && s0[i0] != '\n' && s0[i0])
            i1++;;
    return(i1);
}
*/
_c* word_isolate1(_c** args)
{
    _ul i0, i1, i2, i3;
    _c* sp0, *sp1, *sp2, *sp3;

    _c* word = 0;
    _ul wordlen = 0;

    _c* szin = fget_data(args[1]);
    _ul word_index = baseX_sz( args[2], 10 );

    if(szin[0] == ' '){
        for( sp0 = szin; *sp0 == ' '; ++sp0);
        
    } else {
        sp0 = szin;    
    }

    for( i0 = 0, sp1 = sp0; *sp1 != ' '; ++sp1, ++i0);
    
    for( i1 = 0; i1 < i0; ++i1){
        printf("%c", *sp0);
        
        sp0++;
    }
    printf("\n");

    return(0);
}

// select a whitespace-delimited word from the string contained inside any file.
//accepts integers (_ul w_n) or strings (_c* wordn_sz) representing integers 
//defining the index of the word you want to grab ...
static _c* word_isolat2(_c* filename, _c* wordn_sz, _ul w_n)
{
    _ul i0, i1, i2, i3;
    _c* sp0, *sp1, *sp2, *sp3;

    _c* word = 0;
    _ul wordlen = 0;

    _c* szin = fget_data(filename);
    _ul word_index = (!w_n) ? baseX_sz( wordn_sz, 10 ) : w_n;

    sp3 = szin;

    for( i3 =0; i3 < word_index; ++i3){

        if( *sp3 == ' ' || *sp3 == '\n'){
            for( sp0 = sp3; *sp0 == ' ' || *sp0 == '\n'; ++sp0);
        
        } else {
            sp0 = szin;    
        }

        for( i0 = 0, sp1 = sp0; *sp1 != ' ' && *sp1 != '\n'; ++sp1, ++i0);
    
        if(i3 == word_index-1){
            wordlen = i0;
            word = szalloc( wordlen );
            
            for( i1 = 0; i1 < i0; ){
            
                if( *sp0 != ' ' && *sp0 != '\n'){
                    word[ i1 ] = *sp0;
                    i1++;
                }
                sp0++;
            }
            word[i1] = 0;        
        }else{        
            for( i1 = 0; i1 < i0; ++i1)
                ++sp0;        
        }
        sp3 = sp0;    
    }
    return(word);
}


/*
    no different than above except that word_isolat3(..) does not open the file for you and
    assume indices into otherwise arbitrarily referenced datasets which, in any event, may
    obvioulsy already have been parsed into whatever arrays, structures or node objects and
    hence going back to the beginning unstructured format and picking out single words from
    that stage in processing seems not by far an impressive stride... 
*/
static _c* word_isolat3(_c* szin, _c* wordn_sz, _ul w_n )
{
    _ul i0, i1, i2, i3;
    _c* sp0, *sp1, *sp2, *sp3;

    _c* word = 0;
    _ul wordlen = 0;

    _ul word_index = (!w_n) ? baseX_sz( wordn_sz, 10 ) : w_n;

    sp3 = szin;

    for( i3 =0; i3 < word_index; ++i3){

        if( *sp3 == ' ' || *sp3 == '\n'){
            for( sp0 = sp3; *sp0 == ' ' || *sp0 == '\n'; ++sp0);
        
        } else {
            sp0 = szin;    
        }

        for( i0 = 0, sp1 = sp0; *sp1 != ' ' && *sp1 != '\n'; ++sp1, ++i0);
    
        if(i3 == word_index-1){
            wordlen = i0;
            word = szalloc( wordlen );
            
            for( i1 = 0; i1 < i0; ){
            
                if( *sp0 != ' ' && *sp0 != '\n'){
                    word[ i1 ] = *sp0;
                    i1++;
                }
                sp0++;
            }
            word[i1] = 0;        
        }else{        
            for( i1 = 0; i1 < i0; ++i1)
                ++sp0;        
        }
        sp3 = sp0;    
    }
    return(word);
}

/* just like the other word_isolatX(..) functions except this one does comparisons of separation-delimiters against any
   range of -break- characters, so that it can scan to eliminate spaces, newlines, tabs, whatever in the fuck else you
   want....

*/  // Nov 3, 2018, unixcube-openware (c)

_ul chrcmp_vect(_c* dz, _c ch)
{
    _ul i0;
    _ul cmp;
    _c* dp0 = 0;
    for( dp0 = dz, i0 = cmp = 0; *dp0; ++dp0){
        cmp = (ch == *dp0) ? 1 : 0;
        if(cmp)
            break;
        else
            i0++;
    }
    if(cmp)
        cmp = i0;        
    return(cmp);
}
#define chrcmp_any(A,B)   chrcmp_vect(A,B)


static _c* word_isolat4(_c* szin, _c* wordn_sz, _ul w_n, _c* delms )
{
    _ul i0, i1, i2, i3;
    _c* sp0, *sp1, *sp2, *sp3;

    _c* word = 0;
    _ul wordlen = 0;

    _ul word_index = (!w_n) ? baseX_sz( wordn_sz, 10 ) : w_n;
    _ul dlm_c = szlen2(delms);

    sp3 = szin;

    for( i3 =0; i3 < word_index; ++i3){
        
       
        if( chrcmp_vect( delms, *sp3 ) ){
            for( sp0 = sp3; chrcmp_vect( delms, *sp0 ); ++sp0);
        
        } else {
            sp0 = szin;    
        }

        for( i0 = 0, sp1 = sp0; !chrcmp_vect( delms, *sp1 ); ++sp1, ++i0);
    
        if(i3 == word_index-1){
            wordlen = i0;
            word = szalloc( wordlen );
            
            for( i1 = 0; i1 < i0; ){
            
                if( chrcmp_vect( delms, *sp0 )){
                    word[ i1 ] = *sp0;
                    i1++;
                }
                sp0++;
            }
            word[i1] = 0;        
        }else{        
            for( i1 = 0; i1 < i0; ++i1)
                ++sp0;        
        }
        sp3 = sp0;    
    }
    return(word);
}

static _c* word_isolat5(_c* szin, _c* wordn_sz, _ul w_n, _c dlm )
{
    _ul i0, i1, i2, i3;
    _c* sp0, *sp1, *sp2, *sp3;

    _c* word = 0;
    _ul wordlen = 0;

    _ul word_index = (!w_n) ? baseX_sz( wordn_sz, 10 ) : w_n;

    sp3 = szin;

    for( i3 =0; i3 < word_index; ++i3){

        if( *sp3 == ' ' || *sp3 == '\n' || *sp3 == dlm || *sp3 == '\0' ){
            for( sp0 = sp3; *sp0 == ' ' || *sp0 == '\n'  || *sp3 == dlm || *sp3 == '\0'; ++sp0);

        } else {
            sp0 = szin;
        }

        for( i0 = 0, sp1 = sp0; *sp1 != ' ' && *sp1 != '\n'  && *sp3 == dlm && *sp3 == '\0'; ++sp1, ++i0);

        if(i3 == word_index-1){
            wordlen = i0;
            word = szalloc( wordlen );

            for( i1 = 0; i1 < i0; ){

                if( *sp0 != ' ' && *sp0 != '\n'  && *sp3 == dlm && *sp3 == '\0'){
                    word[ i1 ] = *sp0;
                    i1++;
                }
                sp0++;
            }
            word[i1] = 0;
        }else{
            for( i1 = 0; i1 < i0; ++i1)
                ++sp0;
        }
        sp3 = sp0;
    }
    return(word);
}
_ul onyx_wordtag_isolat_and_ref(_c* s, _c* wordno_sz, _ul wordno_num, _c** tagset, _ul n_tags)
{
 
    _ul i0, i1, i2, i3;
    _ul wordtag_index = 0;
    _c* word_sz = word_isolat3(s, wordno_sz, wordno_num);

    for(i0 = 0; i0 < n_tags; ++i0){
        wordtag_index = ucstrcmp3( word_sz, tagset[i0] );
        if(!wordtag_index){
            wordtag_index = i0;
            break; 
        } else {
            wordtag_index = 0;
            continue;
        }
    }
    return(wordtag_index);
}

#define onyx_word(a,b,c)       word_isolat3(a,b,c)


_c* line_isolate1(_c** args)
{
    _ul i0, i1, i2, i3;
    _c* sp0, *sp1, *sp2, *sp3;

    _c* word = 0;
    _ul wordlen = 0;

    _c* szin = fget_data(args[1]);
    _ul word_index = baseX_sz( args[2], 10 );

    if(szin[0] == '\n'){
        for( sp0 = szin; *sp0 == '\n'; ++sp0);
        
    } else {
        sp0 = szin;    
    }

    for( i0 = 0, sp1 = sp0; *sp1 != '\n'; ++sp1, ++i0);
    
    for( i1 = 0; i1 < i0; ++i1){
        printf("%c", *sp0);
        
        sp0++;
    }
    printf("\n");

    return(0);
}

static _c* line_isolat2(_c* filename, _c* wordn_sz, _ul w_n)
{
    _ul i0, i1, i2, i3;
    _c* sp0, *sp1, *sp2, *sp3;

    _c* word = 0;
    _ul wordlen = 0;

    _c* szin = fget_data(filename);
    _ul word_index = (!w_n) ? baseX_sz( wordn_sz, 10 ) : w_n;

    sp3 = szin;

    for( i3 =0; i3 < word_index; ++i3){

        if(  *sp3 == '\n'){
            for( sp0 = sp3; *sp0 == '\n'; ++sp0);
        
        } else {
            sp0 = szin;    
        }

        for( i0 = 0, sp1 = sp0; *sp1 != '\n'; ++sp1, ++i0);
    
        if(i3 == word_index-1){
            wordlen = i0;
            word = szalloc( wordlen );
            
            for( i1 = 0; i1 < i0; ){
            
                if( *sp0 != '\n'){
                    word[ i1 ] = *sp0;
                    i1++;
                }
                sp0++;
            }
            word[i1] = '\n';        
        }else{        
            for( i1 = 0; i1 < i0; ++i1)
                ++sp0;        
        }
        sp3 = sp0;    
    }
    return(word);
}

static _c* line_isolat3(_c* szin, _c* wordn_sz, _ul w_n, _c delm )
{
    _ul i0, i1, i2, i3;
    _c* sp0, *sp1, *sp2, *sp3;

    _c* word = 0;
    _ul wordlen = 0;

    _ul word_index = (!w_n) ? baseX_sz( wordn_sz, 10 ) : w_n;

    sp3 = szin;

    for( i3 =0; i3 < word_index; ++i3){

        if(  *sp3 == delm){
            for( sp0 = sp3; *sp0 == delm; ++sp0);
        
        } else {
            sp0 = szin;    
        }

        for( i0 = 0, sp1 = sp0; *sp1 != delm; ++sp1, ++i0);
    
        if(i3 == word_index-1){
            wordlen = i0;
            word = szalloc( wordlen );
            
            for( i1 = 0; i1 < i0; ){
            
                if( *sp0 != delm){
                    word[ i1 ] = *sp0;
                    i1++;
                }
                sp0++;
            }
            word[i1] = delm;        
        }else{        
            for( i1 = 0; i1 < i0; ++i1)
                ++sp0;        
        }
        sp3 = sp0;    
    }
    return(word);
}

_ul szlen(_c* s)
{
    _ul i0;
    _c* p = s;
    for(i0=0; *p; ++p) i0++;
    return(i0);
}

_ul szlen_tonex_ch(_c* sz, _c ch)
{
    _ul i0, i1;
    _c* sp0 = sz;
    for(i0=0; *sp0 != ch && *sp0; ++sp0, ++i0);
    return(i0);
}

static _ul skip_whitespaces(_c* sp)
{
    _ul i0 = 0;
    do{
        if(*sp == ' ' || *sp == '\n')
            sp++;
        else
            break;

        i0++;

    }while(*sp);

    return(i0);
}
_v sz0(_c* s, _ul c)
{
    _ul i0=0;
    for(; i0 < c; ++i0) s[i0]=0;
}

_v sz_skip_whtspc(_c* p)
{
    for(; *p == ' ' || *p == '\n'; ++p);
}

_v sz_skip_delim(_c* p, _c d0)
{
    for(; *p == d0; ++p);
}

_v sz_skip_delim2(_c* p, _c d0, _c d1)
{
    for(; *p == d0 || *p == d1; ++p);
}

_ul ucstrcmp3(_c* s0, _c* s1)
{
    _c* p0, *p1;
    _ul match=0;
    for(p0=s0, p1=s1;; ){
        if(*p0 == *p1){
            match++;
            if(*(p0+1) == 0 || *(p1+1) == 0 ){
                break;
            }else{
                p0++;
                p1++;
            }
        }else{
            match = 0;
            break;
        }
    }
    return(match);
}

_ul onyx_is_at_operator(_c* codeseg, _c** tagset)
{
    _ul i0, i1, i2;
    _c* sp0, *sp1;

    _ul opcode;

    for(i0=0; !ucstrcmp2(tagset[i0], "+++" ); ++i0){
        if( ucstrcmp3( tagset[i0], codeseg) ){

            opcode = i0;
            break;
        }else{
            opcode=0;
        }
    }
   // opcode++;// += 1;
    return(opcode);
}

_ul onyx_has_operator(_c* codeseg, _c** tagset)
{
    _c* sp0, *sp1;

    _ul opcode = 0;

    for(sp0 = codeseg; *sp0; ++sp0){
        opcode = onyx_is_at_operator(sp0, tagset);
        if(opcode)
            break;
    }
    return(opcode);
}

static _ul ucstrcmp_seg(_c* s0, _c* s1, _ul s0_a, _ul s0_b, _ul s1_a, _ul s1_b)
{
    
    
    return(0);
}

static _ul ucstrcmp_to( _c* s0, _c* s1, _ul s0_a, _ul s0_b, _ul s1_a, _ul s1_b )
{
   return(0);
}

static _ul ucstrcmp_against(_c* sz, _ul beg, ...)
{
    return(0);
}

_ul ucstrlen(_c* s){
 
    _ul i0;
    _c* p = s;
    for(i0=0; *p; ++p)i0++;
    return(i0);
}

_ul ucstrcpy(_c* s0, _c* s1)
{
    _ul i0;
    _c* p0 = ( s0 != 0 ) ? 
        s0 : 
        ( (_c*)malloc( (sizeof(_c) * ucstrlen(s1)+1)) );
        
    _c* p1 = s1;
    _c* pX = p0;
    for(i0=0; *p1; ++p1, ++pX, ++i0)
        *pX = *p1;
    
    pX[ i0 ] = 0;
    
    return(i0);       
}

//returns 0 on s0==s1, else returns first array index at which variance occurs
_ul ucstrcmp2(_c* s0, _c* s1)
{
    _ul i0, i1;
    _ul slen0 = ucstrlen(s0);
    _ul slen1 = ucstrlen(s1);
   
    if(slen0 != slen1)
        return(0);
    
    for(i1=0, i0=0; i0 <= slen0-1; ++i0){
        if(s0[i0] == s1[i0]){
            i1++;
        }else{ 
            i1=0;
            break;
        }
    }
    
    return(i1);
}

static _c* szchomp(_c* s0, _ul seg_a, _ul seg_b )
{
    _ul i0, i1;
    _ul slen = seg_b-seg_a;
    _c* sR  =(_c*)malloc((sizeof(_c)*slen)+1);
    _c* sP0 = s0;
    
    for(i0=0, i1=seg_a; i1 < seg_b; ++i1, ++i0)
        { sR[i0] = s0[i1]; }
    sR[ slen ] = '\0';

    return(sR);
}
/*
 *  cSz functions::
 * 
 * 
 * 
 * 
 * 
 * 
 */



_ul countch(_c* s0, _c ch)
{
  _ul lin;
  _c* sP = s0;
  for(lin=0, sP=s0; *sP; ++sP){
    if(*sP == ch)
      lin++;
  }
  return(lin);
}

_c* str2splice(_c* s0, _c* s1, _c delim)
{
    _ul i0, i1;
    _ul dmlen = (delim)?1:0;
    _ul s0len = strlen((_cc*)s0);
    _ul s1len = strlen((_cc*)s1);
    _ul sXlen = ((s0len + s1len )+ dmlen) + 1;

    _c* sz = (_c*)malloc((sizeof(_c*) * sXlen) );

    for(i0=0; i0 <= s0len-1; ++i0){ sz[i0] = s0[i0]; }

    if(delim)
        sz[i0++] = delim;

    for(i1=0; i1 <= s1len-1; ++i0, ++i1){ sz[i0] = s1[i1]; }
    sz[i0++] = '\0';

    return(sz);
}

_c* cmdxsplice(_i vN, _c** szV, _c delim)
{
    _ul i0, i1;
    _c* szX = 0;
    //_c* szT = str2splice(szV[0], szV[1], delim);
    _c* szT = str2splice(szV[1], szV[2], delim);
    for(i0=3; i0 <= vN-1; ++i0){
       // szX = str2splice(szV[i0], szT, delim);
        szX = str2splice(szT, szV[i0], delim);
        free(szT);
        szT = szX;
    }
    return(szX);
}

_c* strxsplice_endi(_i vN, _c** szV, _c delim)
{
    _ul i0, i1;
    _c* szX = 0;
    //_c* szT = str2splice(szV[0], szV[1], delim);
    _c* szT = str2splice(szV[1], szV[0], delim);
    for(i0=2; i0 <= vN-1; ++i0){
        szX = str2splice(szV[i0], szT, delim);
       // szX = str2splice(szT, szV[i0], delim);
        free(szT);
        szT = szX;
    }
    return(szX);
}


_c* strchain(_i argc, _c** argv)
{
    _i szc = argc-1;
    _c delm = argv[argc-1][0] != '0'? argv[argc-1][0]: 0;
    if(delm == 'S')
        delm = ' ';
    if(delm == '0')
        delm = 0;
    if(delm == 'N')
        delm = '\n';
    
    _c* r = cmdxsplice(szc, argv, delm);
    //_c* r_lendi = strxsplice_endi(argc, argv, delm);
    //printf("%s\n", r);
    //printf("%s\n\n", r_lendi);
    return(r);
    //free(r_lendi);
}

_c* scatlines(_c* s0, _c excl)
{
  _ul i0;
  _ul n_excl = countch(s0, excl);
  _ul szlen = (_ul)( ( sizeof(_c)*(strlen((_cc*)s0) ) - n_excl) + 1);
  _c* szr = (_c*)malloc( sizeof(_c) * szlen );
  _c* sX = szr;
  _c* sP = s0;
  
  memset(szr, 0, szlen);
  
  for(i0=0; *sP; ++sP){
    if( *sP != excl ){
      *sX = *sP;
      sX++;
    }
  }
  return(szr);  
}

_c* scatlines2(_c* s0, _c* s1)
{
  _c* sP = 0;
  _ul i0, i1, i2;
  _ul s0len = strlen((_cc*) s0);
  _ul s1len = strlen( (_cc*) s1);
  _ul veclen = ( s0len + s1len ) - 
    (countch( s0, '\n' ) + countch( s1, '\n' ) );
  
  _c* r = (_c*)malloc( ( sizeof(_c) * ( veclen ) ) + 1 );
  
  for(i0 = 0, sP=s0; *sP; ++sP)
    { if(*sP != '\n') { r[ i0 ] = *sP; i0++; } }
  
  for(i1 = i0, sP=s1; *sP; ++sP)
    { if(*sP != '\n') { r[ i1 ] = *sP; i1++; } }
    
  r[i1++] = '\0';
  
  return(r);
}

_c* srepllines(_c* s0, _c excl, _c repl)
{
  _ul i0;
  _ul n_excl = countch(s0, excl);
  _ul szlen = (_ul)( ( sizeof(_c)*(strlen((_cc*)s0) ) - n_excl) + 1);
  _c* szr = (_c*)malloc( sizeof(_c) * szlen );
  _c* sX = szr;
  _c* sP = s0;
  
  memset(szr, 0, szlen);
  
  for(i0=0; *sP; ++sP){
    if( *sP != excl )
      *sX = *sP;
    else
      *sX = repl;
    
    sX++;
  }
  return(szr);
}


_ul getbufsize_allfiles(_c** flist, _ul fcount, _ul fstart, _c excl)
block_
  _ul i0;
  _c* sz0;
  _ul r=0;
  for(i0=fstart; i0 <= fcount; ++i0) block_
    sz0= fget_data(flist _a(i0) );
    r += ( ( sizeof(_c)*(strlen((_cc*)sz0) ) - countch(sz0, excl)) + 1);
    free(sz0);
  _block
  rt(r);
_block

_ul scat(_i argc, _c** argv)
block_

  _ul i0, i1, i2;
  _c* sz0, szX;

  for(i0=1; i0 <= argc-1; ++i0){
    sz0= fget_data(argv[i0]);
  
    printf("%s", scatlines(sz0, '\n') );
       
    free(sz0);
  }
  printf("\n");
  return(0);
_block


/*
            sz_t   functions ...            
*/

sz_t* sz_loke(){ return((sz_t*)malloc(sizeof(sz_t)) ); }

sz_t* sz_getend(sz_t* szt)
{
    sz_t* util=0;
    for(util=szt; util->next; util=util->next);;

    return(util);
}

sz_t* sz_allocate(_ul szlen)
{
    _ul len = !szlen ? 0 : szlen;
    sz_t* n = sz_loke();
    if(len)
        n->sz = (_uc*)malloc((sizeof(_c) * len)+1);
    else
        n->sz=0;

    n->sznam=0; n->szlen = len; n->sznamlen=0; n->szindex=0; 
    n->szstate= (!len ? 1 : 2);
    n->head=n; n->next=0;

    return(n);
}

sz_t* sz_spawn(_c* str)
{
    _ul newlen = (_ul)strlen((_cc*)str);
    sz_t* newsz = sz_allocate( newlen );
    strcpy(newsz->sz, str);

    newsz->szstate = 2;
    return(newsz);
}

sz_t* sz_append(sz_t* head, _c* str)
{
    sz_t* end = sz_getend(head);
    sz_t* spawn = sz_spawn(str);
    end->next = spawn;
    return(spawn);
}

/*
         cSz functions section...

*/

_v csz_assign_nodeptrs(cSz* objOp, cSz* master, cSz* head, cSz* prev, cSz* next, cSz* tangent)
{ objOp->master=master; objOp->head=head; objOp->prev=prev;objOp->next=next; objOp->tangent=tangent; }

cSz* csz_alloc_bysiz(_ul slen)
{
    cSz* szObj = csz_alloc_obj();
    szObj->szlen = slen;
    szObj->sz = (_c*)malloc((sizeof(_c*) * szObj->szlen)+1 );

    szObj->ptr[0]=0; szObj->ptr[1]=0;szObj->ptr[2]=0; szObj->ptr[3]=0;
    szObj->ptr_i[0]=0; szObj->ptr_i[1]=0;szObj->ptr_i[2]=0; szObj->ptr_i[3]=0;

    szObj->dat_typ=0; szObj->dat_num=0; szObj->obj_indx=0; szObj->obj_state=1;

    csz_assign_nodeptrs(szObj, 0, 0, 0, 0, 0);

    return(szObj);
}

cSz* csz_alloc_str(_c* str, _ul data_type)
{
    cSz* szObj = csz_alloc_bysiz( (_ul)strlen((_cc*)str) );
    strcpy(szObj->sz, str);

    szObj->dat_typ = data_type;
    if(szObj->dat_typ){
        szObj->dat_num = 0; //cnm_alloc_obj();

        //since type signifies integer here,
        //lets convert the string to its intended
        //format now ...
    }

    return(szObj);
}

cSz* csz_get_end(cSz* node)
{
    cSz* util = 0;
    for(util=node; util->next; util = util->next);;
    return(util);
}

_ul csz_count_nodes(cSz* szn)
{
    _ul n0;
    cSz* util=csz_get_end(szn);
    return(util->obj_indx);    
}

cSz* csz_make_fromseg(_c* sz, _ul seg, _ul off)
{
    _ul i0, i1;
    
    cSz* rsz = csz_alloc_bysiz( off-seg );

    for(i0=0, i1=seg; i1 <= off-1; ++i0, ++i1)
        rsz->sz[i0] = sz[i1];

    return(rsz);
}

_ul csz_count_occrch(_c* sz, _ul seg, _ul off, _c ch)
{
    _ul i0, i1;
    for(i0=seg, i1=0; i0 <= off-1; ++i0){
        if( ch == 'N' || ch == 'Z' || ch == 'S'){
            if(ch == 'N' && sz[i0] == '\n'){
                i1++;
            }
            if(ch == 'Z' && sz[i0] == '\0'){
                i1++;
            }

            if(ch == 'S' && sz[i0] == ' '){
                i1++;
            }
        }else{
            if(sz[i0] == ch) i1++;
        }
    }
    return(i1);
}

_ul csz_count_occrchset(_c* sz, _ul seg, _ul off, _c* chset)
{
    _ul i0;
    _ul excl_tot=0;
    _ul setlen = strlen((_cc*)chset);
   
    for(i0=0; i0 <= setlen-1; ++i0){
        excl_tot += csz_count_occrch(sz, seg, off, chset[i0]);
    }
    return(excl_tot);
}
/*
_ul csz_countch_occr(_sz, _ul seg, _ul off, _c* chset)
{
    _ul i0, i1;
    _ul occr_tot=0;
    
}*/

//exclude char "chr" (ie. whitespaces) from copy
cSz* csz_make_fromseg_excl_ch(_c* sz, _ul seg, _ul off, _c chr)
{
    _ul i0, i1;
    _ul n_excl_ch = csz_count_occrch(sz, seg, off, chr);
    cSz* rsz = csz_alloc_bysiz( (off-seg) - n_excl_ch );

    for(i0=0, i1=seg; i1 <= off-1; ++i0, ++i1){
        if(sz[i0] != chr)
            rsz->sz[i0] = sz[i1];
    }
    return(rsz);
}

cSz* csz_append_str_b(cSz* obj, _c* sz)
{
    cSz* szObj = csz_alloc_str(sz, 0);
    cSz* szEnd = csz_get_end(obj);

    szObj->obj_indx = szEnd->obj_indx + 1;

    szEnd->next = szObj;
    return(szObj);
}

_ul csz_chopout_seg(cSz* head, _c* sz, _c* delim)
{
    _ul match=0;
    _ul i0, i1, i2;

    _ul dlm_len = strlen((_cc*)delim);
    _ul wh_len = strlen((_cc*)sz);

    cSz* nod = csz_get_end(head);
    _ul new_seg = nod->ptr_i[2];
    
    for(i0=nod->ptr_i[2]; i0 <= wh_len-1; ++i0){
        if(sz[i0] == delim[0]){
            for(i1=0, i2=i0; i1 <= dlm_len-1 && i2 <= wh_len-1; ++i1, ++i2){
                if(sz[i2] == delim[i1]){
                    match++;
                }else{
                    match=0;
                    break;
        }   }   }
        if(match){                  
            match=i2;

            nod->next = csz_alloc_bysiz(i0 - new_seg);
            nod->next->ptr_i[0] = new_seg;
            nod->next->ptr_i[1] = i0;
            nod->next->ptr_i[2] = match;
            for(i0=0, i1=new_seg; i0 <= nod->next->szlen-1; ++i0, ++i1)
                nod->next->sz[i0] = sz[i1];;

            nod->next->obj_indx = nod->obj_indx+1;
            nod->next->obj_state=2;

            nod->next->head - nod->head;
            nod->next->prev = nod;

            break;
        }else{
            if(i0 >= wh_len-1)
                break;
        }
    }
   
    return(match);
}

_ul csz_chopout_seg_exclch(cSz* head, _c* sz, _c* delim, _c exclch)
{
    _ul match=0;
    _ul i0, i1, i2;

    _ul excl_len = 0;
    _ul dlm_len = strlen((_cc*)delim);
    _ul wh_len = strlen((_cc*)sz);

    cSz* nod = csz_get_end(head);
    _ul new_seg = nod->ptr_i[2];
    
    for(i0=nod->ptr_i[2]; i0 <= wh_len-1; ++i0){
        if(sz[i0] == delim[0]){
            for(i1=0, i2=i0; i1 <= dlm_len-1 && i2 <= wh_len-1; ++i1, ++i2){
                if(sz[i2] == delim[i1]){
                    match++;
                }else{
                    match=0;
                    break;
        }   }   }
        if(match){                  
            match=i2;
          
            excl_len = csz_count_occrch(sz, new_seg, i0, exclch);
            nod->next = csz_alloc_bysiz(( i0 - new_seg ) - excl_len);

            nod->next->ptr_i[0] = new_seg;
            nod->next->ptr_i[1] = i0;
            nod->next->ptr_i[2] = match;

            for(i0=0, i1=new_seg; i0 <= nod->next->szlen-1; ++i1){
                if( sz[i1] != exclch ){
                    nod->next->sz[i0] = sz[i1];
                    i0++;
                }
            }

            nod->next->obj_indx = nod->obj_indx+1;
            nod->next->obj_state=2;

            nod->next->head - nod->head;
            nod->next->prev = nod;

            break;
        }else{
            if(i0 >= wh_len-1)
                break;
        }
    }
   
    return(match);
}



_ul csz_chopout_seg_exclchset(cSz* head, _c* sz, cSz* delim_set, _c exclch)
{
    _ul ii;
    _ul match=0;
    _ul delim_type=0;
    _ul i0, i1, i2, i3;

    _ul excl_stat=0;
   // _ul spec_excl=0;
    _ul excl_len = 0;
    _c* delim=0;    
    _ul dlm_len = 0; //strlen((_cc*)delim);
   // _ul excl_setlen = strlen((_cc*)exclch);
    
    _ul wh_len = strlen((_cc*)sz);

    cSz* tmp_nod=0;
    cSz* nod = csz_get_end(head);
    _ul new_seg = nod->ptr_i[2];
    
    for(i0=nod->ptr_i[2]; i0 <= wh_len-1; ++i0){
        for(tmp_nod=delim_set, ii=0 ; tmp_nod; ++ii, tmp_nod=tmp_nod->next){
            delim = tmp_nod->sz ? : 0;
            if(delim){
                dlm_len = strlen((_cc*)delim);
                if(sz[i0] == delim[0]){
                    for(i1=0, i2=i0; i1 <= dlm_len-1 && i2 <= wh_len-1; ++i1, ++i2){
                        if(sz[i2] == delim[i1]){
                            match++;
                        }else{
                            match=0;
                            break;
                }   }   }
                if(match){
                    delim_type = ii;
                    break;
                }
            }else{

            }
            if(match){
                delim_type = ii;
                break;
            }
        }
        
        if(match){                  
            match=i2;
          
            excl_len = csz_count_occrch(sz, new_seg, i0, exclch);
            nod->next = csz_alloc_bysiz(( i0 - new_seg ) - excl_len);

            nod->next->dat_typ = delim_type;
            nod->next->ptr_i[0] = new_seg;
            nod->next->ptr_i[1] = i0;
            nod->next->ptr_i[2] = match;

            for(i0=0, i1=new_seg; i0 <= nod->next->szlen-1; ++i1){
                excl_stat=0;
                if(exclch == 'N' || exclch == 'Z' || exclch == 'S'){
                    if(sz[i1] == '\n' || sz[i1] == '\0')
                        excl_stat=1;

                    if(exclch == 'S' && sz[1] == ' ')
                        excl_stat=1;

                }else{
                    if(sz[i1] == exclch)
                        excl_stat=1;
                }
                
                if( !excl_stat ){
                    nod->next->sz[i0] = sz[i1];
                    i0++;
                }
            }
            nod->next->sz[i0+1]='\0';

            nod->next->obj_indx = nod->obj_indx+1;
            nod->next->obj_state=2;

            nod->next->head - nod->head;
            nod->next->prev = nod;

            break;
        }else{
            if(i0 >= wh_len-1)
                break;
        }
    }
   
    return(match);
}

cSz* csz_break_sz(_c* sz, _c* delim_sz, _c exclch)
{
    _ul n=0;
    cSz* szObj = csz_alloc_obj();
    szObj->ptr_i[0]=0;
    szObj->ptr_i[1]=0;
    szObj->ptr_i[2]=0;
    szObj->ptr_i[3]=0;

    szObj->master=0; szObj->tangent=0;
    szObj->head = szObj;
    szObj->next = 0;

    for(; csz_chopout_seg_exclch(szObj, sz, delim_sz, exclch) > 0;);;
    return(szObj);
}
cSz* csz_break_sz_exclset2(_c* fnam, _c* delim_db_fnam, _c* delim_sz, _c exclsz)
{
    _ul n=0;
    _c* sz = fget_data(fnam);
    _c* dlm_sz = fget_data(delim_db_fnam);
    cSz* szDlm =  csz_break_sz(dlm_sz, delim_sz, exclsz);
    cSz* szObj = csz_alloc_obj();

    szObj->ptr_i[0]=0;
    szObj->ptr_i[1]=0;
    szObj->ptr_i[2]=0;
    szObj->ptr_i[3]=0;

    szObj->master=0; szObj->tangent=0;
    szObj->head = szObj;
    szObj->next = 0;
    
    for(; csz_chopout_seg_exclchset(szObj, sz, szDlm, exclsz) > 0;);;
    return(szObj);
}



//way better functioning parser,... works better with a clearly defined
//tokenset that does not contain overlapping or duplicate token segments

_ul csz_chopout_wtoks(cSz* head, _c* sz, cSz* delim_set, mpz_t* flags)
{
    _ul ii;
    _ul match=0;
    _ul delim_type=0;
    _ul i0, i1, i2, i3;
    _ul leav_indx;

    _ul excl_stat=0;
   // _ul spec_excl=0;
    _ul excl_len = 0;
    _c* delim=0;    
    _ul dlm_len = 0; //strlen((_cc*)delim);
   // _ul excl_setlen = strlen((_cc*)exclch);
    
    _ul wh_len = strlen((_cc*)sz);

    cSz* tmp_nod=0;
    cSz* nod = csz_get_end(head);
    _ul new_seg = nod->ptr_i[2];
    
    for(i0=nod->ptr_i[2]; i0 <= wh_len-1; ++i0){
        for(tmp_nod=delim_set, ii=0 ; tmp_nod; ++ii, tmp_nod=tmp_nod->next){
   //         delim = tmp_nod->sz ? : 0;
	    delim = tmp_nod->sz;
            if(delim){
                dlm_len = strlen((_cc*)delim);
                if(sz[i0] == delim[0]){
                    for(i1=0, i2=i0; i1 <= dlm_len-1 && i2 <= wh_len-1; ++i1, ++i2){
                        if(sz[i2] == delim[i1]){
                            match++;
                        }else{
                            match=0;
                            break;
                }   }   }
                if(match){
                    delim_type = ii;
                    break;
                }
            }else{

            }
            if(match){
                delim_type = ii;
                break;
            }
        }
        
        if(match){                  
            match=i2;
	    leav_indx=i0;
          
       //     excl_len = csz_count_occrch(sz, new_seg, i0, exclch);
            nod->next = csz_alloc_bysiz(( (i0+1)-new_seg ) );

            nod->next->dat_typ = delim_type;
            nod->next->ptr_i[0] = new_seg;
            nod->next->ptr_i[1] = i0;
            nod->next->ptr_i[2] = match;

	    if(flags){
	      
	    
                for(i0=0, i1=new_seg; i0 <= nod->next->szlen-1; ++i1){
	        //for(i0=0, i1=new_seg; i0 <= leav_indx-3; ++i1){
                    excl_stat=0;
		    
		    if(mpz_tstbit(flags, 1) && sz[i1] == '\n'){ excl_stat=1; }
		    if(mpz_tstbit(flags, 2) && sz[i1] == '\0'){ excl_stat=1; }
		    if(mpz_tstbit(flags, 3) && sz[i1] == ' '){ excl_stat=1; }
		    if(mpz_tstbit(flags, 4) && sz[i1] == ' '){ excl_stat=1; }
		    
                    if( !excl_stat ){
                        nod->next->sz[i0] = sz[i1];
                        i0++;
                    }
                    
		    excl_stat=0;     
                }
                
	    }else{
	     
	        for(i0=0, i1=new_seg; i0 <= nod->next->szlen-1; ++i1){
                
                    nod->next->sz[i0] = sz[i1];
                    i0++;
                
                }
	    }
            
            nod->next->sz[i0-1]='\0';

            nod->next->obj_indx = nod->obj_indx+1;
            nod->next->obj_state=2;

            nod->next->head - nod->head;
            nod->next->prev = nod;

            break;
        }else{
            if(i0 >= wh_len-1)
                break;
        }
    }
   
    return(match);
}
cSz* csz_parse_wtoks(_c* fnam, _c* delim_db_fnam, _c* delim_sz)
{
    _ul n=0;
    _c* sz = fget_data(fnam);
    _c* dlm_sz = fget_data(delim_db_fnam);
    cSz* szDlm =  csz_break_sz(dlm_sz, delim_sz, 0);
    cSz* szObj = csz_alloc_obj();
    
    mpz_t* flagz = (mpz_t*)malloc(sizeof(mpz_t) );
    mpz_init2(flagz, 8);
    mpz_set_ui(flagz, 0);
  
    
  //  mpz_setbit(flagz, 3);
    mpz_setbit(flagz, 1);
  //  mpz_setbit(flagz, 2);
    

    szObj->ptr_i[0]=0;
    szObj->ptr_i[1]=0;
    szObj->ptr_i[2]=0;
    szObj->ptr_i[3]=0;

    szObj->master=0; szObj->tangent=0;
    szObj->head = szObj;
    szObj->next = 0;
    
    for(; csz_chopout_wtoks(szObj, sz, szDlm, flagz) > 0;);;
    
    mpz_clear(flagz);
    free(flagz);
    
    printf("freed processed mpz_t*\n");
    return(szObj);
}

// now comes the core onyx system functions ..  ...

_ul count_delims(_c* s0, _c* delm)
{
    _ul i0;
    _c* sP = s0;
    _ul n_delim = 0;
    for(i0=0; *sP; ++sP, ++i0){
        if(*sP == delm[0])
            n_delim++;
    }
    return(n_delim);
}

_i szcmp(_c* s0, _c* s1)
{
    _ul match;
    _c* sp0 = s0;
    _c* sp1 = s1;

    for(match=0; *sp0 && *sp1; ){
        if(*sp0 == *sp1){
            match++;
        }else{
            match=0;
            break;
        }
        sp0++;
        sp1++;
    }
    if(match > 0)
        match = 1;
    else
        match = 0;

    return(match);
}



_c* oszalloc(_ul newlen)
{
    return( (_c*)malloc( (sizeof(_c) * newlen) + 1 ) );
}

_ul strlen_from(_c* sz, _ul seg)
{
    _ul i0;
    _c* p1=sz;

    for(p1 += seg, i0 = 0; *p1; ++i0, ++p1);

    return(i0);
}

_ul strcpy_from(_c* s0, _c* s1, _ul seg0)
{    
    _ul i0;
    _c* p1=s1;   
    
    for(p1 += seg0, i0 = 0; *p1; ++i0, ++p1)
        s0[i0] = *p1;
    
    return(i0);
}

_ul strcpy_nth_delim(_c* s0, _c* s1, _ul nth)
{
    _ul i0;
    _ul delim_slen = szlen2(s1);
    _ul cp_seg, cp_off;
    
    _ul dlm_match;
    _ul dlm_match_count=0;
    
    for(dlm_match_count=0; dlm_match_count < nth; ){
        
        
        
    }    
}

_ul strcpy_to(_c* s0, _c* s1, _c* tok_end, _ul seg0)
{
    _ul i0;
    _ul m0;
    _c* p1=s1;   
    _c* p2=0;
    _i end_match=0;
    
    _ul endtok_slen = szlen2(tok_end);
    
    for(p1 += seg0, i0 = 0; *p1; ++i0, ++p1){
        
        if( *p1 == tok_end[0] ){
            
            p2 = s1;
            p2 += i0;
            
            for(m0=end_match=1; (*p2) && (tok_end[m0]) && (*p2 == tok_end[m0]);++p2, ++m0){
                
                if(*p2 == tok_end[m0] ){
                    end_match++;
                }else{
                    end_match=0;
                    break;
                }                    
                if( (m0 >= endtok_slen) && (end_match <= endtok_slen) ){
                    end_match=0;
                    break;
                }
            }       
        }
        
        if(!end_match)
            s0[i0] = *p1;
        else
            break;
    }
    
   s0[i0++] = '\0';
   
   if(!end_match)
       i0=0;
   
    return(i0);
}

_ul strcmp_seg(_c* s0, _c* s1, _ul seg0)
{
    _ul i0;
   
    _c* p1=s1;
    
    _ul match = 0;
    
    _ul slen = szlen2(s0);
    
    for(i0=seg0 ; *p1 && i0 < slen; ++p1){
        
        if( *p1 == s0[i0] ){
            match++;
        }else{
            match=0;
            break;
        }
               
        i0++;
    }
    
    i0 = (match ? i0 : 0);
    //return( ( match ? 1 : 0 ) );   // <<--return bool basically
    return(i0);                   // <<--return offs basically
}

_c* strchop_aft_token(_c* s0, _c* token, _ul seg0)
{
    _ul i0;
    _ul i1;
    _c* p1=token;
    _c* p0 = s0;
    
    _ul match = 0;
    
    _ul slen = szlen2(s0);
    
    _c* chopsz = 0;
    
    for(i0=seg0 ; *p1 && i0 < slen; ++p1){
        
        if( *p1 == s0[i0] ){
            match++;
        }else{
            match=0;
            break;
        }
               
        i0++;
    }   
    if(match){            
        slen -= match;
    
        chopsz = (_c*)malloc( (sizeof(_c) * slen) + 1);
        
        for(i1=0; s0[i0]; ++i0, ++i1)
            chopsz[i1] = s0[i0];
        
    }else{
        chopsz = 0;        
    }
    
    return(chopsz);    
}

/* strchop_tok(a,b,c,d) is exactly [almost] the same as its predecessor, exempting though that it
 * does not give a shit if it finds a second token [tok_b] succeeding tok_a in the str or not, it
 * will still copy out the remnant of the string notwithstanding.
 */
_c* strchop_tok(_c* s0, _c* tok_a, _c* tok_b, _ul seg0)
{
    _ul i0;
    _ul i1;
    _c* p1=tok_a;
    _c* p0 = s0;
    
    _ul match = 0;
    
    _ul slen = szlen2(s0);
    
    _c* chopsz = 0;
    
    for(i0=seg0 ; *p1 && i0 < slen; ++p1){
        
        if( *p1 == s0[i0] ){
            match++;
        }else{
            match=0;
            break;
        }
               
        i0++;
    }   
    if(match){            
        slen -= match;
    
        chopsz = (_c*)malloc( (sizeof(_c) * slen) + 1);
        
        strcpy_to(chopsz, s0, tok_b, i0);
        
    }else{
        chopsz = 0;        
    }
    
    return(chopsz);  
    
}

_ul schop(_c* chopsz, _c* s0, _c* tok_a, _c* tok_b, _ul seg0)
{
    _ul i0;
    _ul i1;
    _c* p1=tok_a;
    _c* p0 = s0;
    
    _ul match = 0;
    
    _ul slen = szlen2(s0);
    
    //_c* chopsz = 0;
    
    for(i0=seg0 ; *p1 && i0 < slen; ++p1){
        
        if( *p1 == s0[i0] ){
            match++;
        }else{
            match=0;
            break;
        }
               
        i0++;
    }   
    if(match){            
        slen -= match;
    
        chopsz = (_c*)malloc( (sizeof(_c) * slen) + 1);
        
        i1 = strcpy_to(chopsz, s0, tok_b, i0);
        
    }else{
        chopsz = 0;        
    }
    
    return(i1);      
}

_ul dlmch_nth(_c* s, _c dlmch)
{
    _ul i0, i1, i2;
    _ul dlmch_n=0;

    _c* sP = s;
    for( ; *sP; ++sP){
        if(*sP == dlmch)
            dlmch_n++;
    }

    return(dlmch_n);
}

//IMPROVISE LIKE FUCK :: needs work still. do not forget to finish this function, BITCH!!!
_ul dlm_nth(_c* s, _c* dlm)
{
    _ul i0, i1, i2;
    _ul dlm_n=0;

    _ul slen_w = szlen2(s);
    _ul slen_d = szlen2(dlm);

    _ul match_a, match_b;

    for( ; ; ){

    }

    return(dlm_n);
}

/*
    p_ref and o_seg functions ....
*/

o_seg* oseg_tail(o_seg* sg)
{
    o_seg* utl = sg;
    for(; utl->next; utl=utl->next);
    re(utl);
}

 p_ref* parseline(_c* packet_sz, _c** tokszz, _c* end_token, _ul n_toks)
{
    _ul indx;
    _ul iC, iI, iI1;
    _ul mat;
    _l r = -1;
    _ul slen = szlen2(packet_sz);
    _c* szp = packet_sz;
    
    p_ref* ret_parse_ref = pref_new();
    p_ref* pset = ret_parse_ref;

    pset->head = ret_parse_ref;
    pset->prev = 0;
    pset->next = 0;

    pset->code = -1;
    pset->match_1=0; pset->match_2=0;
    pset->index = 0;

    mat = 0;

    indx=0;
    for( iI=0; iI <= slen - 1; ++iI ){

        for( iC=0; iC < n_toks; ++iC ){
            
            if(strcmp_seg(packet_sz, tokszz[iC], iI) ){
                r=iC;
                mat=1;                
            }else{
                r=-1;
                mat=0;
            }

            if(mat){
                pset->code = iC;
                pset->beg_seg = iI;
                pset->beg_off = iI+szlen2( tokszz[iC] );
                
                for(iI1=iI; iI1 < slen; ++iI1){
                    if( strcmp_seg(packet_sz, end_token, iI1) ){
                        pset->match_2 = 1;
                        break;
                    }else{
                        pset->match_2 = 0;
                    }
                }
                if(pset->match_2){
                    pset->end_seg = iI1;
                    pset->end_off = iI1 + szlen2(end_token);

                    pset->seglen_whole = (pset->end_off - pset->beg_seg);
                    pset->seglen_sub = (pset->end_seg - pset->beg_off);
                }else{

                    pset->end_seg = 0;
                    pset->end_off = 0;
                }
                pset->index = indx;
                pset->state = 1;
                pset->type = pset->code;

                pset->next = pref_new();
                pset->next->head = ret_parse_ref;
                pset->next->prev = pset;
                pset->next->next = 0;

                pset = pset->next;
              
                indx++;
            }else{

            }
        }
        if(mat){
            break;
        }else{

        }
    }    
    
    return(ret_parse_ref);
}

o_seg* onyx_preparse(_c* s, _c** dlmset, _c* dlm_end, _ul tcount)
{
    _ul excl;
    _ul cpy;
    _ul i0, i1, i2;
    _ul lensz_w = szlen2(s);
    _c* t0=0;
    _c* t1=0;
    o_seg* omap = (o_seg*)malloc(sizeof(o_seg));
    o_seg* omP = omap;
    omap->parse = parseline(s, dlmset, dlm_end, tcount);
    p_ref* pset = omap->parse;

    for(pset=omap->parse; pset->next; pset=pset->next){
        excl=0;
        i2 = pset->end_seg-2 - pset->beg_off;

        omP->runseg = (_c*)malloc((sizeof(_c)*i2)+1);
        omP->runlen = i2;
        omP->runcode = pset->code;

        i1=0;
        for(i0 = pset->beg_off; i0 <= pset->end_seg-1; ++i0 ){
           

            if(s[i0] != ' ' && s[i0] != '\n' && s[i0] != '\r'){
                excl=0;   
            }else{
                excl++;
            }

            if(excl < 2){
                cpy=1;
            }else{
                cpy=0;
            }
            

            if(cpy){
                omP->runseg[i1] = s[i0];
                i1++;
            }
           
        }
        omP->runseg[i1++] = '\0';

        //printf("%s\n",omP->runseg);
 
        omP->next = (o_seg*)malloc(sizeof(o_seg));
        omP = omP->next;

    }
 
 
    return(omap);

}


//NOTE::: FIX THIS LATER!! this only frees the first node in chain, leaving the rest a memory leak like fuck..
static _v onyx_destroy_pref(p_ref* pr)
{
    free(pr);
}

static _ul onyx_destroy_oseg(o_seg* segnode)
{
    _ul dlen=0;
    o_seg* utl0=0;
    o_seg* utl1=0;
    for(utl0=segnode; utl0->next; ){
        dlen += szlen2(utl0->runseg);
        utl1 = utl0;
        if(utl0->next)
            utl0=utl0->next;
        else
            break;
        onyx_destroy_pref(utl1->parse);
        free(utl1->runseg);
        free(utl1);
    }
    dlen += szlen2(utl0->runseg);
    onyx_destroy_pref(utl0->parse);
    free(utl0->runseg);
    free(utl0);
    re(dlen);
}


static _v onyx_sysexecs_ex(o_seg* segs, _ul syscode)
{
    _c* buf_sub=0;
    _c* buf_lin=0;
    o_seg* sP = segs;
    _u slen;
    _ul slen_wh, slen_sb0;

    for( ; sP->next; sP=sP->next){
        
        if(sP->runcode == syscode){
            system(sP->runseg);
        }
        else{
            printf("[[%s]]\n", sP->runseg);
        }
    }
}

static _v onyx_sysexecs(o_seg* segs, _ul syscode)
{
   
    o_seg* sP = segs;
    
    for( ; sP->next; sP=sP->next){
        if(sP->runcode == syscode)
            system(sP->runseg);
    }
}

static _v onyx_dumpsegs(o_seg* segs)
{
   
    o_seg* sP = segs;
    
    for( ; sP->next; sP=sP->next){
        printf("%s\n", sP->runseg);
    }
}

static _v onyx_dumpin2buf(o_seg* segs, _c* buf, _ul maxlen)
{
    _ul i0;
    o_seg* sP = segs;
    _c* bP = buf;
    for( ; sP->next; sP=sP->next){
        sprintf(buf, "%s\n", sP->runseg);
    }
    printf("%s",buf);
}

/*



lineparse functionality



*/


_v line_initptrs(linenode* obj, linenode* head, linenode* prev, linenode* next)
{
    obj->head = head;
    obj->prev = prev;
    obj->next = next;
}

linenode* line_tail(linenode* l)
{
    linenode* utl=l;
    for(; utl->next; utl=utl->next);
    re(utl);
}

_ul line_lenfrom(_c* lsz)
{
    _ul l = 0;
    _c* lP = lsz;
    for(; *lP != '\n' && *lP; ++lP)l++;
    return(l);
}

_ul line_count(_c* lsz)
{
    _ul l = 0;
    _c* lP = lsz;
    for( ; *lP; ++lP){
        if(*lP == '\n'){
            l++;
        }
    }
    return(l);
}

_ul line_destroyset(linenode* lset)
{
    _ul dlen;
    linenode* utl = line_tail(lset);
    linenode* utlcp = 0;
    for(dlen = 0; ; ){
        printf("about to free:%s\n",utl->line);
        utlcp = (utl->prev != 0) ? utl->prev : 0;
        dlen += szlen2(utl->line);
        free(utl->line);
        free(utl);
        if(utlcp)
            utl = utlcp;
        else
            break;
    }
    re(dlen);
}

linenode* line_break(_c* sz)
{
    _ul lc;
    _ul i0, i1, li;
    _ul done=0;

    _c* szP=sz;
    linenode* lines = line_loke();
    linenode* lptr=0;

    lines->index=0; lines->state=1;
    line_initptrs(lines, lines, 0, 0);
    lptr = lines;

    li=0;
    for(;;){

        for(szP=sz+li; *szP != '\n' && *szP; ++szP)i1++;

        lptr->line = (_c*)malloc( (sizeof(_c)*i1)+1);

        for(szP=sz+li, lc=0; lc < i1; ++szP){
            if(*szP == '\n'){ li++; break;}
            if(*szP == '\0'){ done=1; break;}
            lptr->line[lc] = *szP;
            lc++;
         
        }
        if(done){break;}
        lptr->line[lc]=0;
        li += lc;

        lptr->next = line_loke();
        line_initptrs(lptr->next, lines, lptr, 0);
 
        lptr->next->index = lptr->index + 1;
        lptr=lptr->next;
    }

    return(lines);
}


/*

       wordnode     functionality.....

*/





wordn* word_break(linenode* line)
{
    _i done=0;
    _ul i0, i1, wi;

    wordn* w = wordloke();
    wordn* wptr = w;

    _c* wP = line->line;
    
    for(wi=i0=i1=0; ; ){

        for(wP=line->line + wi;  *wP != ' ' && *wP; ++wP)i1++;
        wptr->word = (_c*)malloc((sizeof(_c) * i1)+1);
        for(wP=line->line + wi,i0=0 ; ; ++wP){
            if(*wP == ' '){ break; }
            if(*wP == '\n'){ i0++; break; }
            if(*wP == '\0'){ done=1; break; }
            wptr->word[i0] = *wP;
            i0++;
        }
        wptr->wordlen = i0;

        for(; *wP == ' '; ++i0);
        if(done){break;}

        wptr->next = wordloke();
        wptr->next->head = w;
        wptr->next->prev = wptr;
        wptr->next->next = 0;

        wptr->next->index = wptr->index + 1;
        wptr->next->state = 1;
        wptr->next->type = 1;

        wptr = wptr->next;

        wi += i0;
    }
    re(w);
}

/*
              wordnode functionality [wnode and lnode)

*/
#define onyx_vartype_str    0
#define onyx_vartype_num    1
#define onyx_vartype_mpn    2

#define onyx_vartype_cmd    3
#define onyx_vartype_i64    4

static _c* gl_delims = { ' ', '\n', '\r', '\0', '~' };

static _ul onyx_is_num(_c* s, _ul base)
{
    _c* sp = s;
    _ul is_dig = 0;
    for( ; *sp; ){
        if( *sp > '0' && *sp < base ){
            is_dig++;
            sp++;
        }else{
            is_dig = 0;
            break;
        }
    }
    return(is_dig);
}



/*

        wnode      functionality ...

*/



wnode* wnode_alloc()
{
    return( (wnode*)malloc(sizeof(wnode)) );
}

wnode* wnew()
{
    wnode* wObj = wnode_alloc();
    wnode_attribs(wObj, 0, 0, 0, 0);
    wnode_ptrs(wObj, 0, 0, 0, 0, 0);
    return(wObj);
}

wnode* wnewx(_ul c)
{
    wnode* wObj = (wnode*)malloc(sizeof(wnode));
    wnode* u = wObj;
    
    wnode_ptrs( wObj, wObj, wObj, 0, 0, 0);
        
        for(u->index = 1 ; u->index < c; ){
            u->w = 0;
            u->wlen = 0;
            u->wtype = 1;
            
            wnode_attribs(u, u->index, 0, 0, 1);
            
            u->next = (wnode*)malloc(sizeof(wnode));
            u->next->index = u->index + 1;
            wnode_ptrs(u->next, wObj, wObj, u, 0, 0);    
            u = u->next;
        }

    return(wObj);
}

wnode* onyx_wend(wnode* wObj)
{
    wnode* u = wObj;
    for(; u->next; u=u->next);
    return(u);
}

#define wend(A)    onyx_wend(A)

wnode* wnode_sz(_c* szin )
{
    _ul i0;
    
    _c* t = 0;
    _c* n = 0;
    
    wnode* w = wnew();
    wnode* w_a = w;
    
    wnode_ptrs(w, w_a, w_a, 0, 0, 0);
    wnode_attribs(w, 1, 0, 0, 13);
    
    for( i0 = 1; ; ++i0){
        n = baseX_digitsum( i0, 10 );
        w->w = word_isolat3( szin, n, 0);
        w->wlen = szlen2( w->w );
        free(n);
        
        if(w->w[0] == '\0')
            break;
        
        w->next = wnew();
        
        wnode_ptrs( w->next, w_a, w_a, w, 0, 0);
        wnode_attribs(w->next, w->index+1, 0, 0, 13);

        w = w->next;
    }
    return(w_a);
}

wnode* wnode_sz2(_c* szin )
{
    _ul i0;
    
    _c* t = 0;
    _c* n = 0;
    
    wnode* w = wnew();
    wnode* w_a = w;
    
    wnode_ptrs(w, w_a, w_a, 0, 0, 0);
    wnode_attribs(w, 1, 0, 0, 13);
    
    for( i0 = 1; ; ++i0){
        n = baseX_digitsum( i0, 10 );
        w->w = word_isolat4( szin, n, 0, gl_delims);
        free(n);
        
        if(w->w[0] == '\0')
            break;
        
        w->next = wnew();
        
        wnode_ptrs( w->next, w_a, w_a, w, 0, 0);
        wnode_attribs(w->next, w->index+1, 0, 0, 13);

        w = w->next;
    }
    return(w_a);
}
wnode* wnode_sz5(_c* szin )
{
    _ul i0;

    _c* t = 0;
    _c* n = 0;

    wnode* w = wnew();
    wnode* w_a = w;

    wnode_ptrs(w, w_a, w_a, 0, 0, 0);
    wnode_attribs(w, 1, 0, 0, 13);

    for( i0 = 1; ; ++i0){
        n = baseX_digitsum( i0, 10 );
        w->w = word_isolat5( szin, n, 0, ';');
        free(n);

        if(w->w[0] == '\0')
            break;

        w->next = wnew();

        wnode_ptrs( w->next, w_a, w_a, w, 0, 0);
        wnode_attribs(w->next, w->index+1, 0, 0, 13);

        w = w->next;
    }
    return(w_a);
}
_v wdestroy(wnode* w)
{
    wnode* w_a = w;
    wnode* las = w;
    for( ; w_a->next; w_a = w_a->next){
 //       printf("%i[%i-%i] %s\n", w_a->index, w_a->state, w_a->archetype, w_a->w);
        free(las->w);
        free(las);
        las = w_a;
    }
    free(las->w);
    free(las);
}

_v wprint(wnode* w, _c delm)
{
    wnode* wp = w;
    for(; wp->next; wp=wp->next){
        printf("%s", wp->w);
        if(delm)
            printf("%c",delm);
        else
            printf(" ");
    }
    printf("\n");
}

wnode* wcopy_to( wnode* w0, wnode* w1 )
{
    wnode* wnod = wnew();
    wnode* wcpy = wnod;

    wnode* wutl = w0;
    for( ; w0 != w1; w0 = w0->next ){
        wnod->wlen = w0->wlen ? : szlen( w0->w );
        wnod->wtype = w0->wtype ? : 1;
        wnod->w =  szalloc( (wnod->wlen > 2) ? : 5 );
        szcpy( wnod->w, (w0->w) ? : ("nil\0") );

        if( w0->next ){
            wnod->next = wnew();
            wnod->next->index = wnod->index + 1;
            wnod->next->prev = wnod;
            wnod->next->head = w0->head;
            wnod->next->next = 0;

            wnod = wnod->next;
        }else{
            break;
        }
    }
    return( wcpy );
}

_ul wcount_to(wnode* w0, wnode* w1)
{
    wnode* wutl = w0;
    for( ; w0 != w1; w0 = w0->next );
    return( w0->index );
}


_ul wcount(wnode* w)
{
    return( wend(w)->index );
}



/*

              lnode   and   dnode    stuff...

*/




lnode* lnode_alloc()
{
    return( (lnode*)malloc(sizeof(lnode)) );
}

lnode* lnew()
{
    lnode* wObj = lnode_alloc();
    lnode_attribs(wObj, 0, 0, 0, 0);
    lnode_ptrs(wObj, 0, 0, 0, 0, 0);
    return(wObj);
}

_i main2(_i argc, _c** argv )
{
    _c* szin = fget_data( argv[1] );
    
    wnode* w = wnode_sz(szin);
    
    wprint(w, 0);
    
    wdestroy(w);
    
    return(0);
}

lnode* lnode_sz(_c* insz)
{
    _ul i0, i1;
    _c* n;
    
    lnode* lObj = lnew();
    lnode* l = lObj;
    
    lnode_ptrs(l, lObj, lObj, 0, 0, 0);
    lnode_attribs(l, 1, 0, 0, 13);
    
    for( i0 = 1 ; ; ++i0 ){
    
        n = baseX_digitsum( i0, 10 );
        l->l = line_isolat3( insz, n, 0, '\n');
        free(n);
        
        if(l->l[0] == '\0')
            break;
            
        l->next = lnew();
        
        lnode_ptrs( l->next, lObj, lObj, l, 0, 0);
        lnode_attribs(l->next, l->index+1, 0, 0, 13);
        
        l = l->next;
    }
    
    return(lObj);
       
}

lnode* lnode_sz2(_c* insz, _c dlm)
{
    _ul i0, i1;
    _c* n;

    lnode* lObj = lnew();
    lnode* l = lObj;

    lnode_ptrs(l, lObj, lObj, 0, 0, 0);
    lnode_attribs(l, 1, 0, 0, 13);

    for( i0 = 1 ; ; ++i0 ){

        n = baseX_digitsum( i0, 10 );
        l->l = line_isolat3( insz, n, 0, dlm);
        free(n);

        if(l->l[0] == '\0')
            break;

        l->next = lnew();

        lnode_ptrs( l->next, lObj, lObj, l, 0, 0);
        lnode_attribs(l->next, l->index+1, 0, 0, 13);

        l = l->next;
    }

    return(lObj);

}

lnode* lnode_sz3(_c* insz, _c dlm)
{
    _ul i0, i1;
    _c* n;

    lnode* lObj = lnew();
    lnode* l = lObj;

    lnode_ptrs(l, lObj, lObj, 0, 0, 0);
    lnode_attribs(l, 1, 0, 0, 13);

    for( i0 = 1 ; ; ++i0 ){

        n = baseX_digitsum( i0, 10 );
        l->l = line_isolat3( insz, n, 0, dlm);
        free(n);

        if(l->l[0] == '\0')
            break;

        l->next = lnew();

        lnode_ptrs( l->next, lObj, lObj, l, 0, 0);
        lnode_attribs(l->next, l->index+1, 0, 0, 13);

        l = l->next;
    }

    return(lObj);

}

lnode* onyx_lend(lnode* l)
{
    lnode* lu = l;
    for(; lu->next; lu = lu->next);
    return(lu);
}
#define lend(A)    onyx_lend(A)

_ul onyx_lcount(lnode* l)
{
    return( lend(l)->index );
}
#define lcount(A) onyx_lcount(A)

lnode* lnode_db(_c* din)
{
    lnode* lObj = lnode_sz( din );
    
    lnode* l = lObj;
    wnode* w=0;
    for(; l->next; l=l->next){
        l->wObj = wnode_sz(l->l);
    }
    return(lObj);
}

lnode* lnode_db2(_c* din, _c ln_dlm)
{
    lnode* lObj = lnode_sz2( din, ln_dlm );
    
    lnode* l = lObj;
    wnode* w=0;
    for(; l->next; l=l->next){
        l->wObj = wnode_sz(l->l);
    }
    return(lObj);
}

lnode* lnode_db3(_c* din, _c ln_dlm)
{
    _ul i0;
    _c* lnum;

    lnode* lObj = lnew(); //lnode_sz2( din, ln_dlm );
    lnode* l = lObj;

    lObj->llen = szlen2(din);
    lObj->l = szalloc(lObj->llen);
    szcpy(lObj->l, din);

    for(i0 = 1 ; ; ++i0){

        lnum = baseX_digitsum(i0, 10);

        lObj->next = lnew();
        lObj->next->l = line_isolat3(l->l, lnum, 0, ln_dlm );

        if( lObj->next->l[0] == '\0' || !lObj->next->l )
            break;

        lObj->next->wObj = wnode_sz( lObj->next->l );

        lObj = lObj->next;

        free(lnum);
    }
    return(l);
}

static _v wwalk(wnode* w)
{ w = w->next?w->next:w; }

wnode* wwalkr(wnode* w)
{ return( w->next ?w->next : 0 ); }


static _v lwalk(lnode* l)
{ l = l->next?l->next:l; }

wnode* lwalkr(lnode* l)
{ return( l->next ?l->next : 0 ); }


_ul ldestroy(lnode* l)
{
    _ul x, y, z;
    _ul xl, yl, zl;
    lnode* l_a = l;
    lnode* las = l;
    
    lnode* l0, l1;
    wnode* w0, *w1;
    
    x = y = z = 0;    
    xl=yl=zl=0;
    for( ; l_a->next; ){
        y = x = 0;
        for( w0 = l_a->wObj, w1 = w0; w0->next; ){
            x = szlen2( w0->w );
            
            w1 = w0;
            w1 = wwalkr(w1);
            free(w0);
            w0 = w1;
            
            y += x;
            

        }
        z += y;
        xl = szlen2(l_a->l);
        yl += xl;

        
        l0 = l_a;
        l0 = lwalkr(l0);
        free(l_a->l);
        free(l_a);
        l_a = l0;
        
        zl += yl;
        zl += z;

    }
    return(zl);
}

dnode* dnew(_c* source_text)
{
    dnode* dnod = (dnode*)malloc(sizeof(dnode));
    dnod->s_len = szlen2(source_text);
    dnod->sz = szalloc( dnod->s_len );
    szcpy(dnod->sz, source_text );

    dnod->word = wnode_sz( source_text );

    return(dnod);
}

dnode* db(_c* din, _c dlm)
{
    _ul i0;
    _ul llen;
    dnode* dnod = dnew(din);
    _c* num=0;



    dnod->line = lnew();
    lnode* lnod = dnod->line;

    lnod->index = 1;

    for( i0 = 0; ; ++i0 ){
        num = baseX_digitsum( i0, 10 );
        lnod->l = line_isolat3( dnod->sz, num, 0, dlm);
        free(num);



        if( !lnod->l )
            break;

        lnod->wObj = wnode_sz( lnod->l );

        lnod->next = lnew();
        lnod->next->index = lnod->index + 1;
        lnod = lnod->next;
    }


    return(dnod);
}

_c* chomp_after_delim(_c* sz, _c dlm, _c replc )
{
   _ul i0;
   _c* sz0, *sz1;
   _ul cutoff=  0;
   _ul whole = szlen2( sz );
   for( sz0 = sz, i0 = 0; *sz0 != dlm && *sz0; ++sz0)i0++;

   sz1 = szalloc( i0 + 1 );
   cutoff = i0;
   for( sz0 = sz, i0 = 0; i0 < cutoff; ++i0, ++sz0){
       sz1[ i0 ] = *sz0;
   }

   cutoff = whole - i0;

   *sz0++;
   *sz0 = replc;

   for( i0 = cutoff-1; i0 < whole-1; ++i0 )
       sz1[ i0] = replc;

   return(sz1);
}





























/*-----------__--__-__--_--__-__----End Of File----------------------------*/
