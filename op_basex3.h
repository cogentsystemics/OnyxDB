

//#include "op_baseXZ.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <unixcube/op_types.h>
//#include <unixcube/op_cbase.h>

static _ui n_digz;
static _uc convbuf[1024];

static _uc baseX_table[10] = {'0', '1', '2', '3', '4','5','6','7','8','9'};

static _uc baseX_table_61[] = {'0', '1', '2', '3', '4','5','6','7','8','9',
                               'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
                               'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r',
                               's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A',
                               'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',
                               'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
                               'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '~', '!',
                               '@', '#', '$', '%', '^', '&', '*', '(', ')',
                               '_', '+', '-', '=', '`', '[', ']', ' ', ';',
                               '"', ',', '.', '/', '{', '}', '|', ':', '"',
                               '<', '>', '?'}; 
#define baseX_ch(X)   baseX_table_64[X]
#define baseX_table_64    baseX_table_61
                               
_ui baseX_digit(_uc c, _i base)
{
	_ui num;
	for(num=0; num <= base; ++num){
		if(c == baseX_table_64[num])
		    break;
	}
	return(num);
}

_ul baseX_sz(_uc* n_sz, _i base)
{
    _ui dig;
	_ul num;
	_ui i0[2];
	_ui nsz_len = strlen((_cc*)n_sz);

	for(i0[0] = 0, dig=0, num=0; i0[0] <= nsz_len-1; ++i0[0]){
		dig = (_ul)baseX_digit(n_sz[i0[0]], base);
		num *= base;
		num += dig;
	}
	return(num);
}

_ul baseX_sum_peelback(_ul sum, _ui base)
{
    _i i0=0;
    _ul finl=0;
    for( finl = sum; finl >= base ; finl -= base)
        i0++;

    return(finl);
}

_ul baseX_digitcount(_ul sum, _ui base)
{
    _i i0=0;
    _ul finl=0;
    for( finl = sum; finl >= base ; finl -= base)
        i0++;

    return(i0++);
}

_c* reverseZ(char * str )
{ 
    int i, len ;
      
    // calculating length of the string
    len = strlen(str);
    _c* re = (_c*)malloc((sizeof(_c) * len)+1);
      
    for( i=len ; i>=0 ; i--)
        re[i] = *(str+i);
    
    return(re);
}

_uc* baseX_digitsum(_ul num, _ui base)
{
    _i i0, i1;
    _ui ops;
	_ui count = 0;
	_ul num_tmp;
	_ul num_tot = num;
	_ul num_prev = num_tot;
    _uc* strrep = 0;

    num_tmp = (_ul)baseX_sum_peelback(num_tot, base);

    convbuf[count] = baseX_table_64[num_tmp];
  
	for(ops = 1, count = 1; num_tot >= base; ++ops, ++count){
	        
		num_tot /= base;
	    num_tmp = (_ul)baseX_sum_peelback(num_tot, base);
        convbuf[count] = baseX_table_64[ (_ui)num_tmp];
    }

    n_digz = ops;
    
    strrep = (_uc*)malloc(sizeof(_uc) * (n_digz + 1));
    
    for( i0=0, i1=n_digz-1; i0 <= n_digz; ++i0){
        strrep[i0] = convbuf[i1];
        i1--;
    }
    strrep[n_digz] = '\0';
    reverseZ(strrep);
	return((_uc*)reverseZ((_uc*)strrep));
}

_ul base3(_ul n, _ul newbase)
{
    _uc* nsz = baseX_digitsum(n, 10);
    _ul nbr = baseX_sz(nsz, newbase);

    free(nsz);
    return(nbr);
}

_ul baseXb_N(_ul n, _ul base0, _ul base1)
{
    _uc* nsz = baseX_digitsum(n, base0);
    _ul nbr = baseX_sz(nsz, base1);

    free(nsz);
    return(nbr);
}

_uc* baseXb_ch(_c* nch, _c* base0, _c* base1)
{
    _ul b0_n = baseX_sz(base0, 10);
    _ul b1_n = baseX_sz(base1, 10);
    _ul r0_n = baseX_sz(nch, base0);

    _ul r_n = baseXb_N(r0_n, b0_n, b1_n);

    _uc* rch = baseX_digitsum(r_n, b1_n);

    return(rch);
}

typedef struct basex_t{
    _c* bx_sz;
    _ul bx_n;
    
    _ul bx_szlen_org, bx_szlen_conv;
    _ul bx_bits_org, bx_bits_conv;
    _ul bx_digits_org, bx_digits_conv;
    _ul bx_base_org, bx_base_conv;
    
    struct basex_t* root, *prev, *next, *tail;
    struct basex_t* branch;
    
    _ul bx_indx, bx_stat;
}basex;

basex* basex_alloc()
{
    basex* r = (basex*)malloc(sizeof(basex) );
    memset(r, 0, sizeof(basex));
    r->bx_stat = 1;
    return(r);
}

_v basex_init_ptrs(basex* node, basex* rootnode, basex* pre, basex* nex, basex* end, basex* x)
{
    node->root = rootnode; node->prev=pre; node->next = nex; node->tail=end;
    node->branch = x;
}

basex* basex_tail(basex* node)
{
  
}

basex* basex_n_set(_ul n, _ul base, basex* headnode)
{
    basex* r = basex_alloc();
    basex* util=0;
    r->bx_n = n;
    r->bx_base_org = base;
    r->bx_sz = baseX_digitsum(r->bx_n, r->bx_base_org);
    r->bx_szlen_org = strlen((_cc*)r->bx_sz);
    r->bx_digits_org = baseX_digitcount(r->bx_n, r->bx_base_org);
    
    if(headnode){
        util = basex_tail(headnode);
	r->bx_indx = util->bx_indx + 1;
	basex_init_ptrs(r, headnode, util, 0, r, 0);
	basex_init_ptrs(util, util->root, util->prev, r, r, 0);
    }else{
        r->bx_indx = 0;
        basex_init_ptrs(r, r, 0, 0, r, 0);
    }
    
    r->bx_stat = 2;
    
    return(r);
}

basex* basex_sz_set(_ul n, _c* sz, _ul base, basex* headnode)
{
    basex* r = basex_alloc();
    basex* util=0;
    
    if(!n)
        r->bx_n = baseX_sz(sz, base);
    else
        r->bx_n = n;
    
    r->bx_base_org = base;
    r->bx_sz = baseX_digitsum(r->bx_n, r->bx_base_org);
    r->bx_szlen_org = strlen((_cc*)r->bx_sz);
    r->bx_digits_org = baseX_digitcount(r->bx_n, r->bx_base_org);
    
    if(headnode){
        util = basex_tail(headnode);
	r->bx_indx = util->bx_indx + 1;
	basex_init_ptrs(r, headnode, util, 0, r, 0);
	basex_init_ptrs(util, util->root, util->prev, r, r, 0);
    }else{
        r->bx_indx = 0;
        basex_init_ptrs(r, r, 0, 0, r, 0);
    }
    
    r->bx_stat = 2;
    
    return(r);
}

basex* basex_convert_sz1(_c* szn, _c* szb1, _c* szb2)
{
    _ul base0 = baseX_sz(szb1, 10);
    _ul base1 = baseX_sz(szb2, 10);
 
    basex* bas0 = basex_sz_set(0, szn, base0, 0);
    basex* bas1 = basex_sz_set(bas0->bx_n, 0, base1, 0);
    
    free(bas0->bx_sz);
    free(bas0);

    return(bas1);
}

//00000000000000000000000000000000000000000000000000000000000000000000000000000000000
