#include <string.h>
#include <unistd.h>

//  This is a concatenation of bits from samtools 0.1.19 needed to read sam/bam files.
//
//  The C files are in alpha order, the header files are in (mostly) include order
//
//  cat sam.c bam.c bam_aux.c bam_import.c bgzf.c sam_header.c                > bamcat.c
//  cat bgzf.h bam.h sam.h bam_endian.h kstring.h sam_header.h khash.h kseq.h > bamcat.h
//
//  Then replace all (#include ") with (//#include "), and include just the single bam.h.
//  A few functions were removed to cut down on the amount of crud.
//
//  BPW - Fri Feb 20 16:14:11 EST 2015

#include "bamcat.h"

//#include "sam.h"

#define TYPE_BAM  1
#define TYPE_READ 2

bam_header_t *bam_header_dup(const bam_header_t *h0)
{
	bam_header_t *h;
	int i;
	h = bam_header_init();
	*h = *h0;
	h->hash = h->dict = h->rg2lib = 0;
	h->text = (char*)calloc(h->l_text + 1, 1);
	memcpy(h->text, h0->text, h->l_text);
	h->target_len = (uint32_t*)calloc(h->n_targets, 4);
	h->target_name = (char**)calloc(h->n_targets, sizeof(void*));
	for (i = 0; i < h->n_targets; ++i) {
		h->target_len[i] = h0->target_len[i];
		h->target_name[i] = strdup(h0->target_name[i]);
	}
	return h;
}
static void append_header_text(bam_header_t *header, char* text, int len)
{
	int x = header->l_text + 1;
	int y = header->l_text + len + 1; // 1 byte null
	if (text == 0) return;
	kroundup32(x); 
	kroundup32(y);
	if (x < y) header->text = (char*)realloc(header->text, y);
	strncpy(header->text + header->l_text, text, len); // we cannot use strcpy() here.
	header->l_text += len;
	header->text[header->l_text] = 0;
}

samfile_t *samopen(const char *fn, const char *mode, const void *aux)
{
	samfile_t *fp;
	fp = (samfile_t*)calloc(1, sizeof(samfile_t));
	if (strchr(mode, 'r')) { // read
		fp->type |= TYPE_READ;
		if (strchr(mode, 'b')) { // binary
			fp->type |= TYPE_BAM;
			fp->x.bam = strcmp(fn, "-")? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r");
			if (fp->x.bam == 0) goto open_err_ret;
			fp->header = bam_header_read(fp->x.bam);
		} else { // text
			fp->x.tamr = sam_open(fn);
			if (fp->x.tamr == 0) goto open_err_ret;
			fp->header = sam_header_read(fp->x.tamr);
			if (fp->header->n_targets == 0) { // no @SQ fields
				if (aux) { // check if aux is present
					bam_header_t *textheader = fp->header;
					fp->header = sam_header_read2((const char*)aux);
					if (fp->header == 0) goto open_err_ret;
					append_header_text(fp->header, textheader->text, textheader->l_text);
					bam_header_destroy(textheader);
				}
				if (fp->header->n_targets == 0 && bam_verbose >= 1)
					fprintf(stderr, "[samopen] no @SQ lines in the header.\n");
			} else if (bam_verbose >= 2) fprintf(stderr, "[samopen] SAM header is present: %d sequences.\n", fp->header->n_targets);
		}
	} else if (strchr(mode, 'w')) { // write
    fprintf(stderr, "Write not supported here.\n");
    exit(1);
  }

	return fp;

open_err_ret:
	free(fp);
	return 0;
}

void samclose(samfile_t *fp)
{
	if (fp == 0) return;
	if (fp->header) bam_header_destroy(fp->header);
	if (fp->type & TYPE_BAM) bam_close(fp->x.bam);
	else if (fp->type & TYPE_READ) sam_close(fp->x.tamr);
	else fclose(fp->x.tamw);
	free(fp);
}

int samread(samfile_t *fp, bam1_t *b)
{
	if (fp == 0 || !(fp->type & TYPE_READ)) return -1; // not open for reading
	if (fp->type & TYPE_BAM) return bam_read1(fp->x.bam, b);
	else return sam_read1(fp->x.tamr, fp->header, b);
}
#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
//#include "bam.h"
//#include "bam_endian.h"
//#include "kstring.h"
//#include "sam_header.h"

int bam_is_be = 0, bam_verbose = 2, bam_no_B = 0;
char *bam_flag2char_table = "pPuUrR12sfd\0\0\0\0\0";




/**************************
 * CIGAR related routines *
 **************************/

uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar)
{
	int k, end = c->pos;
	for (k = 0; k < c->n_cigar; ++k) {
		int op  = bam_cigar_op(cigar[k]);
		int len = bam_cigar_oplen(cigar[k]);
		if (op == BAM_CBACK) { // move backward
			int l, u, v;
			if (k == c->n_cigar - 1) break; // skip trailing 'B'
			for (l = k - 1, u = v = 0; l >= 0; --l) {
				int op1  = bam_cigar_op(cigar[l]);
				int len1 = bam_cigar_oplen(cigar[l]);
				if (bam_cigar_type(op1)&1) { // consume query
					if (u + len1 >= len) { // stop
						if (bam_cigar_type(op1)&2) v += len - u;
						break;
					} else u += len1;
				}
				if (bam_cigar_type(op1)&2) v += len1;
			}
			end = l < 0? c->pos : end - v;
		} else if (bam_cigar_type(op)&2) end += bam_cigar_oplen(cigar[k]);
	}
	return end;
}

int32_t bam_cigar2qlen(const bam1_core_t *c, const uint32_t *cigar)
{
	uint32_t k;
	int32_t l = 0;
	for (k = 0; k < c->n_cigar; ++k)
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			l += bam_cigar_oplen(cigar[k]);
	return l;
}

/********************
 * BAM I/O routines *
 ********************/

bam_header_t *bam_header_init()
{
	bam_is_be = bam_is_big_endian();
	return (bam_header_t*)calloc(1, sizeof(bam_header_t));
}

void bam_header_destroy(bam_header_t *header)
{
	int32_t i;
	//extern void bam_destroy_header_hash(bam_header_t *header);
	if (header == 0) return;
	if (header->target_name) {
		for (i = 0; i < header->n_targets; ++i)
			free(header->target_name[i]);
		free(header->target_name);
		free(header->target_len);
	}
	free(header->text);
	if (header->dict) sam_header_free(header->dict);
	if (header->rg2lib) sam_tbl_destroy(header->rg2lib);
	//bam_destroy_header_hash(header);
	free(header);
}

bam_header_t *bam_header_read(bamFile fp)
{
	bam_header_t *header;
	char buf[4];
	int magic_len;
	int32_t i = 1, name_len;
	// check EOF
	i = bgzf_check_EOF(fp);
	if (i < 0) {
		// If the file is a pipe, checking the EOF marker will *always* fail
		// with ESPIPE.  Suppress the error message in this case.
		if (errno != ESPIPE) perror("[bam_header_read] bgzf_check_EOF");
	}
	else if (i == 0) fprintf(stderr, "[bam_header_read] EOF marker is absent. The input is probably truncated.\n");
	// read "BAM1"
	magic_len = bam_read(fp, buf, 4);
	if (magic_len != 4 || strncmp(buf, "BAM\001", 4) != 0) {
		fprintf(stderr, "[bam_header_read] invalid BAM binary header (this is not a BAM file).\n");
		return 0;
	}
	header = bam_header_init();
	// read plain text and the number of reference sequences
	bam_read(fp, &header->l_text, 4);
	if (bam_is_be) bam_swap_endian_4p(&header->l_text);
	header->text = (char*)calloc(header->l_text + 1, 1);
	bam_read(fp, header->text, header->l_text);
	bam_read(fp, &header->n_targets, 4);
	if (bam_is_be) bam_swap_endian_4p(&header->n_targets);
	// read reference sequence names and lengths
	header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
	header->target_len = (uint32_t*)calloc(header->n_targets, 4);
	for (i = 0; i != header->n_targets; ++i) {
		bam_read(fp, &name_len, 4);
		if (bam_is_be) bam_swap_endian_4p(&name_len);
		header->target_name[i] = (char*)calloc(name_len, 1);
		bam_read(fp, header->target_name[i], name_len);
		bam_read(fp, &header->target_len[i], 4);
		if (bam_is_be) bam_swap_endian_4p(&header->target_len[i]);
	}
	return header;
}

static void swap_endian_data(const bam1_core_t *c, int data_len, uint8_t *data)
{
	uint8_t *s;
	uint32_t i, *cigar = (uint32_t*)(data + c->l_qname);
	s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
	for (i = 0; i < c->n_cigar; ++i) bam_swap_endian_4p(&cigar[i]);
	while (s < data + data_len) {
		uint8_t type;
		s += 2; // skip key
		type = toupper(*s); ++s; // skip type
		if (type == 'C' || type == 'A') ++s;
		else if (type == 'S') { bam_swap_endian_2p(s); s += 2; }
		else if (type == 'I' || type == 'F') { bam_swap_endian_4p(s); s += 4; }
		else if (type == 'D') { bam_swap_endian_8p(s); s += 8; }
		else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
		else if (type == 'B') {
			int32_t n, Bsize = bam_aux_type2size(*s);
			memcpy(&n, s + 1, 4);
			if (1 == Bsize) {
			} else if (2 == Bsize) {
				for (i = 0; i < n; i += 2)
					bam_swap_endian_2p(s + 5 + i);
			} else if (4 == Bsize) {
				for (i = 0; i < n; i += 4)
					bam_swap_endian_4p(s + 5 + i);
			}
			bam_swap_endian_4p(s+1); 
		}
	}
}

int bam_read1(bamFile fp, bam1_t *b)
{
	bam1_core_t *c = &b->core;
	int32_t block_len, ret, i;
	uint32_t x[8];

	assert(BAM_CORE_SIZE == 32);
	if ((ret = bam_read(fp, &block_len, 4)) != 4) {
		if (ret == 0) return -1; // normal end-of-file
		else return -2; // truncated
	}
	if (bam_read(fp, x, BAM_CORE_SIZE) != BAM_CORE_SIZE) return -3;
	if (bam_is_be) {
		bam_swap_endian_4p(&block_len);
		for (i = 0; i < 8; ++i) bam_swap_endian_4p(x + i);
	}
	c->tid = x[0]; c->pos = x[1];
	c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
	c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
	c->l_qseq = x[4];
	c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
	b->data_len = block_len - BAM_CORE_SIZE;
	if (b->m_data < b->data_len) {
		b->m_data = b->data_len;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	if (bam_read(fp, b->data, b->data_len) != b->data_len) return -4;
	b->l_aux = b->data_len - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
	if (bam_is_be) swap_endian_data(c, b->data_len, b->data);
	if (bam_no_B) bam_remove_B(b);
	return 4 + block_len;
}


/************
 * Remove B *
 ************/

int bam_remove_B(bam1_t *b)
{
	int i, j, end_j, k, l, no_qual;
	uint32_t *cigar, *new_cigar;
	uint8_t *seq, *qual, *p;
	// test if removal is necessary
	if (b->core.flag & BAM_FUNMAP) return 0; // unmapped; do nothing
	cigar = bam1_cigar(b);
	for (k = 0; k < b->core.n_cigar; ++k)
		if (bam_cigar_op(cigar[k]) == BAM_CBACK) break;
	if (k == b->core.n_cigar) return 0; // no 'B'
	if (bam_cigar_op(cigar[0]) == BAM_CBACK) goto rmB_err; // cannot be removed
	// allocate memory for the new CIGAR
	if (b->data_len + (b->core.n_cigar + 1) * 4 > b->m_data) { // not enough memory
		b->m_data = b->data_len + b->core.n_cigar * 4;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
		cigar = bam1_cigar(b); // after realloc, cigar may be changed
	}
	new_cigar = (uint32_t*)(b->data + (b->m_data - b->core.n_cigar * 4)); // from the end of b->data
	// the core loop
	seq = bam1_seq(b); qual = bam1_qual(b);
	no_qual = (qual[0] == 0xff); // test whether base quality is available
	i = j = 0; end_j = -1;
	for (k = l = 0; k < b->core.n_cigar; ++k) {
		int op  = bam_cigar_op(cigar[k]);
		int len = bam_cigar_oplen(cigar[k]);
		if (op == BAM_CBACK) { // the backward operation
			int t, u;
			if (k == b->core.n_cigar - 1) break; // ignore 'B' at the end of CIGAR
			if (len > j) goto rmB_err; // an excessively long backward
			for (t = l - 1, u = 0; t >= 0; --t) { // look back
				int op1  = bam_cigar_op(new_cigar[t]);
				int len1 = bam_cigar_oplen(new_cigar[t]);
				if (bam_cigar_type(op1)&1) { // consume the query
					if (u + len1 >= len) { // stop
						new_cigar[t] -= (len - u) << BAM_CIGAR_SHIFT;
						break;
					} else u += len1;
				}
			}
			if (bam_cigar_oplen(new_cigar[t]) == 0) --t; // squeeze out the zero-length operation
			l = t + 1;
			end_j = j; j -= len;
		} else { // other CIGAR operations
			new_cigar[l++] = cigar[k];
			if (bam_cigar_type(op)&1) { // consume the query
				if (i != j) { // no need to copy if i == j
					int u, c, c0;
					for (u = 0; u < len; ++u) { // construct the consensus
						c = bam1_seqi(seq, i+u);
						if (j + u < end_j) { // in an overlap
							c0 = bam1_seqi(seq, j+u);
							if (c != c0) { // a mismatch; choose the better base
								if (qual[j+u] < qual[i+u]) { // the base in the 2nd segment is better
									bam1_seq_seti(seq, j+u, c);
									qual[j+u] = qual[i+u] - qual[j+u];
								} else qual[j+u] -= qual[i+u]; // the 1st is better; reduce base quality
							} else qual[j+u] = qual[j+u] > qual[i+u]? qual[j+u] : qual[i+u];
						} else { // not in an overlap; copy over
							bam1_seq_seti(seq, j+u, c);
							qual[j+u] = qual[i+u];
						}
					}
				}
				i += len, j += len;
			}
		}
	}
	if (no_qual) qual[0] = 0xff; // in very rare cases, this may be modified
	// merge adjacent operations if possible
	for (k = 1; k < l; ++k)
		if (bam_cigar_op(new_cigar[k]) == bam_cigar_op(new_cigar[k-1]))
			new_cigar[k] += new_cigar[k-1] >> BAM_CIGAR_SHIFT << BAM_CIGAR_SHIFT, new_cigar[k-1] &= 0xf;
	// kill zero length operations
	for (k = i = 0; k < l; ++k)
		if (new_cigar[k] >> BAM_CIGAR_SHIFT)
			new_cigar[i++] = new_cigar[k];
	l = i;
	// update b
	memcpy(cigar, new_cigar, l * 4); // set CIGAR
	p = b->data + b->core.l_qname + l * 4;
	memmove(p, seq, (j+1)>>1); p += (j+1)>>1; // set SEQ
	memmove(p, qual, j); p += j; // set QUAL
	memmove(p, bam1_aux(b), b->l_aux); p += b->l_aux; // set optional fields
	b->core.n_cigar = l, b->core.l_qseq = j; // update CIGAR length and query length
	b->data_len = p - b->data; // update record length
	return 0;

rmB_err:
	b->core.flag |= BAM_FUNMAP;
	return -1;
}
#include <ctype.h>
//#include "bam.h"
//#include "khash.h"
typedef char *str_p;
KHASH_MAP_INIT_STR(s, int)
KHASH_MAP_INIT_STR(r2l, str_p)

void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data)
{
	int ori_len = b->data_len;
	b->data_len += 3 + len;
	b->l_aux += 3 + len;
	if (b->m_data < b->data_len) {
		b->m_data = b->data_len;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	b->data[ori_len] = tag[0]; b->data[ori_len + 1] = tag[1];
	b->data[ori_len + 2] = type;
	memcpy(b->data + ori_len + 3, data, len);
}

uint8_t *bam_aux_get_core(bam1_t *b, const char tag[2])
{
	return bam_aux_get(b, tag);
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += bam_aux_type2size(type); \
	} while(0)

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return s;
		__skip_tag(s);
	}
	return 0;
}
// s MUST BE returned by bam_aux_get()
int bam_aux_del(bam1_t *b, uint8_t *s)
{
	uint8_t *p, *aux;
	aux = bam1_aux(b);
	p = s - 2;
	__skip_tag(s);
	memmove(p, s, b->l_aux - (s - aux));
	b->data_len -= s - p;
	b->l_aux -= s - p;
	return 0;
}

int bam_aux_drop_other(bam1_t *b, uint8_t *s)
{
	if (s) {
		uint8_t *p, *aux;
		aux = bam1_aux(b);
		p = s - 2;
		__skip_tag(s);
		memmove(aux, p, s - p);
		b->data_len -= b->l_aux - (s - p);
		b->l_aux = s - p;
	} else {
		b->data_len -= b->l_aux;
		b->l_aux = 0;
	}
	return 0;
}

void bam_init_header_hash(bam_header_t *header)
{
	if (header->hash == 0) {
		int ret, i;
		khiter_t iter;
		khash_t(s) *h;
		header->hash = h = kh_init(s);
		for (i = 0; i < header->n_targets; ++i) {
			iter = kh_put(s, h, header->target_name[i], &ret);
			kh_value(h, iter) = i;
		}
	}
}

void bam_destroy_header_hash(bam_header_t *header)
{
	if (header->hash)
		kh_destroy(s, (khash_t(s)*)header->hash);
}

int32_t bam_get_tid(const bam_header_t *header, const char *seq_name)
{
	khint_t k;
	khash_t(s) *h = (khash_t(s)*)header->hash;
	k = kh_get(s, h, seq_name);
	return k == kh_end(h)? -1 : kh_value(h, k);
}

int bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *beg, int *end)
{
	char *s;
	int i, l, k, name_end;
	khiter_t iter;
	khash_t(s) *h;

	bam_init_header_hash(header);
	h = (khash_t(s)*)header->hash;

	*ref_id = *beg = *end = -1;
	name_end = l = strlen(str);
	s = (char*)malloc(l+1);
	// remove space
	for (i = k = 0; i < l; ++i)
		if (!isspace(str[i])) s[k++] = str[i];
	s[k] = 0; l = k;
	// determine the sequence name
	for (i = l - 1; i >= 0; --i) if (s[i] == ':') break; // look for colon from the end
	if (i >= 0) name_end = i;
	if (name_end < l) { // check if this is really the end
		int n_hyphen = 0;
		for (i = name_end + 1; i < l; ++i) {
			if (s[i] == '-') ++n_hyphen;
			else if (!isdigit(s[i]) && s[i] != ',') break;
		}
		if (i < l || n_hyphen > 1) name_end = l; // malformated region string; then take str as the name
		s[name_end] = 0;
		iter = kh_get(s, h, s);
		if (iter == kh_end(h)) { // cannot find the sequence name
			iter = kh_get(s, h, str); // try str as the name
			if (iter == kh_end(h)) {
				if (bam_verbose >= 2) fprintf(stderr, "[%s] fail to determine the sequence name.\n", __func__);
				free(s); return -1;
			} else s[name_end] = ':', name_end = l;
		}
	} else iter = kh_get(s, h, str);
        if (iter == kh_end(h)) {
          free(s); 
          return -1;
        }
	*ref_id = kh_val(h, iter);
	// parse the interval
	if (name_end < l) {
		for (i = k = name_end + 1; i < l; ++i)
			if (s[i] != ',') s[k++] = s[i];
		s[k] = 0;
		*beg = atoi(s + name_end + 1);
		for (i = name_end + 1; i != k; ++i) if (s[i] == '-') break;
		*end = i < k? atoi(s + i + 1) : 1<<29;
		if (*beg > 0) --*beg;
	} else *beg = 0, *end = 1<<29;
	free(s);
	return *beg <= *end? 0 : -1;
}

int32_t bam_aux2i(const uint8_t *s)
{
	int type;
	if (s == 0) return 0;
	type = *s++;
	if (type == 'c') return (int32_t)*(int8_t*)s;
	else if (type == 'C') return (int32_t)*(uint8_t*)s;
	else if (type == 's') return (int32_t)*(int16_t*)s;
	else if (type == 'S') return (int32_t)*(uint16_t*)s;
	else if (type == 'i' || type == 'I') return *(int32_t*)s;
	else return 0;
}

float bam_aux2f(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0.0;
	if (type == 'f') return *(float*)s;
	else return 0.0;
}

double bam_aux2d(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0.0;
	if (type == 'd') return *(double*)s;
	else return 0.0;
}

char bam_aux2A(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0;
	if (type == 'A') return *(char*)s;
	else return 0;
}

char *bam_aux2Z(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0;
	if (type == 'Z' || type == 'H') return (char*)s;
	else return 0;
}

#ifdef _WIN32
double drand48()
{
	return (double)rand() / RAND_MAX;
}
#endif
#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#ifdef _WIN32
#include <fcntl.h>
#endif
//#include "kstring.h"
//#include "bam.h"
//#include "sam_header.h"
//#include "kseq.h"
//#include "khash.h"

KSTREAM_INIT(gzFile, gzread, 16384)
KHASH_MAP_INIT_STR(ref, uint64_t)

void bam_init_header_hash(bam_header_t *header);
void bam_destroy_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);

unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

unsigned short bam_char2flag_table[256] = {
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,BAM_FREAD1,BAM_FREAD2,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	BAM_FPROPER_PAIR,0,BAM_FMREVERSE,0, 0,BAM_FMUNMAP,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, BAM_FDUP,0,BAM_FQCFAIL,0, 0,0,0,0, 0,0,0,0,
	BAM_FPAIRED,0,BAM_FREVERSE,BAM_FSECONDARY, 0,BAM_FUNMAP,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
	0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
};

char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

struct __tamFile_t {
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	uint64_t n_lines;
	int is_first;
};

char **__bam_get_lines(const char *fn, int *_n) // for bam_plcmd.c only
{
	char **list = 0, *s;
	int n = 0, dret, m = 0;
	gzFile fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	kstream_t *ks;
	kstring_t *str;
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	ks = ks_init(fp);
	while (ks_getuntil(ks, '\n', str, &dret) > 0) {
		if (n == m) {
			m = m? m << 1 : 16;
			list = (char**)realloc(list, m * sizeof(char*));
		}
		if (str->s[str->l-1] == '\r')
			str->s[--str->l] = '\0';
		s = list[n++] = (char*)calloc(str->l + 1, 1);
		strcpy(s, str->s);
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	*_n = n;
	return list;
}

static bam_header_t *hash2header(const kh_ref_t *hash)
{
	bam_header_t *header;
	khiter_t k;
	header = bam_header_init();
	header->n_targets = kh_size(hash);
	header->target_name = (char**)calloc(kh_size(hash), sizeof(char*));
	header->target_len = (uint32_t*)calloc(kh_size(hash), 4);
	for (k = kh_begin(hash); k != kh_end(hash); ++k) {
		if (kh_exist(hash, k)) {
			int i = (int)kh_value(hash, k);
			header->target_name[i] = (char*)kh_key(hash, k);
			header->target_len[i] = kh_value(hash, k)>>32;
		}
	}
	bam_init_header_hash(header);
	return header;
}
bam_header_t *sam_header_read2(const char *fn)
{
	bam_header_t *header;
	int c, dret, ret, error = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t *str;
	kh_ref_t *hash;
	khiter_t k;
	if (fn == 0) return 0;
	fp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
	if (fp == 0) return 0;
	hash = kh_init(ref);
	ks = ks_init(fp);
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	while (ks_getuntil(ks, 0, str, &dret) > 0) {
		char *s = strdup(str->s);
		int len, i;
		i = kh_size(hash);
		ks_getuntil(ks, 0, str, &dret);
		len = atoi(str->s);
		k = kh_put(ref, hash, s, &ret);
		if (ret == 0) {
			fprintf(stderr, "[sam_header_read2] duplicated sequence name: %s\n", s);
			error = 1;
		}
		kh_value(hash, k) = (uint64_t)len<<32 | i;
		if (dret != '\n')
			while ((c = ks_getc(ks)) != '\n' && c != -1);
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	fprintf(stderr, "[sam_header_read2] %d sequences loaded.\n", kh_size(hash));
	if (error) return 0;
	header = hash2header(hash);
	kh_destroy(ref, hash);
	return header;
}
static inline uint8_t *alloc_data(bam1_t *b, int size)
{
	if (b->m_data < size) {
		b->m_data = size;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	return b->data;
}
static inline void parse_error(int64_t n_lines, const char * __restrict msg)
{
	fprintf(stderr, "Parse error at line %lld: %s\n", (long long)n_lines, msg);
	abort();
}
static inline void append_text(bam_header_t *header, kstring_t *str)
{
	size_t x = header->l_text, y = header->l_text + str->l + 2; // 2 = 1 byte dret + 1 byte null
	kroundup32(x); kroundup32(y);
	if (x < y) 
    {
        header->n_text = y;
        header->text = (char*)realloc(header->text, y);
        if ( !header->text ) 
        {
            fprintf(stderr,"realloc failed to alloc %ld bytes\n", y);
            abort();
        }
    }
    // Sanity check
    if ( header->l_text+str->l+1 >= header->n_text )
    {
        fprintf(stderr,"append_text FIXME: %ld>=%ld, x=%ld,y=%ld\n",  header->l_text+str->l+1,(long)header->n_text,x,y);
        abort();
    }
	strncpy(header->text + header->l_text, str->s, str->l+1); // we cannot use strcpy() here.
	header->l_text += str->l + 1;
	header->text[header->l_text] = 0;
}

int sam_header_parse(bam_header_t *h)
{
	char **tmp;
	int i;
	free(h->target_len); free(h->target_name);
	h->n_targets = 0; h->target_len = 0; h->target_name = 0;
	if (h->l_text < 3) return 0;
	if (h->dict == 0) h->dict = sam_header_parse2(h->text);
	tmp = sam_header2list(h->dict, "SQ", "SN", &h->n_targets);
	if (h->n_targets == 0) return 0;
	h->target_name = calloc(h->n_targets, sizeof(void*));
	for (i = 0; i < h->n_targets; ++i)
		h->target_name[i] = strdup(tmp[i]);
	free(tmp);
	tmp = sam_header2list(h->dict, "SQ", "LN", &h->n_targets);
	h->target_len = calloc(h->n_targets, 4);
	for (i = 0; i < h->n_targets; ++i)
		h->target_len[i] = atoi(tmp[i]);
	free(tmp);
	return h->n_targets;
}

bam_header_t *sam_header_read(tamFile fp)
{
	int ret, dret;
	bam_header_t *header = bam_header_init();
	kstring_t *str = fp->str;
	while ((ret = ks_getuntil(fp->ks, KS_SEP_TAB, str, &dret)) >= 0 && str->s[0] == '@') { // skip header
		str->s[str->l] = dret; // note that str->s is NOT null terminated!!
		append_text(header, str);
		if (dret != '\n') {
			ret = ks_getuntil(fp->ks, '\n', str, &dret);
			str->s[str->l] = '\n'; // NOT null terminated!!
			append_text(header, str);
		}
		++fp->n_lines;
	}
	sam_header_parse(header);
	bam_init_header_hash(header);
	fp->is_first = 1;
	return header;
}

int sam_read1(tamFile fp, bam_header_t *header, bam1_t *b)
{
	int ret, doff, doff0, dret, z = 0;
	bam1_core_t *c = &b->core;
	kstring_t *str = fp->str;
	kstream_t *ks = fp->ks;

	if (fp->is_first) {
		fp->is_first = 0;
		ret = str->l;
	} else {
		do { // special consideration for empty lines
			ret = ks_getuntil(fp->ks, KS_SEP_TAB, str, &dret);
			if (ret >= 0) z += str->l + 1;
		} while (ret == 0);
	}
	if (ret < 0) return -1;
	++fp->n_lines;
	doff = 0;

	{ // name
		c->l_qname = strlen(str->s) + 1;
		memcpy(alloc_data(b, doff + c->l_qname) + doff, str->s, c->l_qname);
		doff += c->l_qname;
	}
	{ // flag
		long flag;
		char *s;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1;
		flag = strtol((char*)str->s, &s, 0);
		if (*s) { // not the end of the string
			flag = 0;
			for (s = str->s; *s; ++s)
				flag |= bam_char2flag_table[(int)*s];
		}
		c->flag = flag;
	}
	{ // tid, pos, qual
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1; c->tid = bam_get_tid(header, str->s);
		if (c->tid < 0 && strcmp(str->s, "*")) {
			if (header->n_targets == 0) {
				fprintf(stderr, "[sam_read1] missing header? Abort!\n");
				exit(1);
			} else fprintf(stderr, "[sam_read1] reference '%s' is recognized as '*'.\n", str->s);
		}
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1; c->pos = isdigit(str->s[0])? atoi(str->s) - 1 : -1;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1; c->qual = isdigit(str->s[0])? atoi(str->s) : 0;
		if (ret < 0) return -2;
	}
	{ // cigar
		char *s, *t;
		int i, op;
		long x;
		c->n_cigar = 0;
		if (ks_getuntil(ks, KS_SEP_TAB, str, &dret) < 0) return -3;
		z += str->l + 1;
		if (str->s[0] != '*') {
			uint32_t *cigar;
			for (s = str->s; *s; ++s) {
				if ((isalpha(*s)) || (*s=='=')) ++c->n_cigar;
				else if (!isdigit(*s)) parse_error(fp->n_lines, "invalid CIGAR character");
			}
			b->data = alloc_data(b, doff + c->n_cigar * 4);
			cigar = bam1_cigar(b);
			for (i = 0, s = str->s; i != c->n_cigar; ++i) {
				x = strtol(s, &t, 10);
				op = toupper(*t);
				if (op == 'M') op = BAM_CMATCH;
				else if (op == 'I') op = BAM_CINS;
				else if (op == 'D') op = BAM_CDEL;
				else if (op == 'N') op = BAM_CREF_SKIP;
				else if (op == 'S') op = BAM_CSOFT_CLIP;
				else if (op == 'H') op = BAM_CHARD_CLIP;
				else if (op == 'P') op = BAM_CPAD;
				else if (op == '=') op = BAM_CEQUAL;
				else if (op == 'X') op = BAM_CDIFF;
				else if (op == 'B') op = BAM_CBACK;
				else parse_error(fp->n_lines, "invalid CIGAR operation");
				s = t + 1;
				cigar[i] = bam_cigar_gen(x, op);
			}
			if (*s) parse_error(fp->n_lines, "unmatched CIGAR operation");
			c->bin = bam_reg2bin(c->pos, bam_calend(c, cigar));
			doff += c->n_cigar * 4;
		} else {
			if (!(c->flag&BAM_FUNMAP)) {
				fprintf(stderr, "Parse warning at line %lld: mapped sequence without CIGAR\n", (long long)fp->n_lines);
				c->flag |= BAM_FUNMAP;
			}
			c->bin = bam_reg2bin(c->pos, c->pos + 1);
		}
	}
	{ // mtid, mpos, isize
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1;
		c->mtid = strcmp(str->s, "=")? bam_get_tid(header, str->s) : c->tid;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1;
		c->mpos = isdigit(str->s[0])? atoi(str->s) - 1 : -1;
		ret = ks_getuntil(ks, KS_SEP_TAB, str, &dret); z += str->l + 1;
		c->isize = (str->s[0] == '-' || isdigit(str->s[0]))? atoi(str->s) : 0;
		if (ret < 0) return -4;
	}
	{ // seq and qual
		int i;
		uint8_t *p = 0;
		if (ks_getuntil(ks, KS_SEP_TAB, str, &dret) < 0) return -5; // seq
		z += str->l + 1;
		if (strcmp(str->s, "*")) {
			c->l_qseq = strlen(str->s);
			if (c->n_cigar && c->l_qseq != (int32_t)bam_cigar2qlen(c, bam1_cigar(b))) {
				fprintf(stderr, "Line %ld, sequence length %i vs %i from CIGAR\n",
						(long)fp->n_lines, c->l_qseq, (int32_t)bam_cigar2qlen(c, bam1_cigar(b)));
				parse_error(fp->n_lines, "CIGAR and sequence length are inconsistent");
			}
			p = (uint8_t*)alloc_data(b, doff + c->l_qseq + (c->l_qseq+1)/2) + doff;
			memset(p, 0, (c->l_qseq+1)/2);
			for (i = 0; i < c->l_qseq; ++i)
				p[i/2] |= bam_nt16_table[(int)str->s[i]] << 4*(1-i%2);
		} else c->l_qseq = 0;
		if (ks_getuntil(ks, KS_SEP_TAB, str, &dret) < 0) return -6; // qual
		z += str->l + 1;
		if (strcmp(str->s, "*") && c->l_qseq != strlen(str->s))
			parse_error(fp->n_lines, "sequence and quality are inconsistent");
		p += (c->l_qseq+1)/2;
		if (strcmp(str->s, "*") == 0) for (i = 0; i < c->l_qseq; ++i) p[i] = 0xff;
		else for (i = 0; i < c->l_qseq; ++i) p[i] = str->s[i] - 33;
		doff += c->l_qseq + (c->l_qseq+1)/2;
	}
	doff0 = doff;
	if (dret != '\n' && dret != '\r') { // aux
		while (ks_getuntil(ks, KS_SEP_TAB, str, &dret) >= 0) {
			uint8_t *s, type, key[2];
			z += str->l + 1;
			if (str->l < 6 || str->s[2] != ':' || str->s[4] != ':')
				parse_error(fp->n_lines, "missing colon in auxiliary data");
			key[0] = str->s[0]; key[1] = str->s[1];
			type = str->s[3];
			s = alloc_data(b, doff + 3) + doff;
			s[0] = key[0]; s[1] = key[1]; s += 2; doff += 2;
			if (type == 'A' || type == 'a' || type == 'c' || type == 'C') { // c and C for backward compatibility
				s = alloc_data(b, doff + 2) + doff;
				*s++ = 'A'; *s = str->s[5];
				doff += 2;
			} else if (type == 'I' || type == 'i') {
				long long x;
				s = alloc_data(b, doff + 5) + doff;
				x = (long long)atoll(str->s + 5);
				if (x < 0) {
					if (x >= -127) {
						*s++ = 'c'; *(int8_t*)s = (int8_t)x;
						s += 1; doff += 2;
					} else if (x >= -32767) {
						*s++ = 's'; *(int16_t*)s = (int16_t)x;
						s += 2; doff += 3;
					} else {
						*s++ = 'i'; *(int32_t*)s = (int32_t)x;
						s += 4; doff += 5;
						if (x < -2147483648ll)
							fprintf(stderr, "Parse warning at line %lld: integer %lld is out of range.",
									(long long)fp->n_lines, x);
					}
				} else {
					if (x <= 255) {
						*s++ = 'C'; *s++ = (uint8_t)x;
						doff += 2;
					} else if (x <= 65535) {
						*s++ = 'S'; *(uint16_t*)s = (uint16_t)x;
						s += 2; doff += 3;
					} else {
						*s++ = 'I'; *(uint32_t*)s = (uint32_t)x;
						s += 4; doff += 5;
						if (x > 4294967295ll)
							fprintf(stderr, "Parse warning at line %lld: integer %lld is out of range.",
									(long long)fp->n_lines, x);
					}
				}
			} else if (type == 'f') {
				s = alloc_data(b, doff + 5) + doff;
				*s++ = 'f';
				*(float*)s = (float)atof(str->s + 5);
				s += 4; doff += 5;
			} else if (type == 'd') {
				s = alloc_data(b, doff + 9) + doff;
				*s++ = 'd';
				*(float*)s = (float)atof(str->s + 9);
				s += 8; doff += 9;
			} else if (type == 'Z' || type == 'H') {
				int size = 1 + (str->l - 5) + 1;
				if (type == 'H') { // check whether the hex string is valid
					int i;
					if ((str->l - 5) % 2 == 1) parse_error(fp->n_lines, "length of the hex string not even");
					for (i = 0; i < str->l - 5; ++i) {
						int c = toupper(str->s[5 + i]);
						if (!((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F')))
							parse_error(fp->n_lines, "invalid hex character");
					}
				}
				s = alloc_data(b, doff + size) + doff;
				*s++ = type;
				memcpy(s, str->s + 5, str->l - 5);
				s[str->l - 5] = 0;
				doff += size;
			} else if (type == 'B') {
				int32_t n = 0, Bsize, k = 0, size;
				char *p;
				if (str->l < 8) parse_error(fp->n_lines, "too few values in aux type B");
				Bsize = bam_aux_type2size(str->s[5]); // the size of each element
				for (p = (char*)str->s + 6; *p; ++p) // count the number of elements in the array
					if (*p == ',') ++n;
				p = str->s + 7; // now p points to the first number in the array
				size = 6 + Bsize * n; // total number of bytes allocated to this tag
				s = alloc_data(b, doff + 6 * Bsize * n) + doff; // allocate memory
				*s++ = 'B'; *s++ = str->s[5];
				memcpy(s, &n, 4); s += 4; // write the number of elements
				if (str->s[5] == 'c')      while (p < str->s + str->l) ((int8_t*)s)[k++]   = (int8_t)strtol(p, &p, 0),   ++p;
				else if (str->s[5] == 'C') while (p < str->s + str->l) ((uint8_t*)s)[k++]  = (uint8_t)strtol(p, &p, 0),  ++p;
				else if (str->s[5] == 's') while (p < str->s + str->l) ((int16_t*)s)[k++]  = (int16_t)strtol(p, &p, 0),  ++p; // FIXME: avoid unaligned memory
				else if (str->s[5] == 'S') while (p < str->s + str->l) ((uint16_t*)s)[k++] = (uint16_t)strtol(p, &p, 0), ++p;
				else if (str->s[5] == 'i') while (p < str->s + str->l) ((int32_t*)s)[k++]  = (int32_t)strtol(p, &p, 0),  ++p;
				else if (str->s[5] == 'I') while (p < str->s + str->l) ((uint32_t*)s)[k++] = (uint32_t)strtol(p, &p, 0), ++p;
				else if (str->s[5] == 'f') while (p < str->s + str->l) ((float*)s)[k++]    = (float)strtod(p, &p),       ++p;
				else parse_error(fp->n_lines, "unrecognized array type");
				s += Bsize * n; doff += size;
			} else parse_error(fp->n_lines, "unrecognized type");
			if (dret == '\n' || dret == '\r') break;
		}
	}
	b->l_aux = doff - doff0;
	b->data_len = doff;
	if (bam_no_B) bam_remove_B(b);
	return z;
}

tamFile sam_open(const char *fn)
{
	tamFile fp;
	gzFile gzfp = (strcmp(fn, "-") == 0)? gzdopen(fileno(stdin), "rb") : gzopen(fn, "rb");
	if (gzfp == 0) return 0;
	fp = (tamFile)calloc(1, sizeof(struct __tamFile_t));
	fp->str = (kstring_t*)calloc(1, sizeof(kstring_t));
	fp->fp = gzfp;
	fp->ks = ks_init(fp->fp);
	return fp;
}

void sam_close(tamFile fp)
{
	if (fp) {
		ks_destroy(fp->ks);
		gzclose(fp->fp);
		free(fp->str->s); free(fp->str);
		free(fp);
	}
}
/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include <sys/types.h>
//#include "bgzf.h"

#ifdef _USE_KNETFILE
//#include "knetfile.h"
typedef knetFile *_bgzf_file_t;
#define _bgzf_open(fn, mode) knet_open(fn, mode)
#define _bgzf_dopen(fp, mode) knet_dopen(fp, mode)
#define _bgzf_close(fp) knet_close(fp)
#define _bgzf_fileno(fp) ((fp)->fd)
#define _bgzf_tell(fp) knet_tell(fp)
#define _bgzf_seek(fp, offset, whence) knet_seek(fp, offset, whence)
#define _bgzf_read(fp, buf, len) knet_read(fp, buf, len)
#define _bgzf_write(fp, buf, len) knet_write(fp, buf, len)
#else // ~defined(_USE_KNETFILE)
#if defined(_WIN32) || defined(_MSC_VER)
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else // ~defined(_WIN32)
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif // ~defined(_WIN32)
typedef FILE *_bgzf_file_t;
#define _bgzf_open(fn, mode) fopen(fn, mode)
#define _bgzf_dopen(fp, mode) fdopen(fp, mode)
#define _bgzf_close(fp) fclose(fp)
#define _bgzf_fileno(fp) fileno(fp)
#define _bgzf_tell(fp) ftello(fp)
#define _bgzf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define _bgzf_read(fp, buf, len) fread(buf, 1, len, fp)
#define _bgzf_write(fp, buf, len) fwrite(buf, 1, len, fp)
#endif // ~define(_USE_KNETFILE)

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8


/* BGZF/GZIP header (speciallized from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
*/
static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";

#ifdef BGZF_CACHE
typedef struct {
	int size;
	uint8_t *block;
	int64_t end_offset;
} cache_t;
//#include "khash.h"
KHASH_MAP_INIT_INT64(cache, cache_t)
#endif

static inline void packInt16(uint8_t *buffer, uint16_t value)
{
	buffer[0] = value;
	buffer[1] = value >> 8;
}

static inline int unpackInt16(const uint8_t *buffer)
{
	return buffer[0] | buffer[1] << 8;
}

static inline void packInt32(uint8_t *buffer, uint32_t value)
{
	buffer[0] = value;
	buffer[1] = value >> 8;
	buffer[2] = value >> 16;
	buffer[3] = value >> 24;
}

static BGZF *bgzf_read_init()
{
	BGZF *fp;
	fp = calloc(1, sizeof(BGZF));
	fp->is_write = 0;
	fp->uncompressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
	fp->compressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
#ifdef BGZF_CACHE
	fp->cache = kh_init(cache);
#endif
	return fp;
}

static BGZF *bgzf_write_init(int compress_level) // compress_level==-1 for the default level
{
	BGZF *fp;
	fp = calloc(1, sizeof(BGZF));
	fp->is_write = 1;
	fp->uncompressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
	fp->compressed_block = malloc(BGZF_MAX_BLOCK_SIZE);
	fp->compress_level = compress_level < 0? Z_DEFAULT_COMPRESSION : compress_level; // Z_DEFAULT_COMPRESSION==-1
	if (fp->compress_level > 9) fp->compress_level = Z_DEFAULT_COMPRESSION;
	return fp;
}
// get the compress level from the mode string
static int mode2level(const char *__restrict mode)
{
	int i, compress_level = -1;
	for (i = 0; mode[i]; ++i)
		if (mode[i] >= '0' && mode[i] <= '9') break;
	if (mode[i]) compress_level = (int)mode[i] - '0';
	if (strchr(mode, 'u')) compress_level = 0;
	return compress_level;
}

BGZF *bgzf_open(const char *path, const char *mode)
{
	BGZF *fp = 0;
	assert(compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);
	if (strchr(mode, 'r') || strchr(mode, 'R')) {
		_bgzf_file_t fpr;
		if ((fpr = _bgzf_open(path, "r")) == 0) return 0;
		fp = bgzf_read_init();
		fp->fp = fpr;
	} else if (strchr(mode, 'w') || strchr(mode, 'W')) {
		FILE *fpw;
		if ((fpw = fopen(path, "w")) == 0) return 0;
		fp = bgzf_write_init(mode2level(mode));
		fp->fp = fpw;
	}
	return fp;
}

BGZF *bgzf_dopen(int fd, const char *mode)
{
	BGZF *fp = 0;
	assert(compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE);
	if (strchr(mode, 'r') || strchr(mode, 'R')) {
		_bgzf_file_t fpr;
		if ((fpr = _bgzf_dopen(fd, "r")) == 0) return 0;
		fp = bgzf_read_init();
		fp->fp = fpr;
	} else if (strchr(mode, 'w') || strchr(mode, 'W')) {
		FILE *fpw;
		if ((fpw = fdopen(fd, "w")) == 0) return 0;
		fp = bgzf_write_init(mode2level(mode));
		fp->fp = fpw;
	}
	return fp;
}

static int bgzf_compress(void *_dst, int *dlen, void *src, int slen, int level)
{
	uint32_t crc;
	z_stream zs;
	uint8_t *dst = (uint8_t*)_dst;

	// compress the body
	zs.zalloc = NULL; zs.zfree = NULL;
	zs.next_in  = src;
	zs.avail_in = slen;
	zs.next_out = dst + BLOCK_HEADER_LENGTH;
	zs.avail_out = *dlen - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
	if (deflateInit2(&zs, level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) return -1; // -15 to disable zlib header/footer
	if (deflate(&zs, Z_FINISH) != Z_STREAM_END) return -1;
	if (deflateEnd(&zs) != Z_OK) return -1;
	*dlen = zs.total_out + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
	// write the header
	memcpy(dst, g_magic, BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
	packInt16(&dst[16], *dlen - 1); // write the compressed length; -1 to fit 2 bytes
	// write the footer
	crc = crc32(crc32(0L, NULL, 0L), src, slen);
	packInt32((uint8_t*)&dst[*dlen - 8], crc);
	packInt32((uint8_t*)&dst[*dlen - 4], slen);
	return 0;
}

// Deflate the block in fp->uncompressed_block into fp->compressed_block. Also adds an extra field that stores the compressed block length.
static int deflate_block(BGZF *fp, int block_length)
{
	int comp_size = BGZF_MAX_BLOCK_SIZE;
	if (bgzf_compress(fp->compressed_block, &comp_size, fp->uncompressed_block, block_length, fp->compress_level) != 0) {
		fp->errcode |= BGZF_ERR_ZLIB;
		return -1;
	}
	fp->block_offset = 0;
	return comp_size;
}

// Inflate the block in fp->compressed_block into fp->uncompressed_block
static int inflate_block(BGZF* fp, int block_length)
{
	z_stream zs;
	zs.zalloc = NULL;
	zs.zfree = NULL;
	zs.next_in = fp->compressed_block + 18;
	zs.avail_in = block_length - 16;
	zs.next_out = fp->uncompressed_block;
	zs.avail_out = BGZF_MAX_BLOCK_SIZE;

	if (inflateInit2(&zs, -15) != Z_OK) {
		fp->errcode |= BGZF_ERR_ZLIB;
		return -1;
	}
	if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
		inflateEnd(&zs);
		fp->errcode |= BGZF_ERR_ZLIB;
		return -1;
	}
	if (inflateEnd(&zs) != Z_OK) {
		fp->errcode |= BGZF_ERR_ZLIB;
		return -1;
	}
	return zs.total_out;
}

static int check_header(const uint8_t *header)
{
	return (header[0] == 31 && header[1] == 139 && header[2] == 8 && (header[3] & 4) != 0
			&& unpackInt16((uint8_t*)&header[10]) == 6
			&& header[12] == 'B' && header[13] == 'C'
			&& unpackInt16((uint8_t*)&header[14]) == 2);
}

#ifdef BGZF_CACHE
static void free_cache(BGZF *fp)
{
	khint_t k;
	khash_t(cache) *h = (khash_t(cache)*)fp->cache;
	if (fp->is_write) return;
	for (k = kh_begin(h); k < kh_end(h); ++k)
		if (kh_exist(h, k)) free(kh_val(h, k).block);
	kh_destroy(cache, h);
}

static int load_block_from_cache(BGZF *fp, int64_t block_address)
{
	khint_t k;
	cache_t *p;
	khash_t(cache) *h = (khash_t(cache)*)fp->cache;
	k = kh_get(cache, h, block_address);
	if (k == kh_end(h)) return 0;
	p = &kh_val(h, k);
	if (fp->block_length != 0) fp->block_offset = 0;
	fp->block_address = block_address;
	fp->block_length = p->size;
	memcpy(fp->uncompressed_block, p->block, BGZF_MAX_BLOCK_SIZE);
	_bgzf_seek((_bgzf_file_t)fp->fp, p->end_offset, SEEK_SET);
	return p->size;
}

static void cache_block(BGZF *fp, int size)
{
	int ret;
	khint_t k;
	cache_t *p;
	khash_t(cache) *h = (khash_t(cache)*)fp->cache;
	if (BGZF_MAX_BLOCK_SIZE >= fp->cache_size) return;
	if ((kh_size(h) + 1) * BGZF_MAX_BLOCK_SIZE > fp->cache_size) {
		/* A better way would be to remove the oldest block in the
		 * cache, but here we remove a random one for simplicity. This
		 * should not have a big impact on performance. */
		for (k = kh_begin(h); k < kh_end(h); ++k)
			if (kh_exist(h, k)) break;
		if (k < kh_end(h)) {
			free(kh_val(h, k).block);
			kh_del(cache, h, k);
		}
	}
	k = kh_put(cache, h, fp->block_address, &ret);
	if (ret == 0) return; // if this happens, a bug!
	p = &kh_val(h, k);
	p->size = fp->block_length;
	p->end_offset = fp->block_address + size;
	p->block = malloc(BGZF_MAX_BLOCK_SIZE);
	memcpy(kh_val(h, k).block, fp->uncompressed_block, BGZF_MAX_BLOCK_SIZE);
}
#else
static void free_cache(BGZF *fp) {}
static int load_block_from_cache(BGZF *fp, int64_t block_address) {return 0;}
static void cache_block(BGZF *fp, int size) {}
#endif

int bgzf_read_block(BGZF *fp)
{
	uint8_t header[BLOCK_HEADER_LENGTH], *compressed_block;
	int count, size = 0, block_length, remaining;
	int64_t block_address;
	block_address = _bgzf_tell((_bgzf_file_t)fp->fp);
	if (fp->cache_size && load_block_from_cache(fp, block_address)) return 0;
	count = _bgzf_read(fp->fp, header, sizeof(header));
	if (count == 0) { // no data read
		fp->block_length = 0;
		return 0;
	}
	if (count != sizeof(header) || !check_header(header)) {
		fp->errcode |= BGZF_ERR_HEADER;
		return -1;
	}
	size = count;
	block_length = unpackInt16((uint8_t*)&header[16]) + 1; // +1 because when writing this number, we used "-1"
	compressed_block = (uint8_t*)fp->compressed_block;
	memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
	remaining = block_length - BLOCK_HEADER_LENGTH;
	count = _bgzf_read(fp->fp, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
	if (count != remaining) {
		fp->errcode |= BGZF_ERR_IO;
		return -1;
	}
	size += count;
	if ((count = inflate_block(fp, block_length)) < 0) return -1;
	if (fp->block_length != 0) fp->block_offset = 0; // Do not reset offset if this read follows a seek.
	fp->block_address = block_address;
	fp->block_length = count;
	cache_block(fp, size);
	return 0;
}

ssize_t bgzf_read(BGZF *fp, void *data, ssize_t length)
{
	ssize_t bytes_read = 0;
	uint8_t *output = data;
	if (length <= 0) return 0;
	assert(fp->is_write == 0);
	while (bytes_read < length) {
		int copy_length, available = fp->block_length - fp->block_offset;
		uint8_t *buffer;
		if (available <= 0) {
			if (bgzf_read_block(fp) != 0) return -1;
			available = fp->block_length - fp->block_offset;
			if (available <= 0) break;
		}
		copy_length = length - bytes_read < available? length - bytes_read : available;
		buffer = fp->uncompressed_block;
		memcpy(output, buffer + fp->block_offset, copy_length);
		fp->block_offset += copy_length;
		output += copy_length;
		bytes_read += copy_length;
	}
	if (fp->block_offset == fp->block_length) {
		fp->block_address = _bgzf_tell((_bgzf_file_t)fp->fp);
		fp->block_offset = fp->block_length = 0;
	}
	return bytes_read;
}

/***** BEGIN: multi-threading *****/

typedef struct {
	BGZF *fp;
	struct mtaux_t *mt;
	void *buf;
	int i, errcode, toproc;
} worker_t;

typedef struct mtaux_t {
	int n_threads, n_blks, curr, done;
	volatile int proc_cnt;
	void **blk;
	int *len;
	worker_t *w;
	pthread_t *tid;
	pthread_mutex_t lock;
	pthread_cond_t cv;
} mtaux_t;

static int worker_aux(worker_t *w)
{
	int i, tmp, stop = 0;
	// wait for condition: to process or all done
	pthread_mutex_lock(&w->mt->lock);
	while (!w->toproc && !w->mt->done)
		pthread_cond_wait(&w->mt->cv, &w->mt->lock);
	if (w->mt->done) stop = 1;
	w->toproc = 0;
	pthread_mutex_unlock(&w->mt->lock);
	if (stop) return 1; // to quit the thread
	w->errcode = 0;
	for (i = w->i; i < w->mt->curr; i += w->mt->n_threads) {
		int clen = BGZF_MAX_BLOCK_SIZE;
		if (bgzf_compress(w->buf, &clen, w->mt->blk[i], w->mt->len[i], w->fp->compress_level) != 0)
			w->errcode |= BGZF_ERR_ZLIB;
		memcpy(w->mt->blk[i], w->buf, clen);
		w->mt->len[i] = clen;
	}
	tmp = __sync_fetch_and_add(&w->mt->proc_cnt, 1);
	return 0;
}

static void *mt_worker(void *data)
{
	while (worker_aux(data) == 0);
	return 0;
}

int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks)
{
	int i;
	mtaux_t *mt;
	pthread_attr_t attr;
	if (!fp->is_write || fp->mt || n_threads <= 1) return -1;
	mt = calloc(1, sizeof(mtaux_t));
	mt->n_threads = n_threads;
	mt->n_blks = n_threads * n_sub_blks;
	mt->len = calloc(mt->n_blks, sizeof(int));
	mt->blk = calloc(mt->n_blks, sizeof(void*));
	for (i = 0; i < mt->n_blks; ++i)
		mt->blk[i] = malloc(BGZF_MAX_BLOCK_SIZE);
	mt->tid = calloc(mt->n_threads, sizeof(pthread_t)); // tid[0] is not used, as the worker 0 is launched by the master
	mt->w = calloc(mt->n_threads, sizeof(worker_t));
	for (i = 0; i < mt->n_threads; ++i) {
		mt->w[i].i = i;
		mt->w[i].mt = mt;
		mt->w[i].fp = fp;
		mt->w[i].buf = malloc(BGZF_MAX_BLOCK_SIZE);
	}
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_mutex_init(&mt->lock, 0);
	pthread_cond_init(&mt->cv, 0);
	for (i = 1; i < mt->n_threads; ++i) // worker 0 is effectively launched by the master thread
		pthread_create(&mt->tid[i], &attr, mt_worker, &mt->w[i]);
	fp->mt = mt;
	return 0;
}

static void mt_destroy(mtaux_t *mt)
{
	int i;
	// signal all workers to quit
	pthread_mutex_lock(&mt->lock);
	mt->done = 1; mt->proc_cnt = 0;
	pthread_cond_broadcast(&mt->cv);
	pthread_mutex_unlock(&mt->lock);
	for (i = 1; i < mt->n_threads; ++i) pthread_join(mt->tid[i], 0); // worker 0 is effectively launched by the master thread
	// free other data allocated on heap
	for (i = 0; i < mt->n_blks; ++i) free(mt->blk[i]);
	for (i = 0; i < mt->n_threads; ++i) free(mt->w[i].buf);
	free(mt->blk); free(mt->len); free(mt->w); free(mt->tid);
	pthread_cond_destroy(&mt->cv);
	pthread_mutex_destroy(&mt->lock);
	free(mt);
}

static void mt_queue(BGZF *fp)
{
	mtaux_t *mt = (mtaux_t*)fp->mt;
	assert(mt->curr < mt->n_blks); // guaranteed by the caller
	memcpy(mt->blk[mt->curr], fp->uncompressed_block, fp->block_offset);
	mt->len[mt->curr] = fp->block_offset;
	fp->block_offset = 0;
	++mt->curr;
}

static int mt_flush(BGZF *fp)
{
	int i;
	mtaux_t *mt = (mtaux_t*)fp->mt;
	if (fp->block_offset) mt_queue(fp); // guaranteed that assertion does not fail
	// signal all the workers to compress
	pthread_mutex_lock(&mt->lock);
	for (i = 0; i < mt->n_threads; ++i) mt->w[i].toproc = 1;
	mt->proc_cnt = 0;
	pthread_cond_broadcast(&mt->cv);
	pthread_mutex_unlock(&mt->lock);
	// worker 0 is doing things here
	worker_aux(&mt->w[0]);
	// wait for all the threads to complete
	while (mt->proc_cnt < mt->n_threads);
	// dump data to disk
	for (i = 0; i < mt->n_threads; ++i) fp->errcode |= mt->w[i].errcode;
	for (i = 0; i < mt->curr; ++i)
		if (fwrite(mt->blk[i], 1, mt->len[i], fp->fp) != mt->len[i])
			fp->errcode |= BGZF_ERR_IO;
	mt->curr = 0;
	return 0;
}

static int mt_lazy_flush(BGZF *fp)
{
	mtaux_t *mt = (mtaux_t*)fp->mt;
	if (fp->block_offset) mt_queue(fp);
	if (mt->curr == mt->n_blks)
		return mt_flush(fp);
	return -1;
}

static ssize_t mt_write(BGZF *fp, const void *data, ssize_t length)
{
	const uint8_t *input = data;
	ssize_t rest = length;
	while (rest) {
		int copy_length = BGZF_BLOCK_SIZE - fp->block_offset < rest? BGZF_BLOCK_SIZE - fp->block_offset : rest;
		memcpy(fp->uncompressed_block + fp->block_offset, input, copy_length);
		fp->block_offset += copy_length; input += copy_length; rest -= copy_length;
		if (fp->block_offset == BGZF_BLOCK_SIZE) mt_lazy_flush(fp);
	}
	return length - rest;
}

/***** END: multi-threading *****/

int bgzf_flush(BGZF *fp)
{
	if (!fp->is_write) return 0;
	if (fp->mt) return mt_flush(fp);
	while (fp->block_offset > 0) {
		int block_length;
		block_length = deflate_block(fp, fp->block_offset);
		if (block_length < 0) return -1;
		if (fwrite(fp->compressed_block, 1, block_length, fp->fp) != block_length) {
			fp->errcode |= BGZF_ERR_IO; // possibly truncated file
			return -1;
		}
		fp->block_address += block_length;
	}
	return 0;
}

int bgzf_flush_try(BGZF *fp, ssize_t size)
{
	if (fp->block_offset + size > BGZF_BLOCK_SIZE) {
		if (fp->mt) return mt_lazy_flush(fp);
		else return bgzf_flush(fp);
	}
	return -1;
}

ssize_t bgzf_write(BGZF *fp, const void *data, ssize_t length)
{
	const uint8_t *input = data;
	int block_length = BGZF_BLOCK_SIZE, bytes_written = 0;
	assert(fp->is_write);
	if (fp->mt) return mt_write(fp, data, length);
	while (bytes_written < length) {
		uint8_t* buffer = fp->uncompressed_block;
		int copy_length = block_length - fp->block_offset < length - bytes_written? block_length - fp->block_offset : length - bytes_written;
		memcpy(buffer + fp->block_offset, input, copy_length);
		fp->block_offset += copy_length;
		input += copy_length;
		bytes_written += copy_length;
		if (fp->block_offset == block_length && bgzf_flush(fp)) break;
	}
	return bytes_written;
}

int bgzf_close(BGZF* fp)
{
	int ret, count, block_length;
	if (fp == 0) return -1;
	if (fp->is_write) {
		if (bgzf_flush(fp) != 0) return -1;
		fp->compress_level = -1;
		block_length = deflate_block(fp, 0); // write an empty block
		count = fwrite(fp->compressed_block, 1, block_length, fp->fp);
		if (fflush(fp->fp) != 0) {
			fp->errcode |= BGZF_ERR_IO;
			return -1;
		}
		if (fp->mt) mt_destroy(fp->mt);
	}
	ret = fp->is_write? fclose(fp->fp) : _bgzf_close(fp->fp);
	if (ret != 0) return -1;
	free(fp->uncompressed_block);
	free(fp->compressed_block);
	free_cache(fp);
	free(fp);
	return 0;
}

void bgzf_set_cache_size(BGZF *fp, int cache_size)
{
	if (fp) fp->cache_size = cache_size;
}

int bgzf_check_EOF(BGZF *fp)
{
	static uint8_t magic[28] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
	uint8_t buf[28];
	off_t offset;
	offset = _bgzf_tell((_bgzf_file_t)fp->fp);
	if (_bgzf_seek(fp->fp, -28, SEEK_END) < 0) return 0;
	_bgzf_read(fp->fp, buf, 28);
	_bgzf_seek(fp->fp, offset, SEEK_SET);
	return (memcmp(magic, buf, 28) == 0)? 1 : 0;
}

int64_t bgzf_seek(BGZF* fp, int64_t pos, int where)
{
	int block_offset;
	int64_t block_address;

	if (fp->is_write || where != SEEK_SET) {
		fp->errcode |= BGZF_ERR_MISUSE;
		return -1;
	}
	block_offset = pos & 0xFFFF;
	block_address = pos >> 16;
	if (_bgzf_seek(fp->fp, block_address, SEEK_SET) < 0) {
		fp->errcode |= BGZF_ERR_IO;
		return -1;
	}
	fp->block_length = 0;  // indicates current block has not been loaded
	fp->block_address = block_address;
	fp->block_offset = block_offset;
	return 0;
}

int bgzf_is_bgzf(const char *fn)
{
	uint8_t buf[16];
	int n;
	_bgzf_file_t fp;
	if ((fp = _bgzf_open(fn, "r")) == 0) return 0;
	n = _bgzf_read(fp, buf, 16);
	_bgzf_close(fp);
	if (n != 16) return 0;
	return memcmp(g_magic, buf, 16) == 0? 1 : 0;
}

int bgzf_getc(BGZF *fp)
{
	int c;
	if (fp->block_offset >= fp->block_length) {
		if (bgzf_read_block(fp) != 0) return -2; /* error */
		if (fp->block_length == 0) return -1; /* end-of-file */
	}
	c = ((unsigned char*)fp->uncompressed_block)[fp->block_offset++];
    if (fp->block_offset == fp->block_length) {
        fp->block_address = _bgzf_tell((_bgzf_file_t)fp->fp);
        fp->block_offset = 0;
        fp->block_length = 0;
    }
	return c;
}

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

int bgzf_getline(BGZF *fp, int delim, kstring_t *str)
{
	int l, state = 0;
	unsigned char *buf = (unsigned char*)fp->uncompressed_block;
	str->l = 0;
	do {
		if (fp->block_offset >= fp->block_length) {
			if (bgzf_read_block(fp) != 0) { state = -2; break; }
			if (fp->block_length == 0) { state = -1; break; }
		}
		for (l = fp->block_offset; l < fp->block_length && buf[l] != delim; ++l);
		if (l < fp->block_length) state = 1;
		l -= fp->block_offset;
		if (str->l + l + 1 >= str->m) {
			str->m = str->l + l + 2;
			kroundup32(str->m);
			str->s = (char*)realloc(str->s, str->m);
		}
		memcpy(str->s + str->l, buf + fp->block_offset, l);
		str->l += l;
		fp->block_offset += l + 1;
		if (fp->block_offset >= fp->block_length) {
			fp->block_address = _bgzf_tell((_bgzf_file_t)fp->fp);
			fp->block_offset = 0;
			fp->block_length = 0;
		} 
	} while (state == 0);
	if (str->l == 0 && state < 0) return state;
	str->s[str->l] = 0;
	return str->l;
}
//#include "sam_header.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdarg.h>

//#include "khash.h"
KHASH_MAP_INIT_STR(str, const char *)

struct _HeaderList
{
    struct _HeaderList *last;   // Hack: Used and maintained only by list_append_to_end. Maintained in the root node only.
    struct _HeaderList *next;
    void *data;
};
typedef struct _HeaderList list_t;
typedef list_t HeaderDict;

typedef struct
{
    char key[2];
    char *value;
}
HeaderTag;

typedef struct
{
    char type[2];
    list_t *tags;
}
HeaderLine;

const char *o_hd_tags[] = {"SO","GO",NULL};
const char *r_hd_tags[] = {"VN",NULL};

const char *o_sq_tags[] = {"AS","M5","UR","SP",NULL};
const char *r_sq_tags[] = {"SN","LN",NULL};
const char *u_sq_tags[] = {"SN",NULL};

const char *o_rg_tags[] = {"CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM",NULL};
const char *r_rg_tags[] = {"ID",NULL};
const char *u_rg_tags[] = {"ID",NULL};

const char *o_pg_tags[] = {"VN","CL",NULL};
const char *r_pg_tags[] = {"ID",NULL};

const char *types[]          = {"HD","SQ","RG","PG","CO",NULL};
const char **optional_tags[] = {o_hd_tags,o_sq_tags,o_rg_tags,o_pg_tags,NULL,NULL};
const char **required_tags[] = {r_hd_tags,r_sq_tags,r_rg_tags,r_pg_tags,NULL,NULL};
const char **unique_tags[]   = {NULL,     u_sq_tags,u_rg_tags,NULL,NULL,NULL};


static void debug(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}

#if 0
// Replaced by list_append_to_end
static list_t *list_prepend(list_t *root, void *data)
{
    list_t *l = malloc(sizeof(list_t));
    l->next = root;
    l->data = data;
    return l;
}
#endif

// Relies on the root->last being correct. Do not use with the other list_*
//  routines unless they are fixed to modify root->last as well.
static list_t *list_append_to_end(list_t *root, void *data)
{
    list_t *l = malloc(sizeof(list_t));
    l->last = l;
    l->next = NULL;
    l->data = data;

    if ( !root )
        return l;

    root->last->next = l;
    root->last = l;
    return root;
}

static list_t *list_append(list_t *root, void *data)
{
    list_t *l = root;
    while (l && l->next)
        l = l->next;
    if ( l ) 
    {
        l->next = malloc(sizeof(list_t));
        l = l->next;
    }
    else
    {
        l = malloc(sizeof(list_t));
        root = l;
    }
    l->data = data;
    l->next = NULL;
    return root;
}

static void list_free(list_t *root)
{
    list_t *l = root;
    while (root)
    {
        l = root;
        root = root->next;
        free(l);
    }
}



// Look for a tag "XY" in a predefined const char *[] array.
static int tag_exists(const char *tag, const char **tags)
{
    int itag=0;
    if ( !tags ) return -1;
    while ( tags[itag] )
    {
        if ( tags[itag][0]==tag[0] && tags[itag][1]==tag[1] ) return itag; 
        itag++;
    }
    return -1;
}



// Mimics the behaviour of getline, except it returns pointer to the next chunk of the text
//  or NULL if everything has been read. The lineptr should be freed by the caller. The
//  newline character is stripped.
static const char *nextline(char **lineptr, size_t *n, const char *text)
{
    int len;
    const char *to = text;

    if ( !*to ) return NULL;

    while ( *to && *to!='\n' && *to!='\r' ) to++;
    len = to - text + 1;

    if ( *to )
    {
        // Advance the pointer for the next call
        if ( *to=='\n' ) to++;
        else if ( *to=='\r' && *(to+1)=='\n' ) to+=2;
    }
    if ( !len )
        return to;

    if ( !*lineptr ) 
    {
        *lineptr = malloc(len);
        *n = len;
    }
    else if ( *n<len ) 
    {
        *lineptr = realloc(*lineptr, len);
        *n = len;
    }
    if ( !*lineptr ) {
		debug("[nextline] Insufficient memory!\n");
		return 0;
	}

    memcpy(*lineptr,text,len);
    (*lineptr)[len-1] = 0;

    return to;
}

// name points to "XY", value_from points to the first character of the value string and
//  value_to points to the last character of the value string.
static HeaderTag *new_tag(const char *name, const char *value_from, const char *value_to)
{
    HeaderTag *tag = malloc(sizeof(HeaderTag));
    int len = value_to-value_from+1;

    tag->key[0] = name[0];
    tag->key[1] = name[1];
    tag->value = malloc(len+1);
    memcpy(tag->value,value_from,len+1);
    tag->value[len] = 0;
    return tag;
}

static HeaderTag *header_line_has_tag(HeaderLine *hline, const char *key)
{
    list_t *tags = hline->tags;
    while (tags)
    {
        HeaderTag *tag = tags->data;
        if ( tag->key[0]==key[0] && tag->key[1]==key[1] ) return tag;
        tags = tags->next;
    }
    return NULL;
}


// Return codes:
//   0 .. different types or unique tags differ or conflicting tags, cannot be merged
//   1 .. all tags identical -> no need to merge, drop one
//   2 .. the unique tags match and there are some conflicting tags (same tag, different value) -> error, cannot be merged nor duplicated
//   3 .. there are some missing complementary tags and no unique conflict -> can be merged into a single line
static int sam_header_compare_lines(HeaderLine *hline1, HeaderLine *hline2)
{
    HeaderTag *t1, *t2;

    if ( hline1->type[0]!=hline2->type[0] || hline1->type[1]!=hline2->type[1] )
        return 0;

    int itype = tag_exists(hline1->type,types);
    if ( itype==-1 ) {
		debug("[sam_header_compare_lines] Unknown type [%c%c]\n", hline1->type[0],hline1->type[1]);
		return -1; // FIXME (lh3): error; I do not know how this will be handled in Petr's code
	}

    if ( unique_tags[itype] )
    {
        t1 = header_line_has_tag(hline1,unique_tags[itype][0]);
        t2 = header_line_has_tag(hline2,unique_tags[itype][0]);
        if ( !t1 || !t2 ) // this should never happen, the unique tags are required
            return 2;

        if ( strcmp(t1->value,t2->value) )
            return 0;   // the unique tags differ, cannot be merged
    }
    if ( !required_tags[itype] && !optional_tags[itype] )
    {
        t1 = hline1->tags->data;
        t2 = hline2->tags->data;
        if ( !strcmp(t1->value,t2->value) ) return 1; // identical comments
        return 0;
    }

    int missing=0, itag=0;
    while ( required_tags[itype] && required_tags[itype][itag] )
    {
        t1 = header_line_has_tag(hline1,required_tags[itype][itag]);
        t2 = header_line_has_tag(hline2,required_tags[itype][itag]);
        if ( !t1 && !t2 )
            return 2;       // this should never happen
        else if ( !t1 || !t2 )
            missing = 1;    // there is some tag missing in one of the hlines
        else if ( strcmp(t1->value,t2->value) )
        {
            if ( unique_tags[itype] )
                return 2;   // the lines have a matching unique tag but have a conflicting tag
                    
            return 0;    // the lines contain conflicting tags, cannot be merged
        }
        itag++;
    }
    itag = 0;
    while ( optional_tags[itype] && optional_tags[itype][itag] )
    {
        t1 = header_line_has_tag(hline1,optional_tags[itype][itag]);
        t2 = header_line_has_tag(hline2,optional_tags[itype][itag]);
        if ( !t1 && !t2 )
        {
            itag++;
            continue;
        }
        if ( !t1 || !t2 )
            missing = 1;    // there is some tag missing in one of the hlines
        else if ( strcmp(t1->value,t2->value) )
        {
            if ( unique_tags[itype] )
                return 2;   // the lines have a matching unique tag but have a conflicting tag

            return 0;   // the lines contain conflicting tags, cannot be merged
        }
        itag++;
    }
    if ( missing ) return 3;    // there are some missing complementary tags with no conflicts, can be merged
    return 1;
}


static HeaderLine *sam_header_line_clone(const HeaderLine *hline)
{
    list_t *tags;
    HeaderLine *out = malloc(sizeof(HeaderLine));
    out->type[0] = hline->type[0];
    out->type[1] = hline->type[1];
    out->tags = NULL;

    tags = hline->tags;
    while (tags)
    {
        HeaderTag *old = tags->data;

        HeaderTag *new = malloc(sizeof(HeaderTag));
        new->key[0] = old->key[0];
        new->key[1] = old->key[1];
        new->value  = strdup(old->value);
        out->tags = list_append(out->tags, new);

        tags = tags->next;
    }
    return out;
}

static int sam_header_line_merge_with(HeaderLine *out_hline, const HeaderLine *tmpl_hline)
{
    list_t *tmpl_tags;

    if ( out_hline->type[0]!=tmpl_hline->type[0] || out_hline->type[1]!=tmpl_hline->type[1] )
        return 0;
    
    tmpl_tags = tmpl_hline->tags;
    while (tmpl_tags)
    {
        HeaderTag *tmpl_tag = tmpl_tags->data;
        HeaderTag *out_tag  = header_line_has_tag(out_hline, tmpl_tag->key);
        if ( !out_tag )
        {
            HeaderTag *tag = malloc(sizeof(HeaderTag));
            tag->key[0] = tmpl_tag->key[0];
            tag->key[1] = tmpl_tag->key[1];
            tag->value  = strdup(tmpl_tag->value);
            out_hline->tags = list_append(out_hline->tags,tag);
        }
        tmpl_tags = tmpl_tags->next;
    }
    return 1;
}


static HeaderLine *sam_header_line_parse(const char *headerLine)
{
    HeaderLine *hline;
    HeaderTag *tag;
    const char *from, *to;
    from = headerLine;

    if ( *from != '@' ) {
		debug("[sam_header_line_parse] expected '@', got [%s]\n", headerLine);
		return 0;
	}
    to = ++from;

    while (*to && *to!='\t') to++;
    if ( to-from != 2 ) {
		debug("[sam_header_line_parse] expected '@XY', got [%s]\nHint: The header tags must be tab-separated.\n", headerLine);
		return 0;
	}
    
    hline = malloc(sizeof(HeaderLine));
    hline->type[0] = from[0];
    hline->type[1] = from[1];
    hline->tags = NULL;

    int itype = tag_exists(hline->type, types);
    
    from = to;
    while (*to && *to=='\t') to++;
    if ( to-from != 1 ) {
        debug("[sam_header_line_parse] multiple tabs on line [%s] (%d)\n", headerLine,(int)(to-from));
        free(hline);
		return 0;
	}
    from = to;
    while (*from)
    {
        while (*to && *to!='\t') to++;

        if ( !required_tags[itype] && !optional_tags[itype] )
        {
            // CO is a special case, it can contain anything, including tabs
            if ( *to ) { to++; continue; }
            tag = new_tag("  ",from,to-1);
        }
        else
            tag = new_tag(from,from+3,to-1);

        if ( header_line_has_tag(hline,tag->key) ) 
                debug("The tag '%c%c' present (at least) twice on line [%s]\n", tag->key[0],tag->key[1], headerLine);
        hline->tags = list_append(hline->tags, tag);

        from = to;
        while (*to && *to=='\t') to++;
        if ( *to && to-from != 1 ) {
			debug("[sam_header_line_parse] multiple tabs on line [%s] (%d)\n", headerLine,(int)(to-from));
			return 0;
		}

        from = to;
    }
    return hline;
}


// Must be of an existing type, all tags must be recognised and all required tags must be present
static int sam_header_line_validate(HeaderLine *hline)
{
    list_t *tags;
    HeaderTag *tag;
    int itype, itag;
    
    // Is the type correct?
    itype = tag_exists(hline->type, types);
    if ( itype==-1 ) 
    {
        debug("The type [%c%c] not recognised.\n", hline->type[0],hline->type[1]);
        return 0;
    }

    // Has all required tags?
    itag = 0;
    while ( required_tags[itype] && required_tags[itype][itag] )
    {
        if ( !header_line_has_tag(hline,required_tags[itype][itag]) )
        {
            debug("The tag [%c%c] required for [%c%c] not present.\n", required_tags[itype][itag][0],required_tags[itype][itag][1],
                hline->type[0],hline->type[1]);
            return 0;
        }
        itag++;
    }

    // Are all tags recognised?
    tags = hline->tags;
    while ( tags )
    {
        tag = tags->data;
        if ( !tag_exists(tag->key,required_tags[itype]) && !tag_exists(tag->key,optional_tags[itype]) )
        {
            // Lower case tags are user-defined values.
            if( !(islower(tag->key[0]) || islower(tag->key[1])) )
            {
                // Neither is lower case, but tag was not recognized.
                debug("Unknown tag [%c%c] for [%c%c].\n", tag->key[0],tag->key[1], hline->type[0],hline->type[1]);
                // return 0; // Even unknown tags are allowed - for forward compatibility with new attributes
            }
            // else - allow user defined tag
        }
        tags = tags->next;
    }

    return 1;
}


static void print_header_line(FILE *fp, HeaderLine *hline)
{
    list_t *tags = hline->tags;
    HeaderTag *tag;

    fprintf(fp, "@%c%c", hline->type[0],hline->type[1]);
    while (tags)
    {
        tag = tags->data;

        fprintf(fp, "\t");
        if ( tag->key[0]!=' ' || tag->key[1]!=' ' )
            fprintf(fp, "%c%c:", tag->key[0],tag->key[1]);
        fprintf(fp, "%s", tag->value);

        tags = tags->next;
    }
    fprintf(fp,"\n");
}


static void sam_header_line_free(HeaderLine *hline)
{
    list_t *tags = hline->tags;
    while (tags)
    {
        HeaderTag *tag = tags->data;
        free(tag->value);
        free(tag);
        tags = tags->next;
    }
    list_free(hline->tags);
    free(hline);
}

void sam_header_free(void *_header)
{
	HeaderDict *header = (HeaderDict*)_header;
    list_t *hlines = header;
    while (hlines)
    {
        sam_header_line_free(hlines->data);
        hlines = hlines->next;
    }
    list_free(header);
}

HeaderDict *sam_header_clone(const HeaderDict *dict)
{
    HeaderDict *out = NULL;
    while (dict)
    {
        HeaderLine *hline = dict->data;
        out = list_append(out, sam_header_line_clone(hline));
        dict = dict->next;
    }
    return out;
}

// Returns a newly allocated string
char *sam_header_write(const void *_header)
{
	const HeaderDict *header = (const HeaderDict*)_header;
    char *out = NULL;
    int len=0, nout=0;
    const list_t *hlines;

    // Calculate the length of the string to allocate
    hlines = header;
    while (hlines)
    {
        len += 4;   // @XY and \n

        HeaderLine *hline = hlines->data;
        list_t *tags = hline->tags;
        while (tags)
        {
            HeaderTag *tag = tags->data;
            len += strlen(tag->value) + 1;                  // \t
            if ( tag->key[0]!=' ' || tag->key[1]!=' ' )
                len += strlen(tag->value) + 3;              // XY:
            tags = tags->next;
        }
        hlines = hlines->next;
    }

    nout = 0;
    out  = malloc(len+1);
    hlines = header;
    while (hlines)
    {
        HeaderLine *hline = hlines->data;

        nout += sprintf(out+nout,"@%c%c",hline->type[0],hline->type[1]);

        list_t *tags = hline->tags;
        while (tags)
        {
            HeaderTag *tag = tags->data;
            nout += sprintf(out+nout,"\t");
            if ( tag->key[0]!=' ' || tag->key[1]!=' ' )
                nout += sprintf(out+nout,"%c%c:", tag->key[0],tag->key[1]);
            nout += sprintf(out+nout,"%s", tag->value);
            tags = tags->next;
        }
        hlines = hlines->next;
        nout += sprintf(out+nout,"\n");
    }
    out[len] = 0;
    return out;
}

void *sam_header_parse2(const char *headerText)
{
    list_t *hlines = NULL;
    HeaderLine *hline;
    const char *text;
    char *buf=NULL;
    size_t nbuf = 0;
	int tovalidate = 0;

    if ( !headerText )
		return 0;

    text = headerText;
    while ( (text=nextline(&buf, &nbuf, text)) )
    {
        hline = sam_header_line_parse(buf);
        if ( hline && (!tovalidate || sam_header_line_validate(hline)) )
            // With too many (~250,000) reference sequences the header parsing was too slow with list_append.
            hlines = list_append_to_end(hlines, hline);
        else
        {
			if (hline) sam_header_line_free(hline);
			sam_header_free(hlines);
            if ( buf ) free(buf);
            return NULL;
        }
    }
    if ( buf ) free(buf);

    return hlines;
}

void *sam_header2tbl(const void *_dict, char type[2], char key_tag[2], char value_tag[2])
{
	const HeaderDict *dict = (const HeaderDict*)_dict;
    const list_t *l   = dict;
    khash_t(str) *tbl = kh_init(str);
    khiter_t k;
    int ret;

	if (_dict == 0) return tbl; // return an empty (not null) hash table
    while (l)
    {
        HeaderLine *hline = l->data;
        if ( hline->type[0]!=type[0] || hline->type[1]!=type[1] ) 
        {
            l = l->next;
            continue;
        }
        
        HeaderTag *key, *value;
        key   = header_line_has_tag(hline,key_tag);
        value = header_line_has_tag(hline,value_tag); 
        if ( !key || !value )
        {
            l = l->next;
            continue;
        }
        
        k = kh_get(str, tbl, key->value);
        if ( k != kh_end(tbl) )
            debug("[sam_header_lookup_table] They key %s not unique.\n", key->value);
        k = kh_put(str, tbl, key->value, &ret);
        kh_value(tbl, k) = value->value;

        l = l->next;
    }
    return tbl;
}

char **sam_header2list(const void *_dict, char type[2], char key_tag[2], int *_n)
{
	const HeaderDict *dict = (const HeaderDict*)_dict;
    const list_t *l   = dict;
    int max, n;
	char **ret;

	ret = 0; *_n = max = n = 0;
    while (l)
    {
        HeaderLine *hline = l->data;
        if ( hline->type[0]!=type[0] || hline->type[1]!=type[1] ) 
        {
            l = l->next;
            continue;
        }
        
        HeaderTag *key;
        key   = header_line_has_tag(hline,key_tag);
        if ( !key )
        {
            l = l->next;
            continue;
        }

		if (n == max) {
			max = max? max<<1 : 4;
			ret = realloc(ret, max * sizeof(void*));
		}
		ret[n++] = key->value;

        l = l->next;
    }
	*_n = n;
    return ret;
}

void *sam_header2key_val(void *iter, const char type[2], const char key_tag[2], const char value_tag[2], const char **_key, const char **_value)
{
    list_t *l = iter;
    if ( !l ) return NULL;

    while (l)
    {
        HeaderLine *hline = l->data;
        if ( hline->type[0]!=type[0] || hline->type[1]!=type[1] )
        {
            l = l->next;
            continue;
        }

        HeaderTag *key, *value;
        key   = header_line_has_tag(hline,key_tag);
        value = header_line_has_tag(hline,value_tag);
        if ( !key && !value ) 
        {
            l = l->next;
            continue;
        }

        *_key = key->value;
        *_value = value->value;
        return l->next;
    }
    return l;
}

const char *sam_tbl_get(void *h, const char *key)
{
	khash_t(str) *tbl = (khash_t(str)*)h;
	khint_t k;
	k = kh_get(str, tbl, key);
	return k == kh_end(tbl)? 0 : kh_val(tbl, k);
}

int sam_tbl_size(void *h)
{
	khash_t(str) *tbl = (khash_t(str)*)h;
	return h? kh_size(tbl) : 0;
}

void sam_tbl_destroy(void *h)
{
	khash_t(str) *tbl = (khash_t(str)*)h;
	kh_destroy(str, tbl);
}

void *sam_header_merge(int n, const void **_dicts)
{
	const HeaderDict **dicts = (const HeaderDict**)_dicts;
    HeaderDict *out_dict;
    int idict, status;

    if ( n<2 ) return NULL;

    out_dict = sam_header_clone(dicts[0]);

    for (idict=1; idict<n; idict++)
    {
        const list_t *tmpl_hlines = dicts[idict];

        while ( tmpl_hlines )
        {
            list_t *out_hlines = out_dict;
            int inserted = 0;
            while ( out_hlines )
            {
                status = sam_header_compare_lines(tmpl_hlines->data, out_hlines->data);
                if ( status==0 )
                {
                    out_hlines = out_hlines->next;
                    continue;
                }
                
                if ( status==2 ) 
                {
                    print_header_line(stderr,tmpl_hlines->data);
                    print_header_line(stderr,out_hlines->data);
                    debug("Conflicting lines, cannot merge the headers.\n");
					return 0;
                }
                if ( status==3 )
                    sam_header_line_merge_with(out_hlines->data, tmpl_hlines->data);

                inserted = 1;
                break;
            }
            if ( !inserted )
                out_dict = list_append(out_dict, sam_header_line_clone(tmpl_hlines->data));

            tmpl_hlines = tmpl_hlines->next;
        }
    }

    return out_dict;
}

char **sam_header2tbl_n(const void *dict, const char type[2], const char *tags[], int *n)
{
    int nout = 0;
    char **out = NULL;

    *n = 0;
    list_t *l = (list_t *)dict;
    if ( !l ) return NULL;

    int i, ntags = 0;
    while ( tags[ntags] ) ntags++;

    while (l)
    {
        HeaderLine *hline = l->data;
        if ( hline->type[0]!=type[0] || hline->type[1]!=type[1] )
        {
            l = l->next;
            continue;
        }
        out = (char**) realloc(out, sizeof(char*)*(nout+1)*ntags);
        for (i=0; i<ntags; i++)
        {
            HeaderTag *key = header_line_has_tag(hline, tags[i]);
            if ( !key ) 
            {
                out[nout*ntags+i] = NULL;
                continue;
            }
            out[nout*ntags+i] = key->value;
        }
        nout++;
        l = l->next;
    }
    *n = nout;
    return out;
}

