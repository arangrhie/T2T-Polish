/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>

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

/* The BGZF library was originally written by Bob Handsaker from the Broad
 * Institute. It was later improved by the SAMtools developers. */

#ifndef __BGZF_H
#define __BGZF_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <sys/types.h>

#define BGZF_BLOCK_SIZE     0xff00 // make sure compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE 0x10000

#define BGZF_ERR_ZLIB   1
#define BGZF_ERR_HEADER 2
#define BGZF_ERR_IO     4
#define BGZF_ERR_MISUSE 8

typedef struct {
	int errcode:16, is_write:2, compress_level:14;
	int cache_size;
    int block_length, block_offset;
    int64_t block_address;
    void *uncompressed_block, *compressed_block;
	void *cache; // a pointer to a hash table
	void *fp; // actual file handler; FILE* on writing; FILE* or knetFile* on reading
	void *mt; // only used for multi-threading
} BGZF;

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

	/******************
	 * Basic routines *
	 ******************/

	/**
	 * Open an existing file descriptor for reading or writing.
	 *
	 * @param fd    file descriptor
	 * @param mode  mode matching /[rwu0-9]+/: 'r' for reading, 'w' for writing and a digit specifies
	 *              the zlib compression level; if both 'r' and 'w' are present, 'w' is ignored.
	 * @return      BGZF file handler; 0 on error
	 */
	BGZF* bgzf_dopen(int fd, const char *mode);

	#define bgzf_fdopen(fd, mode) bgzf_dopen((fd), (mode)) // for backward compatibility

	/**
	 * Open the specified file for reading or writing.
	 */
	BGZF* bgzf_open(const char* path, const char *mode);

	/**
	 * Close the BGZF and free all associated resources.
	 *
	 * @param fp    BGZF file handler
	 * @return      0 on success and -1 on error
	 */
	int bgzf_close(BGZF *fp);

	/**
	 * Read up to _length_ bytes from the file storing into _data_.
	 *
	 * @param fp     BGZF file handler
	 * @param data   data array to read into
	 * @param length size of data to read
	 * @return       number of bytes actually read; 0 on end-of-file and -1 on error
	 */
	ssize_t bgzf_read(BGZF *fp, void *data, ssize_t length);

	/**
	 * Write _length_ bytes from _data_ to the file.
	 *
	 * @param fp     BGZF file handler
	 * @param data   data array to write
	 * @param length size of data to write
	 * @return       number of bytes actually written; -1 on error
	 */
	ssize_t bgzf_write(BGZF *fp, const void *data, ssize_t length);

	/**
	 * Write the data in the buffer to the file.
	 */
	int bgzf_flush(BGZF *fp);

	/**
	 * Return a virtual file pointer to the current location in the file.
	 * No interpetation of the value should be made, other than a subsequent
	 * call to bgzf_seek can be used to position the file at the same point.
	 * Return value is non-negative on success.
	 */
	#define bgzf_tell(fp) ((fp->block_address << 16) | (fp->block_offset & 0xFFFF))

	/**
	 * Set the file to read from the location specified by _pos_.
	 *
	 * @param fp     BGZF file handler
	 * @param pos    virtual file offset returned by bgzf_tell()
	 * @param whence must be SEEK_SET
	 * @return       0 on success and -1 on error
	 */
	int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence);

	/**
	 * Check if the BGZF end-of-file (EOF) marker is present
	 *
	 * @param fp    BGZF file handler opened for reading
	 * @return      1 if EOF is present; 0 if not or on I/O error
	 */
	int bgzf_check_EOF(BGZF *fp);

	/**
	 * Check if a file is in the BGZF format
	 *
	 * @param fn    file name
	 * @return      1 if _fn_ is BGZF; 0 if not or on I/O error
	 */
	 int bgzf_is_bgzf(const char *fn);

	/*********************
	 * Advanced routines *
	 *********************/

	/**
	 * Set the cache size. Only effective when compiled with -DBGZF_CACHE.
	 *
	 * @param fp    BGZF file handler
	 * @param size  size of cache in bytes; 0 to disable caching (default)
	 */
	void bgzf_set_cache_size(BGZF *fp, int size);

	/**
	 * Flush the file if the remaining buffer size is smaller than _size_ 
	 */
	int bgzf_flush_try(BGZF *fp, ssize_t size);

	/**
	 * Read one byte from a BGZF file. It is faster than bgzf_read()
	 * @param fp     BGZF file handler
	 * @return       byte read; -1 on end-of-file or error
	 */
	int bgzf_getc(BGZF *fp);

	/**
	 * Read one line from a BGZF file. It is faster than bgzf_getc()
	 *
	 * @param fp     BGZF file handler
	 * @param delim  delimitor
	 * @param str    string to write to; must be initialized
	 * @return       length of the string; 0 on end-of-file; negative on error
	 */
	int bgzf_getline(BGZF *fp, int delim, kstring_t *str);

	/**
	 * Read the next BGZF block.
	 */
	int bgzf_read_block(BGZF *fp);

	/**
	 * Enable multi-threading (only effective on writing)
	 *
	 * @param fp          BGZF file handler; must be opened for writing
	 * @param n_threads   #threads used for writing
	 * @param n_sub_blks  #blocks processed by each thread; a value 64-256 is recommended
	 */
	int bgzf_mt(BGZF *fp, int n_threads, int n_sub_blks);

#ifdef __cplusplus
}
#endif

#endif
/* The MIT License

   Copyright (c) 2008-2010 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#ifndef BAM_BAM_H
#define BAM_BAM_H

/*!
  @header

  BAM library provides I/O and various operations on manipulating files
  in the BAM (Binary Alignment/Mapping) or SAM (Sequence Alignment/Map)
  format. It now supports importing from or exporting to SAM, sorting,
  merging, generating pileup, and quickly retrieval of reads overlapped
  with a specified region.

  @copyright Genome Research Ltd.
 */

#define BAM_VERSION "0.1.19-96b5f2294a"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef BAM_LITE
#define BAM_VIRTUAL_OFFSET16
//#include "bgzf.h"
/*! @abstract BAM file handler */
typedef BGZF *bamFile;
#define bam_open(fn, mode) bgzf_open(fn, mode)
#define bam_dopen(fd, mode) bgzf_fdopen(fd, mode)
#define bam_close(fp) bgzf_close(fp)
#define bam_read(fp, buf, size) bgzf_read(fp, buf, size)
#define bam_write(fp, buf, size) bgzf_write(fp, buf, size)
#define bam_tell(fp) bgzf_tell(fp)
#define bam_seek(fp, pos, dir) bgzf_seek(fp, pos, dir)
#else
#define BAM_TRUE_OFFSET
#include <zlib.h>
typedef gzFile bamFile;
#define bam_open(fn, mode) gzopen(fn, mode)
#define bam_dopen(fd, mode) gzdopen(fd, mode)
#define bam_close(fp) gzclose(fp)
#define bam_read(fp, buf, size) gzread(fp, buf, size)
/* no bam_write/bam_tell/bam_seek() here */
#endif

/*! @typedef
  @abstract Structure for the alignment header.
  @field n_targets   number of reference sequences
  @field target_name names of the reference sequences
  @field target_len  lengths of the referene sequences
  @field dict        header dictionary
  @field hash        hash table for fast name lookup
  @field rg2lib      hash table for @RG-ID -> LB lookup
  @field l_text      length of the plain text in the header
  @field text        plain text

  @discussion Field hash points to null by default. It is a private
  member.
 */
typedef struct {
	int32_t n_targets;
	char **target_name;
	uint32_t *target_len;
	void *dict, *hash, *rg2lib;
	uint32_t l_text, n_text;
	char *text;
} bam_header_t;

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024

#define BAM_OFDEC          0
#define BAM_OFHEX          1
#define BAM_OFSTR          2

/*! @abstract defautl mask for pileup */
#define BAM_DEF_MASK (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

#define BAM_CORE_SIZE   sizeof(bam1_core_t)

/**
 * Describing how CIGAR operation/length is packed in a 32-bit integer.
 */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

/*
  CIGAR operations.
 */
/*! @abstract CIGAR: M = match or mismatch*/
#define BAM_CMATCH      0
/*! @abstract CIGAR: I = insertion to the reference */
#define BAM_CINS        1
/*! @abstract CIGAR: D = deletion from the reference */
#define BAM_CDEL        2
/*! @abstract CIGAR: N = skip on the reference (e.g. spliced alignment) */
#define BAM_CREF_SKIP   3
/*! @abstract CIGAR: S = clip on the read with clipped sequence
  present in qseq */
#define BAM_CSOFT_CLIP  4
/*! @abstract CIGAR: H = clip on the read with clipped sequence trimmed off */
#define BAM_CHARD_CLIP  5
/*! @abstract CIGAR: P = padding */
#define BAM_CPAD        6
/*! @abstract CIGAR: equals = match */
#define BAM_CEQUAL      7
/*! @abstract CIGAR: X = mismatch */
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR  "MIDNSHP=XB"
#define BAM_CIGAR_TYPE 0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference

/*! @typedef
  @abstract Structure for core alignment information.
  @field  tid     chromosome ID, defined by bam_header_t
  @field  pos     0-based leftmost coordinate
  @field  bin     bin calculated by bam_reg2bin()
  @field  qual    mapping quality
  @field  l_qname length of the query name
  @field  flag    bitwise flag
  @field  n_cigar number of CIGAR operations
  @field  l_qseq  length of the query sequence (read)
 */
typedef struct {
	int32_t tid;
	int32_t pos;
	uint32_t bin:16, qual:8, l_qname:8;
	uint32_t flag:16, n_cigar:32;
	int32_t l_qseq;
	int32_t mtid;
	int32_t mpos;
	int32_t isize;
} bam1_core_t;

/*! @typedef
  @abstract Structure for one alignment.
  @field  core       core information about the alignment
  @field  l_aux      length of auxiliary data
  @field  data_len   current length of bam1_t::data
  @field  m_data     maximum length of bam1_t::data
  @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux

  @discussion Notes:
 
   1. qname is zero tailing and core.l_qname includes the tailing '\0'.
   2. l_qseq is calculated from the total length of an alignment block
      on reading or from CIGAR.
   3. cigar data is encoded 4 bytes per CIGAR operation.
   4. seq is nybble-encoded according to bam_nt16_table.
 */
typedef struct {
	bam1_core_t core;
	int l_aux, data_len, m_data;
	uint8_t *data;
} bam1_t;

typedef struct __bam_iter_t *bam_iter_t;

#define bam1_strand(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define bam1_mstrand(b) (((b)->core.flag&BAM_FMREVERSE) != 0)

/*! @function
  @abstract  Get the CIGAR array
  @param  b  pointer to an alignment
  @return    pointer to the CIGAR array

  @discussion In the CIGAR array, each element is a 32-bit integer. The
  lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
  length of a CIGAR.
 */
#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))

/*! @function
  @abstract  Get the name of the query
  @param  b  pointer to an alignment
  @return    pointer to the name string, null terminated
 */
#define bam1_qname(b) ((char*)((b)->data))

/*! @function
  @abstract  Get query sequence
  @param  b  pointer to an alignment
  @return    pointer to sequence

  @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
  8 for T and 15 for N. Two bases are packed in one byte with the base
  at the higher 4 bits having smaller coordinate on the read. It is
  recommended to use bam1_seqi() macro to get the base.
 */
#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)

/*! @function
  @abstract  Get query quality
  @param  b  pointer to an alignment
  @return    pointer to quality string
 */
#define bam1_qual(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))

/*! @function
  @abstract  Get a base on read
  @param  s  Query sequence returned by bam1_seq()
  @param  i  The i-th position, 0-based
  @return    4-bit integer representing the base.
 */
//#define bam1_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define bam1_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

/*! @function
  @abstract  Get query sequence and quality
  @param  b  pointer to an alignment
  @return    pointer to the concatenated auxiliary data
 */
#define bam1_aux(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2)

#ifndef kroundup32
/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*!
  @abstract Whether the machine is big-endian; modified only in
  bam_header_init().
 */
extern int bam_is_be;

/*!
  @abstract Verbose level between 0 and 3; 0 is supposed to disable all
  debugging information, though this may not have been implemented.
 */
extern int bam_verbose;

extern int bam_no_B;

/*! @abstract Table for converting a nucleotide character to the 4-bit encoding. */
extern unsigned char bam_nt16_table[256];

/*! @abstract Table for converting a 4-bit encoded nucleotide to a letter. */
extern char *bam_nt16_rev_table;

extern char bam_nt16_nt4_table[];

#ifdef __cplusplus
extern "C" {
#endif

	/*********************
	 * Low-level SAM I/O *
	 *********************/

	/*! @abstract TAM file handler */
	typedef struct __tamFile_t *tamFile;

	/*!
	  @abstract   Open a SAM file for reading, either uncompressed or compressed by gzip/zlib.
	  @param  fn  SAM file name
	  @return     SAM file handler
	 */
	tamFile sam_open(const char *fn);

	/*!
	  @abstract   Close a SAM file handler
	  @param  fp  SAM file handler
	 */
	void sam_close(tamFile fp);

	/*!
	  @abstract      Read one alignment from a SAM file handler
	  @param  fp     SAM file handler
	  @param  header header information (ordered names of chromosomes)
	  @param  b      read alignment; all members in b will be updated
	  @return        0 if successful; otherwise negative
	 */
	int sam_read1(tamFile fp, bam_header_t *header, bam1_t *b);

	/*!
	  @abstract       Read header information from a TAB-delimited list file.
	  @param  fn_list file name for the list
	  @return         a pointer to the header structure

	  @discussion Each line in this file consists of chromosome name and
	  the length of chromosome.
	 */
	bam_header_t *sam_header_read2(const char *fn_list);

	/*!
	  @abstract       Read header from a SAM file (if present)
	  @param  fp      SAM file handler
	  @return         pointer to header struct; 0 if no @SQ lines available
	 */
	bam_header_t *sam_header_read(tamFile fp);

	/*!
	  @abstract       Parse @SQ lines a update a header struct
	  @param  h       pointer to the header struct to be updated
	  @return         number of target sequences

	  @discussion bam_header_t::{n_targets,target_len,target_name} will
	  be destroyed in the first place.
	 */
	int sam_header_parse(bam_header_t *h);
	int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);

	/*!
	  @abstract       Parse @RG lines a update a header struct
	  @param  h       pointer to the header struct to be updated
	  @return         number of @RG lines

	  @discussion bam_header_t::rg2lib will be destroyed in the first
	  place.
	 */
	int sam_header_parse_rg(bam_header_t *h);

#define sam_write1(header, b) bam_view1(header, b)


	/********************************
	 * APIs for string dictionaries *
	 ********************************/

	int bam_strmap_put(void *strmap, const char *rg, const char *lib);
	const char *bam_strmap_get(const void *strmap, const char *rg);
	void *bam_strmap_dup(const void*);
	void *bam_strmap_init();
	void bam_strmap_destroy(void *strmap);


	/*********************
	 * Low-level BAM I/O *
	 *********************/

	/*!
	  @abstract Initialize a header structure.
	  @return   the pointer to the header structure

	  @discussion This function also modifies the global variable
	  bam_is_be.
	 */
	bam_header_t *bam_header_init();

	/*!
	  @abstract        Destroy a header structure.
	  @param  header  pointer to the header
	 */
	void bam_header_destroy(bam_header_t *header);

	/*!
	  @abstract   Read a header structure from BAM.
	  @param  fp  BAM file handler, opened by bam_open()
	  @return     pointer to the header structure

	  @discussion The file position indicator must be placed at the
	  beginning of the file. Upon success, the position indicator will
	  be set at the start of the first alignment.
	 */
	bam_header_t *bam_header_read(bamFile fp);

	/*!
	  @abstract      Write a header structure to BAM.
	  @param  fp     BAM file handler
	  @param  header pointer to the header structure
	  @return        always 0 currently
	 */
	int bam_header_write(bamFile fp, const bam_header_t *header);

	/*!
	  @abstract   Read an alignment from BAM.
	  @param  fp  BAM file handler
	  @param  b   read alignment; all members are updated.
	  @return     number of bytes read from the file

	  @discussion The file position indicator must be
	  placed right before an alignment. Upon success, this function
	  will set the position indicator to the start of the next
	  alignment. This function is not affected by the machine
	  endianness.
	 */
	int bam_read1(bamFile fp, bam1_t *b);

	int bam_remove_B(bam1_t *b);

	/*!
	  @abstract Write an alignment to BAM.
	  @param  fp       BAM file handler
	  @param  c        pointer to the bam1_core_t structure
	  @param  data_len total length of variable size data related to
	                   the alignment
	  @param  data     pointer to the concatenated data
	  @return          number of bytes written to the file

	  @discussion This function is not affected by the machine
	  endianness.
	 */
	int bam_write1_core(bamFile fp, const bam1_core_t *c, int data_len, uint8_t *data);

	/*!
	  @abstract   Write an alignment to BAM.
	  @param  fp  BAM file handler
	  @param  b   alignment to write
	  @return     number of bytes written to the file

	  @abstract It is equivalent to:
	    bam_write1_core(fp, &b->core, b->data_len, b->data)
	 */
	int bam_write1(bamFile fp, const bam1_t *b);

	/*! @function
	  @abstract  Initiate a pointer to bam1_t struct
	 */
#if 1
#define bam_init1() ((bam1_t*)calloc(1, sizeof(bam1_t)))
#else
  bam1_t  *bam_init1(void) {
    bam1_t *b = new bam1_t;
    memset(b, 0, sizeof(bam1_t));
    return(b);
  }
#endif


	/*! @function
	  @abstract  Free the memory allocated for an alignment.
	  @param  b  pointer to an alignment
	 */
#define bam_destroy1(b) do {					\
		if (b) { free((b)->data); free(b); }	\
	} while (0)

	/*!
	  @abstract       Format a BAM record in the SAM format
	  @param  header  pointer to the header structure
	  @param  b       alignment to print
	  @return         a pointer to the SAM string
	 */
	char *bam_format1(const bam_header_t *header, const bam1_t *b);

	char *bam_format1_core(const bam_header_t *header, const bam1_t *b, int of);

	/*!
	  @abstract       Check whether a BAM record is plausibly valid
	  @param  header  associated header structure, or NULL if unavailable
	  @param  b       alignment to validate
	  @return         0 if the alignment is invalid; non-zero otherwise

	  @discussion  Simple consistency check of some of the fields of the
	  alignment record.  If the header is provided, several additional checks
	  are made.  Not all fields are checked, so a non-zero result is not a
	  guarantee that the record is valid.  However it is usually good enough
	  to detect when bam_seek() has been called with a virtual file offset
	  that is not the offset of an alignment record.
	 */
	int bam_validate1(const bam_header_t *header, const bam1_t *b);

	const char *bam_get_library(bam_header_t *header, const bam1_t *b);


	/***************
	 * pileup APIs *
	 ***************/

	/*! @typedef
	  @abstract Structure for one alignment covering the pileup position.
	  @field  b      pointer to the alignment
	  @field  qpos   position of the read base at the pileup site, 0-based
	  @field  indel  indel length; 0 for no indel, positive for ins and negative for del
	  @field  is_del 1 iff the base on the padded read is a deletion
	  @field  level  the level of the read in the "viewer" mode

	  @discussion See also bam_plbuf_push() and bam_lplbuf_push(). The
	  difference between the two functions is that the former does not
	  set bam_pileup1_t::level, while the later does. Level helps the
	  implementation of alignment viewers, but calculating this has some
	  overhead.
	 */
	typedef struct {
		bam1_t *b;
		int32_t qpos;
		int indel, level;
		uint32_t is_del:1, is_head:1, is_tail:1, is_refskip:1, aux:28;
	} bam_pileup1_t;

	typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);

	struct __bam_plp_t;
	typedef struct __bam_plp_t *bam_plp_t;

	bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data);
	int bam_plp_push(bam_plp_t iter, const bam1_t *b);
	const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
	const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
	void bam_plp_set_mask(bam_plp_t iter, int mask);
	void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt);
	void bam_plp_reset(bam_plp_t iter);
	void bam_plp_destroy(bam_plp_t iter);

	struct __bam_mplp_t;
	typedef struct __bam_mplp_t *bam_mplp_t;

	bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data);
	void bam_mplp_destroy(bam_mplp_t iter);
	void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt);
	int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp);

	/*! @typedef
	  @abstract    Type of function to be called by bam_plbuf_push().
	  @param  tid  chromosome ID as is defined in the header
	  @param  pos  start coordinate of the alignment, 0-based
	  @param  n    number of elements in pl array
	  @param  pl   array of alignments
	  @param  data user provided data
	  @discussion  See also bam_plbuf_push(), bam_plbuf_init() and bam_pileup1_t.
	 */
	typedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);

	typedef struct {
		bam_plp_t iter;
		bam_pileup_f func;
		void *data;
	} bam_plbuf_t;

	void bam_plbuf_set_mask(bam_plbuf_t *buf, int mask);
	void bam_plbuf_reset(bam_plbuf_t *buf);
	bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data);
	void bam_plbuf_destroy(bam_plbuf_t *buf);
	int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf);

	int bam_pileup_file(bamFile fp, int mask, bam_pileup_f func, void *func_data);

	struct __bam_lplbuf_t;
	typedef struct __bam_lplbuf_t bam_lplbuf_t;

	void bam_lplbuf_reset(bam_lplbuf_t *buf);

	/*! @abstract  bam_plbuf_init() equivalent with level calculated. */
	bam_lplbuf_t *bam_lplbuf_init(bam_pileup_f func, void *data);

	/*! @abstract  bam_plbuf_destroy() equivalent with level calculated. */
	void bam_lplbuf_destroy(bam_lplbuf_t *tv);

	/*! @abstract  bam_plbuf_push() equivalent with level calculated. */
	int bam_lplbuf_push(const bam1_t *b, bam_lplbuf_t *buf);


	/*********************
	 * BAM indexing APIs *
	 *********************/

	struct __bam_index_t;
	typedef struct __bam_index_t bam_index_t;

	/*!
	  @abstract   Build index for a BAM file.
	  @discussion Index file "fn.bai" will be created.
	  @param  fn  name of the BAM file
	  @return     always 0 currently
	 */
	int bam_index_build(const char *fn);

	/*!
	  @abstract   Load index from file "fn.bai".
	  @param  fn  name of the BAM file (NOT the index file)
	  @return     pointer to the index structure
	 */
	bam_index_t *bam_index_load(const char *fn);

	/*!
	  @abstract    Destroy an index structure.
	  @param  idx  pointer to the index structure
	 */
	void bam_index_destroy(bam_index_t *idx);

	/*! @typedef
	  @abstract      Type of function to be called by bam_fetch().
	  @param  b     the alignment
	  @param  data  user provided data
	 */
	typedef int (*bam_fetch_f)(const bam1_t *b, void *data);

	/*!
	  @abstract Retrieve the alignments that are overlapped with the
	  specified region.

	  @discussion A user defined function will be called for each
	  retrieved alignment ordered by its start position.

	  @param  fp    BAM file handler
	  @param  idx   pointer to the alignment index
	  @param  tid   chromosome ID as is defined in the header
	  @param  beg   start coordinate, 0-based
	  @param  end   end coordinate, 0-based
	  @param  data  user provided data (will be transferred to func)
	  @param  func  user defined function
	 */
	int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func);

	bam_iter_t bam_iter_query(const bam_index_t *idx, int tid, int beg, int end);
	int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b);
	void bam_iter_destroy(bam_iter_t iter);

	/*!
	  @abstract       Parse a region in the format: "chr2:100,000-200,000".
	  @discussion     bam_header_t::hash will be initialized if empty.
	  @param  header  pointer to the header structure
	  @param  str     string to be parsed
	  @param  ref_id  the returned chromosome ID
	  @param  begin   the returned start coordinate
	  @param  end     the returned end coordinate
	  @return         0 on success; -1 on failure
	 */
	int bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *begin, int *end);


	/**************************
	 * APIs for optional tags *
	 **************************/

	/*!
	  @abstract       Retrieve data of a tag
	  @param  b       pointer to an alignment struct
	  @param  tag     two-character tag to be retrieved

	  @return  pointer to the type and data. The first character is the
	  type that can be 'iIsScCdfAZH'.

	  @discussion  Use bam_aux2?() series to convert the returned data to
	  the corresponding type.
	*/
	uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);

	int32_t bam_aux2i(const uint8_t *s);
	float bam_aux2f(const uint8_t *s);
	double bam_aux2d(const uint8_t *s);
	char bam_aux2A(const uint8_t *s);
	char *bam_aux2Z(const uint8_t *s);

	int bam_aux_del(bam1_t *b, uint8_t *s);
	void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data);
	uint8_t *bam_aux_get_core(bam1_t *b, const char tag[2]); // an alias of bam_aux_get()


	/*****************
	 * Miscellaneous *
	 *****************/

	/*!  
	  @abstract Calculate the rightmost coordinate of an alignment on the
	  reference genome.

	  @param  c      pointer to the bam1_core_t structure
	  @param  cigar  the corresponding CIGAR array (from bam1_t::cigar)
	  @return        the rightmost coordinate, 0-based
	*/
	uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar);

	/*!
	  @abstract      Calculate the length of the query sequence from CIGAR.
	  @param  c      pointer to the bam1_core_t structure
	  @param  cigar  the corresponding CIGAR array (from bam1_t::cigar)
	  @return        length of the query sequence
	*/
	int32_t bam_cigar2qlen(const bam1_core_t *c, const uint32_t *cigar);

#ifdef __cplusplus
}
#endif

/*!
  @abstract    Calculate the minimum bin that contains a region [beg,end).
  @param  beg  start of the region, 0-based
  @param  end  end of the region, 0-based
  @return      bin
 */
static inline int bam_reg2bin(uint32_t beg, uint32_t end)
{
	--end;
	if (beg>>14 == end>>14) return 4681 + (beg>>14);
	if (beg>>17 == end>>17) return  585 + (beg>>17);
	if (beg>>20 == end>>20) return   73 + (beg>>20);
	if (beg>>23 == end>>23) return    9 + (beg>>23);
	if (beg>>26 == end>>26) return    1 + (beg>>26);
	return 0;
}

/*!
  @abstract     Copy an alignment
  @param  bdst  destination alignment struct
  @param  bsrc  source alignment struct
  @return       pointer to the destination alignment struct
 */
#if 0
static inline bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
{
	uint8_t *data = bdst->data;
	int m_data = bdst->m_data;   // backup data and m_data
	if (m_data < bsrc->data_len) { // double the capacity
		m_data = bsrc->data_len; kroundup32(m_data);
		data = (uint8_t*)realloc(data, m_data);
	}
	memcpy(data, bsrc->data, bsrc->data_len); // copy var-len data
	*bdst = *bsrc; // copy the rest
	// restore the backup
	bdst->m_data = m_data;
	bdst->data = data;
	return bdst;
}

/*!
  @abstract     Duplicate an alignment
  @param  src   source alignment struct
  @return       pointer to the destination alignment struct
 */
static inline bam1_t *bam_dup1(const bam1_t *src)
{
	bam1_t *b;
	b = bam_init1();
	*b = *src;
	b->m_data = b->data_len;
	b->data = (uint8_t*)calloc(b->data_len, 1);
	memcpy(b->data, src->data, b->data_len);
	return b;
}
#endif

static inline int bam_aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f' || x == 'F') return 4;
	else return 0;
}

/*********************************
 *** Compatibility with htslib ***
 *********************************/

typedef bam_header_t bam_hdr_t;

#define bam_get_qname(b) bam1_qname(b)
#define bam_get_cigar(b) bam1_cigar(b)

#define bam_hdr_read(fp) bam_header_read(fp)
#define bam_hdr_write(fp, h) bam_header_write(fp, h)
#define bam_hdr_destroy(fp) bam_header_destroy(fp)

#endif
#ifndef BAM_SAM_H
#define BAM_SAM_H

//#include "bam.h"

/*!
  @header

  This file provides higher level of I/O routines and unifies the APIs
  for SAM and BAM formats. These APIs are more convenient and
  recommended.

  @copyright Genome Research Ltd.
 */

/*! @typedef
  @abstract SAM/BAM file handler
  @field  type    type of the handler; bit 1 for BAM, 2 for reading and bit 3-4 for flag format
  @field  bam   BAM file handler; valid if (type&1) == 1
  @field  tamr  SAM file handler for reading; valid if type == 2
  @field  tamw  SAM file handler for writing; valid if type == 0
  @field  header  header struct
 */
typedef struct {
	int type;
	union {
		tamFile tamr;
		bamFile bam;
		FILE *tamw;
	} x;
	bam_header_t *header;
} samfile_t;

#ifdef __cplusplus
extern "C" {
#endif

	/*!
	  @abstract     Open a SAM/BAM file

	  @param fn SAM/BAM file name; "-" is recognized as stdin (for
	  reading) or stdout (for writing).

	  @param mode open mode /[rw](b?)(u?)(h?)([xX]?)/: 'r' for reading,
	  'w' for writing, 'b' for BAM I/O, 'u' for uncompressed BAM output,
	  'h' for outputing header in SAM, 'x' for HEX flag and 'X' for
	  string flag. If 'b' present, it must immediately follow 'r' or
	  'w'. Valid modes are "r", "w", "wh", "wx", "whx", "wX", "whX",
	  "rb", "wb" and "wbu" exclusively.

	  @param aux auxiliary data; if mode[0]=='w', aux points to
	  bam_header_t; if strcmp(mode, "rb")!=0 and @SQ header lines in SAM
	  are absent, aux points the file name of the list of the reference;
	  aux is not used otherwise. If @SQ header lines are present in SAM,
	  aux is not used, either.

	  @return       SAM/BAM file handler
	 */
	samfile_t *samopen(const char *fn, const char *mode, const void *aux);

	/*!
	  @abstract     Close a SAM/BAM handler
	  @param  fp    file handler to be closed
	 */
	void samclose(samfile_t *fp);

	/*!
	  @abstract     Read one alignment
	  @param  fp    file handler
	  @param  b     alignment
	  @return       bytes read
	 */
	int samread(samfile_t *fp, bam1_t *b);

	/*!
	  @abstract     Write one alignment
	  @param  fp    file handler
	  @param  b     alignment
	  @return       bytes written
	 */
	int samwrite(samfile_t *fp, const bam1_t *b);

#if 0
	/*!
	  @abstract     Get the pileup for a whole alignment file
	  @param  fp    file handler
	  @param  mask  mask transferred to bam_plbuf_set_mask()
	  @param  func  user defined function called in the pileup process
	  #param  data  user provided data for func()
	 */
	int sampileup(samfile_t *fp, int mask, bam_pileup_f func, void *data);
#endif

	char *samfaipath(const char *fn_ref);
	int samthreads(samfile_t *fp, int n_threads, int n_sub_blks);

#ifdef __cplusplus
}
#endif

#endif
#ifndef BAM_ENDIAN_H
#define BAM_ENDIAN_H

#include <stdint.h>

static inline int bam_is_big_endian()
{
	long one= 1;
	return !(*((char *)(&one)));
}
static inline uint16_t bam_swap_endian_2(uint16_t v)
{
	return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *bam_swap_endian_2p(void *x)
{
	*(uint16_t*)x = bam_swap_endian_2(*(uint16_t*)x);
	return x;
}
static inline uint32_t bam_swap_endian_4(uint32_t v)
{
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *bam_swap_endian_4p(void *x)
{
	*(uint32_t*)x = bam_swap_endian_4(*(uint32_t*)x);
	return x;
}
static inline uint64_t bam_swap_endian_8(uint64_t v)
{
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *bam_swap_endian_8p(void *x)
{
	*(uint64_t*)x = bam_swap_endian_8(*(uint64_t*)x);
	return x;
}

#endif
/* The MIT License

   Copyright (c) by Attractive Chaos <attractor@live.co.uk> 

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef KSTRING_H
#define KSTRING_H

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

typedef struct {
	uint64_t tab[4];
	int sep, finished;
	const char *p; // end of the current token
} ks_tokaux_t;

#ifdef __cplusplus
extern "C" {
#endif

	int ksprintf(kstring_t *s, const char *fmt, ...);
	int ksplit_core(char *s, int delimiter, int *_max, int **_offsets);
	char *kstrstr(const char *str, const char *pat, int **_prep);
	char *kstrnstr(const char *str, const char *pat, int n, int **_prep);
	void *kmemmem(const void *_str, int n, const void *_pat, int m, int **_prep);

	/* kstrtok() is similar to strtok_r() except that str is not
	 * modified and both str and sep can be NULL. For efficiency, it is
	 * actually recommended to set both to NULL in the subsequent calls
	 * if sep is not changed. */
	char *kstrtok(const char *str, const char *sep, ks_tokaux_t *aux);

#ifdef __cplusplus
}
#endif

static inline void ks_resize(kstring_t *s, size_t size)
{
	if (s->m < size) {
		s->m = size;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	memcpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}

static inline int kputs(const char *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

static inline int kputc(int c, kstring_t *s)
{
	if (s->l + 1 >= s->m) {
		s->m = s->l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	s->s[s->l++] = c;
	s->s[s->l] = 0;
	return c;
}

static inline int kputw(int c, kstring_t *s)
{
	char buf[16];
	int l, x;
	if (c == 0) return kputc('0', s);
        if(c < 0) for (l = 0, x = c; x < 0; x /= 10) buf[l++] = '0' - (x%10);
        else for (l = 0, x = c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputuw(unsigned c, kstring_t *s)
{
	char buf[16];
	int l, i;
	unsigned x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputl(long c, kstring_t *s)
{
	char buf[32];
	long l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

static inline int *ksplit(kstring_t *s, int delimiter, int *n)
{
	int max = 0, *offsets = 0;
	*n = ksplit_core(s->s, delimiter, &max, &offsets);
	return offsets;
}

#endif
#ifndef __SAM_HEADER_H__
#define __SAM_HEADER_H__

#ifdef __cplusplus
extern "C" {
#endif

	void *sam_header_parse2(const char *headerText);
	void *sam_header_merge(int n, const void **dicts);
	void sam_header_free(void *header);
	char *sam_header_write(const void *headerDict);   // returns a newly allocated string

    /*
        // Usage example 
        const char *key, *val; 
        void *iter = sam_header_parse2(bam->header->text);
        while ( iter = sam_header_key_val(iter, "RG","ID","SM" &key,&val) ) printf("%s\t%s\n", key,val);
    */
    void *sam_header2key_val(void *iter, const char type[2], const char key_tag[2], const char value_tag[2], const char **key, const char **value);
	char **sam_header2list(const void *_dict, char type[2], char key_tag[2], int *_n);

    /*
        // Usage example
        int i, j, n;
        const char *tags[] = {"SN","LN","UR","M5",NULL}; 
        void *dict = sam_header_parse2(bam->header->text);
        char **tbl = sam_header2tbl_n(h->dict, "SQ", tags, &n);
        for (i=0; i<n; i++)
        {
            for (j=0; j<4; j++) 
                if ( tbl[4*i+j] ) printf("\t%s", tbl[4*i+j]); 
                else printf("-");
            printf("\n");
        }
        if (tbl) free(tbl);
     */
    char **sam_header2tbl_n(const void *dict, const char type[2], const char *tags[], int *n);

	void *sam_header2tbl(const void *dict, char type[2], char key_tag[2], char value_tag[2]);
	const char *sam_tbl_get(void *h, const char *key);
	int sam_tbl_size(void *h);
	void sam_tbl_destroy(void *h);

#ifdef __cplusplus
}
#endif

#endif
/* The MIT License

   Copyright (c) 2008, 2009, 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
  An example:

#include "khash.h"
KHASH_MAP_INIT_INT(32, char)
int main() {
	int ret, is_missing;
	khiter_t k;
	khash_t(32) *h = kh_init(32);
	k = kh_put(32, h, 5, &ret);
	if (!ret) kh_del(32, h, k);
	kh_value(h, k) = 10;
	k = kh_get(32, h, 10);
	is_missing = (k == kh_end(h));
	k = kh_get(32, h, 5);
	kh_del(32, h, k);
	for (k = kh_begin(h); k != kh_end(h); ++k)
		if (kh_exist(h, k)) kh_value(h, k) = 1;
	kh_destroy(32, h);
	return 0;
}
*/

/*
  2011-02-14 (0.2.5):

    * Allow to declare global functions.

  2009-09-26 (0.2.4):

    * Improve portability

  2008-09-19 (0.2.3):

	* Corrected the example
	* Improved interfaces

  2008-09-11 (0.2.2):

	* Improved speed a little in kh_put()

  2008-09-10 (0.2.1):

	* Added kh_clear()
	* Fixed a compiling error

  2008-09-02 (0.2.0):

	* Changed to token concatenation which increases flexibility.

  2008-08-31 (0.1.2):

	* Fixed a bug in kh_get(), which has not been tested previously.

  2008-08-31 (0.1.1):

	* Added destructor
*/


#ifndef __AC_KHASH_H
#define __AC_KHASH_H

/*!
  @header

  Generic hash table library.

  @copyright Heng Li
 */

#define AC_VERSION_KHASH_H "0.2.5"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* compipler specific configuration */

#if UINT_MAX == 0xffffffffu
typedef unsigned int khint32_t;
#elif ULONG_MAX == 0xffffffffu
typedef unsigned long khint32_t;
#endif

#if ULONG_MAX == ULLONG_MAX
typedef unsigned long khint64_t;
#else
typedef unsigned long long khint64_t;
#endif

#ifdef _MSC_VER
#define inline __inline
#endif

typedef khint32_t khint_t;
typedef khint_t khiter_t;

#define __ac_HASH_PRIME_SIZE 32
static const khint32_t __ac_prime_list[__ac_HASH_PRIME_SIZE] =
{
  0ul,          3ul,          11ul,         23ul,         53ul,
  97ul,         193ul,        389ul,        769ul,        1543ul,
  3079ul,       6151ul,       12289ul,      24593ul,      49157ul,
  98317ul,      196613ul,     393241ul,     786433ul,     1572869ul,
  3145739ul,    6291469ul,    12582917ul,   25165843ul,   50331653ul,
  100663319ul,  201326611ul,  402653189ul,  805306457ul,  1610612741ul,
  3221225473ul, 4294967291ul
};

#define __ac_isempty(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&2)
#define __ac_isdel(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&1)
#define __ac_iseither(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&3)
#define __ac_set_isdel_false(flag, i) (flag[i>>4]&=~(1ul<<((i&0xfU)<<1)))
#define __ac_set_isempty_false(flag, i) (flag[i>>4]&=~(2ul<<((i&0xfU)<<1)))
#define __ac_set_isboth_false(flag, i) (flag[i>>4]&=~(3ul<<((i&0xfU)<<1)))
#define __ac_set_isdel_true(flag, i) (flag[i>>4]|=1ul<<((i&0xfU)<<1))

static const double __ac_HASH_UPPER = 0.77;

#define KHASH_DECLARE(name, khkey_t, khval_t)		 					\
	typedef struct {													\
		khint_t n_buckets, size, n_occupied, upper_bound;				\
		khint32_t *flags;												\
		khkey_t *keys;													\
		khval_t *vals;													\
	} kh_##name##_t;													\
	extern kh_##name##_t *kh_init_##name();								\
	extern void kh_destroy_##name(kh_##name##_t *h);					\
	extern void kh_clear_##name(kh_##name##_t *h);						\
	extern khint_t kh_get_##name(const kh_##name##_t *h, khkey_t key); 	\
	extern void kh_resize_##name(kh_##name##_t *h, khint_t new_n_buckets); \
	extern khint_t kh_put_##name(kh_##name##_t *h, khkey_t key, int *ret); \
	extern void kh_del_##name(kh_##name##_t *h, khint_t x);

#define KHASH_INIT2(name, SCOPE, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal) \
	typedef struct {													\
		khint_t n_buckets, size, n_occupied, upper_bound;				\
		khint32_t *flags;												\
		khkey_t *keys;													\
		khval_t *vals;													\
	} kh_##name##_t;													\
	SCOPE kh_##name##_t *kh_init_##name() {								\
		return (kh_##name##_t*)calloc(1, sizeof(kh_##name##_t));		\
	}																	\
	SCOPE void kh_destroy_##name(kh_##name##_t *h)						\
	{																	\
		if (h) {														\
			free(h->keys); free(h->flags);								\
			free(h->vals);												\
			free(h);													\
		}																\
	}																	\
	SCOPE void kh_clear_##name(kh_##name##_t *h)						\
	{																	\
		if (h && h->flags) {											\
			memset(h->flags, 0xaa, ((h->n_buckets>>4) + 1) * sizeof(khint32_t)); \
			h->size = h->n_occupied = 0;								\
		}																\
	}																	\
	SCOPE khint_t kh_get_##name(const kh_##name##_t *h, khkey_t key) 	\
	{																	\
		if (h->n_buckets) {												\
			khint_t inc, k, i, last;									\
			k = __hash_func(key); i = k % h->n_buckets;					\
			inc = 1 + k % (h->n_buckets - 1); last = i;					\
			while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
				if (i + inc >= h->n_buckets) i = i + inc - h->n_buckets; \
				else i += inc;											\
				if (i == last) return h->n_buckets;						\
			}															\
			return __ac_iseither(h->flags, i)? h->n_buckets : i;		\
		} else return 0;												\
	}																	\
	SCOPE void kh_resize_##name(kh_##name##_t *h, khint_t new_n_buckets) \
	{																	\
		khint32_t *new_flags = 0;										\
		khint_t j = 1;													\
		{																\
			khint_t t = __ac_HASH_PRIME_SIZE - 1;						\
			while (__ac_prime_list[t] > new_n_buckets) --t;				\
			new_n_buckets = __ac_prime_list[t+1];						\
			if (h->size >= (khint_t)(new_n_buckets * __ac_HASH_UPPER + 0.5)) j = 0;	\
			else {														\
				new_flags = (khint32_t*)malloc(((new_n_buckets>>4) + 1) * sizeof(khint32_t));	\
				memset(new_flags, 0xaa, ((new_n_buckets>>4) + 1) * sizeof(khint32_t)); \
				if (h->n_buckets < new_n_buckets) {						\
					h->keys = (khkey_t*)realloc(h->keys, new_n_buckets * sizeof(khkey_t)); \
					if (kh_is_map)										\
						h->vals = (khval_t*)realloc(h->vals, new_n_buckets * sizeof(khval_t)); \
				}														\
			}															\
		}																\
		if (j) {														\
			for (j = 0; j != h->n_buckets; ++j) {						\
				if (__ac_iseither(h->flags, j) == 0) {					\
					khkey_t key = h->keys[j];							\
					khval_t val;										\
					if (kh_is_map) val = h->vals[j];					\
					__ac_set_isdel_true(h->flags, j);					\
					while (1) {											\
						khint_t inc, k, i;								\
						k = __hash_func(key);							\
						i = k % new_n_buckets;							\
						inc = 1 + k % (new_n_buckets - 1);				\
						while (!__ac_isempty(new_flags, i)) {			\
							if (i + inc >= new_n_buckets) i = i + inc - new_n_buckets; \
							else i += inc;								\
						}												\
						__ac_set_isempty_false(new_flags, i);			\
						if (i < h->n_buckets && __ac_iseither(h->flags, i) == 0) { \
							{ khkey_t tmp = h->keys[i]; h->keys[i] = key; key = tmp; } \
							if (kh_is_map) { khval_t tmp = h->vals[i]; h->vals[i] = val; val = tmp; } \
							__ac_set_isdel_true(h->flags, i);			\
						} else {										\
							h->keys[i] = key;							\
							if (kh_is_map) h->vals[i] = val;			\
							break;										\
						}												\
					}													\
				}														\
			}															\
			if (h->n_buckets > new_n_buckets) {							\
				h->keys = (khkey_t*)realloc(h->keys, new_n_buckets * sizeof(khkey_t)); \
				if (kh_is_map)											\
					h->vals = (khval_t*)realloc(h->vals, new_n_buckets * sizeof(khval_t)); \
			}															\
			free(h->flags);												\
			h->flags = new_flags;										\
			h->n_buckets = new_n_buckets;								\
			h->n_occupied = h->size;									\
			h->upper_bound = (khint_t)(h->n_buckets * __ac_HASH_UPPER + 0.5); \
		}																\
	}																	\
	SCOPE khint_t kh_put_##name(kh_##name##_t *h, khkey_t key, int *ret) \
	{																	\
		khint_t x;														\
		if (h->n_occupied >= h->upper_bound) {							\
			if (h->n_buckets > (h->size<<1)) kh_resize_##name(h, h->n_buckets - 1); \
			else kh_resize_##name(h, h->n_buckets + 1);					\
		}																\
		{																\
			khint_t inc, k, i, site, last;								\
			x = site = h->n_buckets; k = __hash_func(key); i = k % h->n_buckets; \
			if (__ac_isempty(h->flags, i)) x = i;						\
			else {														\
				inc = 1 + k % (h->n_buckets - 1); last = i;				\
				while (!__ac_isempty(h->flags, i) && (__ac_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
					if (__ac_isdel(h->flags, i)) site = i;				\
					if (i + inc >= h->n_buckets) i = i + inc - h->n_buckets; \
					else i += inc;										\
					if (i == last) { x = site; break; }					\
				}														\
				if (x == h->n_buckets) {								\
					if (__ac_isempty(h->flags, i) && site != h->n_buckets) x = site; \
					else x = i;											\
				}														\
			}															\
		}																\
		if (__ac_isempty(h->flags, x)) {								\
			h->keys[x] = key;											\
			__ac_set_isboth_false(h->flags, x);							\
			++h->size; ++h->n_occupied;									\
			*ret = 1;													\
		} else if (__ac_isdel(h->flags, x)) {							\
			h->keys[x] = key;											\
			__ac_set_isboth_false(h->flags, x);							\
			++h->size;													\
			*ret = 2;													\
		} else *ret = 0;												\
		return x;														\
	}																	\
	SCOPE void kh_del_##name(kh_##name##_t *h, khint_t x)				\
	{																	\
		if (x != h->n_buckets && !__ac_iseither(h->flags, x)) {			\
			__ac_set_isdel_true(h->flags, x);							\
			--h->size;													\
		}																\
	}

#define KHASH_INIT(name, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal) \
	KHASH_INIT2(name, static inline, khkey_t, khval_t, kh_is_map, __hash_func, __hash_equal)

/* --- BEGIN OF HASH FUNCTIONS --- */

/*! @function
  @abstract     Integer hash function
  @param  key   The integer [khint32_t]
  @return       The hash value [khint_t]
 */
#define kh_int_hash_func(key) (khint32_t)(key)
/*! @function
  @abstract     Integer comparison function
 */
#define kh_int_hash_equal(a, b) ((a) == (b))
/*! @function
  @abstract     64-bit integer hash function
  @param  key   The integer [khint64_t]
  @return       The hash value [khint_t]
 */
#define kh_int64_hash_func(key) (khint32_t)((key)>>33^(key)^(key)<<11)
/*! @function
  @abstract     64-bit integer comparison function
 */
#define kh_int64_hash_equal(a, b) ((a) == (b))
/*! @function
  @abstract     const char* hash function
  @param  s     Pointer to a null terminated string
  @return       The hash value
 */
static inline khint_t __ac_X31_hash_string(const char *s)
{
	khint_t h = *s;
	if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
	return h;
}
/*! @function
  @abstract     Another interface to const char* hash function
  @param  key   Pointer to a null terminated string [const char*]
  @return       The hash value [khint_t]
 */
#define kh_str_hash_func(key) __ac_X31_hash_string(key)
/*! @function
  @abstract     Const char* comparison function
 */
#define kh_str_hash_equal(a, b) (strcmp(a, b) == 0)

/* --- END OF HASH FUNCTIONS --- */

/* Other necessary macros... */

/*!
  @abstract Type of the hash table.
  @param  name  Name of the hash table [symbol]
 */
#define khash_t(name) kh_##name##_t

/*! @function
  @abstract     Initiate a hash table.
  @param  name  Name of the hash table [symbol]
  @return       Pointer to the hash table [khash_t(name)*]
 */
#define kh_init(name) kh_init_##name()

/*! @function
  @abstract     Destroy a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
 */
#define kh_destroy(name, h) kh_destroy_##name(h)

/*! @function
  @abstract     Reset a hash table without deallocating memory.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
 */
#define kh_clear(name, h) kh_clear_##name(h)

/*! @function
  @abstract     Resize a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  s     New size [khint_t]
 */
#define kh_resize(name, h, s) kh_resize_##name(h, s)

/*! @function
  @abstract     Insert a key to the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Key [type of keys]
  @param  r     Extra return code: 0 if the key is present in the hash table;
                1 if the bucket is empty (never used); 2 if the element in
				the bucket has been deleted [int*]
  @return       Iterator to the inserted element [khint_t]
 */
#define kh_put(name, h, k, r) kh_put_##name(h, k, r)

/*! @function
  @abstract     Retrieve a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Key [type of keys]
  @return       Iterator to the found element, or kh_end(h) is the element is absent [khint_t]
 */
#define kh_get(name, h, k) kh_get_##name(h, k)

/*! @function
  @abstract     Remove a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  k     Iterator to the element to be deleted [khint_t]
 */
#define kh_del(name, h, k) kh_del_##name(h, k)


/*! @function
  @abstract     Test whether a bucket contains data.
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       1 if containing data; 0 otherwise [int]
 */
#define kh_exist(h, x) (!__ac_iseither((h)->flags, (x)))

/*! @function
  @abstract     Get key given an iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       Key [type of keys]
 */
#define kh_key(h, x) ((h)->keys[x])

/*! @function
  @abstract     Get value given an iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @param  x     Iterator to the bucket [khint_t]
  @return       Value [type of values]
  @discussion   For hash sets, calling this results in segfault.
 */
#define kh_val(h, x) ((h)->vals[x])

/*! @function
  @abstract     Alias of kh_val()
 */
#define kh_value(h, x) ((h)->vals[x])

/*! @function
  @abstract     Get the start iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       The start iterator [khint_t]
 */
#define kh_begin(h) (khint_t)(0)

/*! @function
  @abstract     Get the end iterator
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       The end iterator [khint_t]
 */
#define kh_end(h) ((h)->n_buckets)

/*! @function
  @abstract     Get the number of elements in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       Number of elements in the hash table [khint_t]
 */
#define kh_size(h) ((h)->size)

/*! @function
  @abstract     Get the number of buckets in the hash table
  @param  h     Pointer to the hash table [khash_t(name)*]
  @return       Number of buckets in the hash table [khint_t]
 */
#define kh_n_buckets(h) ((h)->n_buckets)

/* More conenient interfaces */

/*! @function
  @abstract     Instantiate a hash set containing integer keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_INT(name)										\
	KHASH_INIT(name, khint32_t, char, 0, kh_int_hash_func, kh_int_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing integer keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_INT(name, khval_t)								\
	KHASH_INIT(name, khint32_t, khval_t, 1, kh_int_hash_func, kh_int_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_INT64(name)										\
	KHASH_INIT(name, khint64_t, char, 0, kh_int64_hash_func, kh_int64_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_INT64(name, khval_t)								\
	KHASH_INIT(name, khint64_t, khval_t, 1, kh_int64_hash_func, kh_int64_hash_equal)

typedef const char *kh_cstr_t;
/*! @function
  @abstract     Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
 */
#define KHASH_SET_INIT_STR(name)										\
	KHASH_INIT(name, kh_cstr_t, char, 0, kh_str_hash_func, kh_str_hash_equal)

/*! @function
  @abstract     Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
  @param  khval_t  Type of values [type]
 */
#define KHASH_MAP_INIT_STR(name, khval_t)								\
	KHASH_INIT(name, kh_cstr_t, khval_t, 1, kh_str_hash_func, kh_str_hash_equal)

#endif /* __AC_KHASH_H */
/* The MIT License

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Last Modified: 05MAR2012 */

#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

#define __KS_TYPE(type_t)						\
	typedef struct __kstream_t {				\
		unsigned char *buf;						\
		int begin, end, is_eof;					\
		type_t f;								\
	} kstream_t;

#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

#define __KS_BASIC(type_t, __bufsize)								\
	static inline kstream_t *ks_init(type_t f)						\
	{																\
		kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));	\
		ks->f = f;													\
		ks->buf = (unsigned char*)malloc(__bufsize);				\
		return ks;													\
	}																\
	static inline void ks_destroy(kstream_t *ks)					\
	{																\
		if (ks) {													\
			free(ks->buf);											\
			free(ks);												\
		}															\
	}

#define __KS_GETC(__read, __bufsize)						\
	static inline int ks_getc(kstream_t *ks)				\
	{														\
		if (ks->is_eof && ks->begin >= ks->end) return -1;	\
		if (ks->begin >= ks->end) {							\
			ks->begin = 0;									\
			ks->end = __read(ks->f, ks->buf, __bufsize);	\
			if (ks->end < __bufsize) ks->is_eof = 1;		\
			if (ks->end == 0) return -1;					\
		}													\
		return (int)ks->buf[ks->begin++];					\
	}

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define __KS_GETUNTIL(__read, __bufsize)								\
	static int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append) \
	{																	\
		if (dret) *dret = 0;											\
		str->l = append? str->l : 0;									\
		if (ks->begin >= ks->end && ks->is_eof) return -1;				\
		for (;;) {														\
			int i;														\
			if (ks->begin >= ks->end) {									\
				if (!ks->is_eof) {										\
					ks->begin = 0;										\
					ks->end = __read(ks->f, ks->buf, __bufsize);		\
					if (ks->end < __bufsize) ks->is_eof = 1;			\
					if (ks->end == 0) break;							\
				} else break;											\
			}															\
			if (delimiter == KS_SEP_LINE) { \
				for (i = ks->begin; i < ks->end; ++i) \
					if (ks->buf[i] == '\n') break; \
			} else if (delimiter > KS_SEP_MAX) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (ks->buf[i] == delimiter) break;					\
			} else if (delimiter == KS_SEP_SPACE) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i])) break;						\
			} else if (delimiter == KS_SEP_TAB) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; \
			} else i = 0; /* never come to here! */						\
			if (str->m - str->l < (size_t)(i - ks->begin + 1)) {		\
				str->m = str->l + (i - ks->begin) + 1;					\
				kroundup32(str->m);										\
				str->s = (char*)realloc(str->s, str->m);				\
			}															\
			memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); \
			str->l = str->l + (i - ks->begin);							\
			ks->begin = i + 1;											\
			if (i < ks->end) {											\
				if (dret) *dret = ks->buf[i];							\
				break;													\
			}															\
		}																\
		if (str->s == 0) {												\
			str->m = 1;													\
			str->s = (char*)calloc(1, 1);								\
		} else if (delimiter == KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; \
		str->s[str->l] = '\0';											\
		return str->l;													\
	} \
	static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) \
	{ return ks_getuntil2(ks, delimiter, str, dret, 0); }

#define KSTREAM_INIT(type_t, __read, __bufsize) \
	__KS_TYPE(type_t)							\
	__KS_BASIC(type_t, __bufsize)				\
	__KS_GETC(__read, __bufsize)				\
	__KS_GETUNTIL(__read, __bufsize)

#define kseq_rewind(ks) ((ks)->last_char = (ks)->f->is_eof = (ks)->f->begin = (ks)->f->end = 0)

#define __KSEQ_BASIC(SCOPE, type_t)										\
	SCOPE kseq_t *kseq_init(type_t fd)									\
	{																	\
		kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));					\
		s->f = ks_init(fd);												\
		return s;														\
	}																	\
	SCOPE void kseq_destroy(kseq_t *ks)									\
	{																	\
		if (!ks) return;												\
		free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s); \
		ks_destroy(ks->f);												\
		free(ks);														\
	}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
#define __KSEQ_READ(SCOPE) \
	SCOPE int kseq_read(kseq_t *seq) \
	{ \
		int c; \
		kstream_t *ks = seq->f; \
		if (seq->last_char == 0) { /* then jump to the next header line */ \
			while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@'); \
			if (c == -1) return -1; /* end of file */ \
			seq->last_char = c; \
		} /* else: the first header char has been read in the previous call */ \
		seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ \
		if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1; /* normal exit: EOF */ \
		if (c != '\n') ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */ \
		if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */ \
			seq->seq.m = 256; \
			seq->seq.s = (char*)malloc(seq->seq.m); \
		} \
		while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') { \
			if (c == '\n') continue; /* skip empty lines */ \
			seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ \
			ks_getuntil2(ks, KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */ \
		} \
		if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	\
		if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ \
			seq->seq.m = seq->seq.l + 2; \
			kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ \
			seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
		} \
		seq->seq.s[seq->seq.l] = 0;	/* null terminated string */ \
		if (c != '+') return seq->seq.l; /* FASTA */ \
		if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ \
			seq->qual.m = seq->seq.m; \
			seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
		} \
		while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */ \
		if (c == -1) return -2; /* error: no quality string */ \
		while (ks_getuntil2(ks, KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l); \
		seq->last_char = 0;	/* we have not come to the next header line */ \
		if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ \
		return seq->seq.l; \
	}

#define __KSEQ_TYPE(type_t)						\
	typedef struct {							\
		kstring_t name, comment, seq, qual;		\
		int last_char;							\
		kstream_t *f;							\
	} kseq_t;

#define KSEQ_INIT2(SCOPE, type_t, __read)		\
	KSTREAM_INIT(type_t, __read, 16384)			\
	__KSEQ_TYPE(type_t)							\
	__KSEQ_BASIC(SCOPE, type_t)					\
	__KSEQ_READ(SCOPE)

#define KSEQ_INIT(type_t, __read) KSEQ_INIT2(static, type_t, __read)

#define KSEQ_DECLARE(type_t) \
	__KS_TYPE(type_t) \
	__KSEQ_TYPE(type_t) \
	extern kseq_t *kseq_init(type_t fd); \
	void kseq_destroy(kseq_t *ks); \
	int kseq_read(kseq_t *seq);

#endif
