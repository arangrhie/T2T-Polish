#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>       /* ceil */

#include <zlib.h>
#include "kseq.h"
#include "bamcat.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

/*
  basemap[] works by storing a very small array that maps a base to 
  its complement, by dereferencing the array with the ASCII char's 
  decimal value as the index
  (int) 'A' = 65;
  (int) 'C' = 67;
  (int) 'G' = 71;
  (int) 'T' = 84;
  (int) 'a' = 97;
  (int) 'c' = 99;
  (int) 'g' = 103;
  (int) 't' = 116;
  (int) 'N' = 78;
  (int) 'U' = 85;
  (int) 'u' = 117;
  for example: basemap['A'] => basemap[65] => 'T' etc.
*/

static const char basemap[255] =
    {
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
        '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
        '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
        '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', /*  90 -  99 */
        '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
        '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a', '\0', '\0', /* 110 - 119 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
        '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
    };

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static bool is_mapped(const bam1_core_t *core) {

    if (core->flag & BAM_FUNMAP) {
        return false;
    }

    return true;
}

char *rc(char *input, int len) {
   char *seq = new char[len+1];
   input+=len-1;
   for (int i = 0; i < len; i++) {
      seq[i]=basemap[*input];
      input--;
   }
   return seq;
}

samfile_t * open_alignment_file(std::string path) {
    samfile_t * fp = NULL;
    std::string flag = "r";
    if (path.substr(path.size() - 3).compare("bam") == 0) {
        //BAM file!
        flag += "b";
    }
    if ((fp = samopen(path.c_str(), flag.c_str(), 0)) == 0) {
        fprintf(stderr, "qaCompute: Failed to open file %s\n", path.c_str());
    }
    return fp;
}

/***
 * Usage: samToAlignment in.bam ref.fa
 * Output: 
 */
int
main (int argc, char **argv) {
    map<string, string> reference;
    // loadFasta
    FILE *in = fopen(argv[2], "r");
    gzFile f = gzdopen(fileno(in), "r");
    kseq_t *seq = kseq_init(f);
    int l = 0;
    while ((l=kseq_read(seq)) >= 0) {
       reference[seq->name.s] = seq->seq.s;
       string name = seq->name.s;
       name.append("rc");
       reference[name] = rc(seq->seq.s, strlen(seq->seq.s));
    }
    kseq_destroy(seq);
    gzclose(f);
    fclose(in);

    bam1_t *b = bam_init1();
    samfile_t * fp = open_alignment_file(argv[1]);
    if (fp == NULL) {
        exit(1);
    }
    bam_header_t* head = fp->header; // sam header
    char *qrySeq = NULL;
    char *orig   = NULL;
    string lastid;

    while (samread(fp, b) >= 0) {
        //Get bam core.
        const bam1_core_t *core = &b->core;
        if (core == NULL) {
            printf("Input file is corrupt!");
            exit(1);
        }
        if (!is_mapped(core)) {
            continue;
        }

        string id = bam1_qname(b);
        string ref = head->target_name[core->tid];
        uint32_t* cigar = bam1_cigar(b);
        uint32_t refLen = head->target_len[core->tid];
        uint32_t refLo = core->pos + 1;
        uint32_t refHigh = bam_calend(core, cigar);
        uint32_t seqLow = 0;
        uint32_t seqHi = core->l_qseq; 
        uint32_t seqLen = core->l_qseq;
        bool isFwd = !(core->flag & BAM_FREVERSE);

        // now parse, we need both sequecnes for this
        if (strlen(reference[ref].c_str()) <= 0) { 
           fprintf(stderr, "Error: unknown reference sequence %s\n", ref.c_str());
           exit(1);
        }
        char *refSeq = (char*)reference[ref].c_str()+refLo-1;
        uint32_t qryLen = core->l_qseq;
        if (core->l_qseq == 0)  { // alt alignment, look for it from previous reads
           if (lastid != id) {
              fprintf(stderr, "Error: no sequence defined for %s, it doesn't match previous %s and has no sequence!\n", id.c_str(), lastid.c_str());
              continue;
           }
           if (orig == NULL) {
              fprintf(stderr, "Error: no sequence defined for %s and I haven't seen it before to reuse!\n", id.c_str());
              exit(1);
            }
           qryLen = strlen(orig);
           qrySeq = orig;
        } else {
           if (orig != NULL) { delete[] orig; qrySeq=NULL; }
           qrySeq = new char[core->l_qseq+1];
           orig = qrySeq;
           for (int i = 0; i < core->l_qseq; i++) {
             qrySeq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)];
           }
           qrySeq[core->l_qseq]='\0';
           lastid = id;
        }
        string alnStr;
        uint32_t pos = 0;

        int errors = 0;
        int matches = 0;
        uint32_t len = 0;
        for (int k = 0; k < core->n_cigar; ++k) {
            int cop = cigar[k] & BAM_CIGAR_MASK; // operation
            int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
            switch (cop) {
            // we don't care about clipping, not part of length
            case BAM_CSOFT_CLIP:
                 if (k == 0) seqLow+=cl;
                 else seqHi-=cl;
                 qrySeq+=cl;
                 for (int i = 0; i < cl; i++) {
                    alnStr += "N"; //'-';
                    pos++;
                 }
                 break;
            case BAM_CHARD_CLIP:
                 if (k == 0) seqLow+=cl;
                 else seqHi-=cl;
                 seqLen+=cl;
                 seqHi+=cl;
                 break;
            // we don't care about matches either, not an error
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF: 
              for (int i = 0; i < cl; i++) {
                 if (toupper(*refSeq) == toupper(*qrySeq)) {
                    matches++;
                    alnStr += toupper(*qrySeq); //'.';
                 } else {
                    alnStr += 'N';
                    errors++;
                 }
                 pos++;
                 refSeq++;
                 qrySeq++;
             }
             len+=cl;
             break;
            case BAM_CINS:
               errors+=cl;
               qrySeq+=cl;
               len+=cl;
               break;
            case BAM_CREF_SKIP:
               for (int i = 0; i < cl; i++) {
                  alnStr += 'N'; //'-';
                  pos++;
               }
               refSeq+=cl;
               break;
            case BAM_CDEL:
               for (int i = 0; i < cl; i++) {
                  alnStr += 'N'; //'-';
                  pos++;
               }
               errors+=cl;
               refSeq+=cl;
               len+=cl;
               break;
            default:
                cerr << "Error: unknown base " << cop << " of len " << cl << endl;
                exit(1);
                break;
            }
        }

        double idy = (1-((double) errors / len /*(seqHi-seqLow)*/ /*(refHigh - refLo + 1)*/)) * 100;
        cout << id << "\t" << ref << "\t" << int(ceil(-1*idy/100*len)) << "\t" << idy << "\t0\t" << (isFwd == true ? seqLow : seqLen-seqHi) << "\t" << (isFwd == true ? seqHi : seqLen-seqLow) << "\t" << seqLen << "\t" <<  (isFwd == true ? 0 : 1) << "\t" << (isFwd == true ? refLo : refLen-refHigh) << "\t" << (isFwd == true ? refHigh : refLen-refLo) << "\t" << refLen << "\t" << len << "\t" << (seqHi-seqLow) << "\t" << (refHigh-refLo) << "\t" << matches << "\t" << errors << "\t" << alnStr.length() << "\t" << alnStr << endl;
    }
    samclose(fp);
return 0;
}
