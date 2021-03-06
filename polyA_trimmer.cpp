/**
 * Copyright INRIA, CNRS, ENS, main contributor: Gustavo Sacomoto.
 * sacomoto@gmail.com
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * http://www.cecill.info.

 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.

 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 */

#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "kseq.h"
#define MAX 1000
#define IS_FASTQ(seq) (seq->qual.s != NULL)
KSEQ_INIT(gzFile, gzread)

int trimm_end(char *seq, int size, int min_length, int max_error, char base)
{
  int i, error = 0;

  for (i = size-1; i >= 0; i--)
  { 
    error += (seq[i] != base);
    if (error > max_error)
      break;
  }
  if ((size - i - 1) < min_length)
    return size;
  else 
    return i;
}

int trimm_begin(char *seq, int size, int min_length, int max_error, char base)
{
  int i, error = 0;
  
  for (i = 0; i < size; i++)
  { 
    error += (seq[i] != base);
    if (error > max_error)
      break;
  }
  if (i < min_length)
    return 0;
  else 
    return i-1;
}

void print_substr(FILE *stream, char *seq, int begin, int end)
{
  for (int i= begin; i < end; i++)
    fprintf(stream, "%c", seq[i]);
  fprintf(stream, "\n");
}

void print_usage_and_exit(char *name)
{
  fprintf(stderr, "\nusage: %s [-i input file fasta/q] [-o output file] [-l min length (20)] [-e max error (3)] [-m min read length (40)] \n", name);
  exit(0);
}

int main(int argc, char **argv)
{
  char *input_fname, *output_fname;
  int error_code, total_trimmed_bases = 0, min_length = 20, max_error = 3, min_read_length = 40;
 
  gzFile fp = NULL;
  FILE *out = NULL;

  int c;
  while ((c = getopt(argc, argv, "i:o:l:e:")) != -1)
  {
    switch (c)
    {
       case 'i':
	 fp = gzopen(optarg, "r"); 
	 if (fp == NULL)
	 {
	   fprintf(stderr, "Problem opening: %s\n", optarg);
	   exit(0);
	 }  
	 break;

       case 'o':
	 out = fopen(optarg, "w");
	 if (out == NULL)
	 {
	   fprintf(stderr, "Problem opening: %s\n", optarg);
	   exit(0);
	 }  
	 break;	 
	 
      case 'l':
	min_length = atoi(optarg);
	break;
	
      case 'e':
	max_error = atoi(optarg);
	break;
	
      default:
	print_usage_and_exit(argv[0]);
      
    }
  }
  if (fp == NULL || out == NULL)
    print_usage_and_exit(argv[0]);
   
  kseq_t *seq = kseq_init(fp); 
  while ((error_code = kseq_read(seq)) >= 0) 
  { 
    int size = strlen(seq->seq.s), begin, end;

    end = trimm_end(seq->seq.s, size, min_length, max_error, 'A');
    if (end == size)
      end = trimm_end(seq->seq.s, size, min_length, max_error, 'T');
    
    begin = trimm_begin(seq->seq.s, size, min_length, max_error, 'A');
    if (begin == 0)
      begin = trimm_begin(seq->seq.s, size, min_length, max_error, 'T');

    if ((end - begin) >= min_read_length)
    {  
      fprintf(out, "%c%s\n", IS_FASTQ(seq) ? '@' : '>', seq->name.s);
      print_substr(out, seq->seq.s, begin, end);
      if (IS_FASTQ(seq))
      {
	fprintf(out, "+%s\n", seq->name.s);
	print_substr(out, seq->qual.s, begin, end);
      }
      total_trimmed_bases += (size - (end-begin));
    }
    else
      total_trimmed_bases += size;
  }
  if (error_code == -2)
    fprintf(stderr, "Truncated quality string!\n");
  printf("Total trimmed bases: %d\n", total_trimmed_bases);
  
  kseq_destroy(seq); 
  gzclose(fp);
  fclose(out);
  return 0;
}
