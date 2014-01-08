/*
Motivation:
	test wether one or more BAM file is sorted.
	returns 0 on success
Author:
	Pierre Lindenbaum PhD
WWW:
	http://plindenbaum.blogspot.com
Contact:
	plindenbaum@yahoo.fr
Reference:
	http://plindenbaum.blogspot.com/2011/02/testing-if-bam-file-is-sorted-using.html
Compilation:
	FLAG64= -m64
	gcc ${FLAG64} -O3 -I ${SAMDIR} -L ${SAMDIR} bamsorted.c -lbam -lz
Reference:
	http://sourceforge.net/mailarchive/message.php?msg_id=26996499
API:
	http://samtools.sourceforge.net/samtools/sam/

*/
#include <stdio.h>
#include "bam.h"
#include "sam.h"

#define NO_REFERENCE -99999

/** returns EXIT_SUCCESS  if the file 'filename' is a sorted BAM file. */
static int test_sort(const char* filename)
 {
 /** returned status */
 int status=EXIT_SUCCESS;
 /** previous genomic position in the reference */
 int32_t prev_pos=-1;
 /** previous reference index */
 int32_t prev_reference=NO_REFERENCE;
 /** did we find an unmapped read ? */
 int unmapped_flag=0;
 samfile_t *fp_in = NULL;
 bam1_t *b=NULL;
 fprintf(stdout,"%s\t",filename);
 fflush(stdout);
 
 fp_in = samopen(filename, "rb", 0);
 if(NULL == fp_in)
	  {
	  fprintf(stdout,"Could not open file.\n");
	  return EXIT_FAILURE;
	  }
 b = bam_init1();
 while(samread(fp_in, b) > 0) /* loop over the records */
	{
        /** current reference is not the previous reference */
        if(b->core.tid != prev_reference)
		{
		/* it is NOT the very first read */
		if(prev_reference!=NO_REFERENCE)
			{
			/* it's a mapped read and we previously found an unmapped read */
			if(b->core.tid!=-1 && unmapped_flag==1)
				{
				fprintf(stdout,"Found Reference[%d]%s after unmapped reads.",
					b->core.tid,
					fp_in->header->target_name[b->core.tid]
					);
				status=EXIT_FAILURE;
			    	break;
				}
			/* it's an unmapped read */
			else if(b->core.tid==-1)
				{
				unmapped_flag=1;
				}
			/* current reference index is lower than the previous reference index */
			else if(b->core.tid < prev_reference )
			    {
			    fprintf(stdout,"Unsorted: ");
			    fprintf(stdout,"Reference[%d]",prev_reference);
		            if(prev_reference>=0)  fprintf(stdout,"=\"%s\"",fp_in->header->target_name[prev_reference]);
			    fprintf(stdout," followed by ");
	 		    fprintf(stdout,"Reference[%d]",b->core.tid);
		            if(b->core.tid>=0)  fprintf(stdout,"=\"%s\"",fp_in->header->target_name[b->core.tid]);
			    status=EXIT_FAILURE;
			    break;
			    }
			}
		prev_reference= b->core.tid;
		prev_pos=-1;
		}
     if(b->core.tid==-1)
		{
		unmapped_flag=1;
		}
     /** current genomic position is lower than the previous genomic position */
     else if(b->core.pos < prev_pos)
		{
		fprintf(stdout,"Unsorted: On Reference[%d]=\"%s\" position=%d after %d.",
			b->core.tid,
			fp_in->header->target_name[b->core.tid],
			b->core.pos,
			prev_pos
			);
		status=EXIT_FAILURE;
		break;
		}
      prev_pos = b->core.pos;
      bam_destroy1(b);
      b = bam_init1();
      }
 if(status==EXIT_SUCCESS) fprintf(stdout,"OK");
 fprintf(stdout,"\n");
 bam_destroy1(b);
 samclose(fp_in);
 return status;
 }

int main(int argc, char *args[])
  {
  int optind=1;
  int status=EXIT_SUCCESS;
  /* loop over the arguments */
  while(optind<argc)
	    {
	    if(strcmp(args[optind],"-h")==0)
		    {
		    fprintf(stdout, "Author: Pierre Lindenbaum PHD. 2011.\n");
		    fprintf(stdout, "Last compilation:%s %s\n",__DATE__,__TIME__);
		    return EXIT_SUCCESS;
		    }
	    else if(strcmp(args[optind],"--")==0)
		    {
		    optind++;
		    break;
		    }
	    else if(args[optind][0]=='-')
		    {
		    fprintf(stderr,"Unnown option: %s\n",args[optind]);
		    return EXIT_FAILURE;
		    }
	    else
		    {
		    break;
		    }
	    ++optind;
	    }

  if(optind==argc)
      {
      fprintf(stderr,"Illegal number of arguments\n");
      return EXIT_FAILURE;
      }
  /* loop over the files */
  while(optind<argc)
	{
	if(test_sort(args[optind++])!=EXIT_SUCCESS)
		{
		status=EXIT_FAILURE;
		}
	}
  return status;
  }
