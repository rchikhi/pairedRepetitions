/* C version of pairs-v7.py
 *
 * computes paired repetitions within a genome sequence 
 *
 * for a genome length of n,
 * space complexity: C*n*delta where C is a multiple of 8 (since sizeof(int)=64 bits)
 * time complexity: O(n*delta), observed a 17x speedup over python
 *
 * Note: we discard undetermined nucleotides ('N'). 
 *
 * Further optimization: replace v by not_v and use a Judy1 array
 *
 * */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

#include "my_sarray/sarray.h"
#include "33bits.h"
#include "short2judy.h"
#include "p-judy.h"
#include "argtable2-11/src/argtable2.h"

/* Default parameters */
int print_length=16;
int span=300;
int delta=10;
int sigma_delta_unicity=1;
int hd=0;
int verbose=0;
int mkesa=0;
int nopairs=0;

//char *filename="test_paired_seq_s50_2,4.fna";
char *filename="NC_001416.fna";
//char *filename="../bigdata/drerio_ref_chr25.fa";
//char *filename="NC_000913.fna";
//char *filename="../bigdata/oan_ref_chr10.fa";
//char *filename="randomseq.fna";

unsigned char *seq; //is defined in p-judy.h
my_int length; //same
unsigned int no_r_below; // same
// my_int is defined in sarray.h
my_int *sa;
my_lc *lc2;
my_sj *lc;
#define LC(idx) (get_sj(lc,(idx)))
my_int len_seq;
my_int len_sa;
my_33bits **h;
my_int paired_dupes=0;
my_int *len_h;
my_int **p;
my_r *r;
my_33bits *invsa;
#define INVSA(idx) (get_33bits(invsa,(idx)))
my_int paired_dupes;
my_int left;
my_int right;
my_int maxl;
bool *d;
bool *v;
#define vINVSA(idx) (hd?v[(idx)]:v[INVSA((idx))])
my_int dupes;
#define MAXL_BOUND 165535

void error(char * s)
{
   printf("%s \n",s);
   exit(1);
}

void verif_alloc(void *ptr)
{
   if (ptr==NULL)
	error("Failed to allocate memory");
}


void logging(char * s)
{
      if (!verbose) return;
      time_t now;
      time(&now);
      struct tm *tim;
      tim=localtime(&now);
      printf("%2d:%02d - %s..\n",tim->tm_hour,tim->tm_min,s);
      fflush(stdout);
}

my_int open_file()
{

   // opens the file
   FILE *f=fopen(filename,"r");
   if (f==NULL) return 0;
   fseek(f, 0, SEEK_END);
   my_int size = ftell(f); 
   fseek(f, 0, SEEK_SET);

   seq=(char *)malloc((2*size+1)*sizeof(unsigned char));
   char *block=malloc(size);
   if (seq==NULL || block==NULL)
      error("Failed to allocate memory (already?)");

   // read the file (in one block, is that good?)
   len_seq=0;
   int block_read;
   block_read=fread(block,1,size,f);
   int j=0;
   unsigned char current;
  if (block_read==0) return 0;
   while (j<block_read)
   {
      while (block[j]=='>')
      {
	 while (block[j]!='\n') j++;
	 j++;
      }
      current=block[j];
      if (current=='A' || current=='C' || current=='G' || current=='T') //|| current=='N')
	 seq[len_seq++]=current;
      if (current=='a')
	 seq[len_seq++]='A';
      else if (current=='c')
	 seq[len_seq++]='C';
      else if (current=='g')
	 seq[len_seq++]='G';
      else if (current=='t')
	 seq[len_seq++]='T';
      j++;
   }
   free(block);
   seq[len_seq]=0;

   // compute the reverse complement
   strcat(seq,"N");

   my_int i;
   for (i=0;i<len_seq;i++)
   {
      current=seq[len_seq-1-i];
      if (current=='A' || current=='a') seq[len_seq+1+i]='T';
      else if (current=='C' || current=='c') seq[len_seq+1+i]='G';
      else if (current=='T' || current=='t') seq[len_seq+1+i]='A';
      else if (current=='G' || current=='g') seq[len_seq+1+i]='C';
//      else seq[len_seq+1+i]='N';
   }
   seq[2*len_seq+1]=0;

   fclose(f);

   if (mkesa)
   {
     char filename2[100];
     sprintf(filename2,"mkesafiles/%s.both",filename);
     f=fopen(filename2,"r");
     if (f==NULL)
     {
       logging("building sequence with both strands for mkesa");
       f=fopen(filename2,"w");
       fputs("> ",f);fputs(filename,f);fputs("\n",f);
       for (i=0;i<2*len_seq+1;i+=80)
        {
	 fwrite(seq+i,1,(i+80>=2*len_seq+1)?(2*len_seq+1-i):80,f);
         fputs("\n",f);
        }
       fclose(f);
     }
   }

   return (2*len_seq+1);
}

void print_time()
{
   printf("\n TIME: %d\n",clock());fflush(stdout);
}


void build_h()
{
   my_int *count;
   count=malloc(MAXL_BOUND*sizeof(my_int));
   verif_alloc(count);
   memset(count,0,MAXL_BOUND*sizeof(my_int));
   my_int vmax=0;
   my_int i;
   for (i=0;i<len_sa;i++)
   {
      if (LC(i)>vmax)
      {
	 vmax=LC(i);
	 if (vmax>=MAXL_BOUND) 
	 {
	    printf("try a larger MAXL_BOUND (%d)",vmax);
	    exit(1);
	 }
      }
      count[LC(i)]++;
   }
   maxl=vmax;
   h=(my_33bits**)malloc((maxl+1)*sizeof(my_33bits*));
   len_h=(my_int *)malloc((maxl+1)*sizeof(my_int));
   verif_alloc(h);
   verif_alloc(len_h);

   for (i=0;i<=maxl;i++)
   {
      h[i]=(my_33bits*)alloc_33bits(count[i]);
      len_h[i]=count[i];
      count[i]=0;
   }
   for (i=0;i<len_sa;i++)
      set_33bits(h[LC(i)],count[LC(i)]++,i);
   free(count);
}

void free_h()
{
   free(len_h);
   my_int i;
   free(h[0]);
   free(h);
}


void update_v(int length)
{

   my_int j;
   my_int i;
   for (i=0;i<len_h[length];i++)
   {
      j=get_33bits(h[length],i);
      if (!v[j])
      {
	 dupes+=1;
	 // major change from .py version: the d array is never used here
	 // d[sa[j]]=1;
	 v[j]=1;
      }
      if (!v[j-1])
      {
	 dupes+=1;
	 // d[sa[j-1]]=1;
	 v[j-1]=1;
      }
   }
   free(h[length]); // h[length] wont actually be used anymore
}


// r is an array that links every identical subword of length=length to a single one
// trick: for large subwords, it stores a suffix array offset 
// for small subwords, a counter
#define SWITCH_R 8 // 9 may not work for 16bits, monitor maxr to see that
void build_r(int length)
{
   my_int lastr=0;
   my_int i;
   my_int r_value; 
   my_int maxr=0;

   for (i=0;i<len_sa;i++)
   {
	   if (length>=SWITCH_R)
	   {
	   	if (LC(i)<length) 
		   lastr=i;
	        r_value=i-lastr;
	   } 
	   else
	   {
	   	if (LC(i)<length) 
		   lastr++;
	        r_value=lastr;
	   }
	   r[i]=(my_r)r_value;
	   if (r_value>maxr) maxr=r_value;
	   if (r_value>=UINT32_MAX)
	   {
		   printf("try a larger my_r type for the r array\n");
		   printf("length: %d\tmaxr: %d\t SWITCH_R: %d",length,maxr,SWITCH_R);
		   exit(1);
	   }
   }
//	printf("maxr %d\n",maxr);
}
#define R(idx) length>no_r_below?(\
	       (length<SWITCH_R?r[INVSA((idx))]:(my_int)INVSA((idx))-(my_int)r[INVSA((idx))]))\
		:(idx)


void get_range(my_int *left,my_int *right,my_int i,int length)  
{
   if (!delta)
   {
      *left=i+length+span;
      *right=i+length+span;
      return;
   }
   *left=i+length+span-delta;
   if ((i+length+span)<delta)
      *left=length+span-delta;
   else if (*left<(len_seq+1+length+span) && *left>len_seq)
      *left=len_seq+1+length+span-delta;
   *right=i+length+span+delta;
   if (*right>(len_sa-length))
      *right=len_sa-length;
   else if (*right>(len_seq-length) && *right<(len_seq+length+span))
      *right=len_seq-length;
}


void init_p(my_int len_sa);
void init_p_entries(my_int len_sa);
my_int p_size(my_int len_sa);
void free_p();
void free_p_entries(my_int len_sa);
void set_p(my_int j,my_int k,int v);
int return_p(my_int j,my_int k);

#define ALLOC_SA ((len_sa+1)*sizeof(my_int))
#define ALLOC_LC ((len_sa+1)*sizeof(my_lc))
#define ALLOC_R ((len_sa+1)*sizeof(my_r))



// vfiles: use only v[], which is stored on disk
// highly memory-effective
// but works only when we use p-judy-seq.c and length<=16

bool check_vfiles()
{
  char vfilename[100];
  bool vfiles=0;
  sprintf(vfilename,"vfiles/%s.1.v",filename);
  FILE *f=fopen(vfilename,"r");
  if (f!=NULL)
  { 
    vfiles=1;
    fclose(f);
  }
  return vfiles;
}

void create_vfiles()
{
  char vfilename[100];
  my_int i;
  unsigned char *boolv;
  boolv=(unsigned char *)malloc((len_sa+1)/8);
  verif_alloc(boolv);
    for (length=maxl;length>0;length--)
    {
      update_v(length);

      if (length>0 && length<=print_length)
      {
	sprintf(vfilename,"vfiles/%s.%d.v",filename,length);
	FILE *f=fopen(vfilename,"w");
	if (f==NULL) error("Cannot open v-file");
	memset(boolv,0,(len_sa+1)/8);
	for (i=0;i<len_sa;i++)
	  boolv[i/8]+=(v[INVSA(i)]&1)<<(i%8);
	fwrite(boolv,1,(len_sa+1)/8,f);
	fclose(f);
      }
    }
  free(boolv);
}

void load_vfiles(length)
{
  char vfilename[100];
   unsigned char *boolv;
  boolv=(unsigned char *)malloc((len_sa+1)/8);
  verif_alloc(boolv);
  my_int i;
      sprintf(vfilename,"vfiles/%s.%d.v",filename,length);
      FILE *f=fopen(vfilename,"r");
      if (f==NULL) error("Cannot open v-file");
      int block_read;
      block_read=fread(boolv,1,(len_sa+1)/8,f);
      if (block_read!=(len_sa+1)/8) error("Corrupt v-file!");
      for (i=0;i<len_sa;i++)
	v[i]=((boolv[i/8]>>(i%8))&1);
      fclose(f);
   free(boolv);
}

void load_mkesa_lcp()
{
  my_int i;
  char filename2[100];
  sprintf(filename2,"mkesafiles/%s.lcp",filename);
  FILE *f=fopen(filename2,"r");
  if (f==NULL)
  {  error("Cannot open lcp file from mkesa");}
  else
  {
    lc = (my_sj*)alloc_sj(len_sa);
    unsigned char currentlcp[1];
    for (i=0;i<len_sa;i++)
    {
      fread(currentlcp,1,1,f);
      set_sj(lc,i,currentlcp[0]);
    }
    fclose(f);
    sprintf(filename2,"mkesafiles/%s.llv",filename);
    FILE *llv=fopen(filename2,"r");
    if (llv==NULL)
    { error("Cannot open large lcp file from mkesa");}
    else
      while (!feof(llv))
    {
      unsigned char elt[8];
      fread(elt,1,8,llv);
      my_int pos=(uint64_t)elt[0]+((uint64_t)elt[1]<<8)+((uint64_t)elt[2]<<16)+\
		 ((uint64_t)elt[3]<<24)+((uint64_t)elt[4]<<32)+((uint64_t)elt[5]<<40);
      fread(elt,1,8,llv);
      my_int val=(uint64_t)elt[0]+((uint64_t)elt[1]<<8)+((uint64_t)elt[2]<<16)+\
		 ((uint64_t)elt[3]<<24)+((uint64_t)elt[4]<<32)+((uint64_t)elt[5]<<40);
      set_sj(lc,pos,val);
    }
    fclose(llv);
  }
}

void load_mkesa_invsa()
{
  my_int i=0;
  char filename2[100];
  sprintf(filename2,"mkesafiles/%s.suf",filename);
  FILE *f=fopen(filename2,"r");
  if (f==NULL)
  {  error("Cannot open suffix array (.suf) file from mkesa");}
  else
  {
    while (!feof(f))
    {
      unsigned char elt[8];
      fread(elt,1,8,f);
      my_int val=(uint64_t)elt[0]+((uint64_t)elt[1]<<8)+((uint64_t)elt[2]<<16)+\
		 ((uint64_t)elt[3]<<24)+((uint64_t)elt[4]<<32)+((uint64_t)elt[5]<<40);
      set_33bits(invsa,val,i++);
    } 
  }
  fclose(f);
}


void print_mem(int length)
{
       if (!verbose) return;
	printf("---------- memory usage:   ");
	my_int total=0;

        if (!hd)
	{ 
	  printf("lc: %d Mb\t", size_sj(lc,len_sa) / 1000000);
	  total+=size_sj(lc,len_sa);
	  printf("invsa: %d Mb\t",size_33bits(invsa,len_sa) / 1000000);
	  total+=size_33bits(invsa,len_sa); 

	  my_int h_size=0;
	  h_size+=(maxl+1)*sizeof(my_int); //len_h
	  my_int i;
	  for (i=0;i<length;i++) // not to maxl, because we free h[i] progressively
	    h_size+=size_33bits(h[i],len_h[i]); // size of h[i]
	  printf("h: %d Mb\t",h_size / 1000000);
	  total+=h_size;

	}

        my_int psize=p_size(len_sa);
	printf("p: %d Mb\t",psize / 1000000);
	total+=psize;

	if (!hd && length>no_r_below)
	{	
		printf("r: %d Mb\t",(ALLOC_R) / 1000000);
		total+=ALLOC_R;
	}
	if (no_r_below>0)
	{
		printf("seq: %d Mb\t",strlen(seq) / 1000000);
		total+=strlen(seq);
	}

	printf("v: %d Mb\t",(len_sa)*sizeof(bool) / 1000000);
	total+=((len_sa)*sizeof(bool));

	printf("\n-----------Total: %d Mb (sizeof(my_int)=%d)\n",total/1000000,sizeof(my_int));
}

/* 
 *
 *           _                     ______ 
            (_)                   |____  |
 _ __   __ _ _ _ __ ___ ________   __ / / 
| '_ \ / _` | | '__/ __|______\ \ / // /  
| |_) | (_| | | |  \__ \       \ V // /   
| .__/ \__,_|_|_|  |___/        \_//_/    
| |                                       
|_|          
 *
 *
 */

int main(int argc, char** argv)
{

   // cmdline processing
    struct arg_int *cmd_sigma  = arg_int0("s","sigma",NULL,"define sigma value (default is 300)");
    struct arg_int *cmd_delta  = arg_int0("d","delta",NULL,"define delta value (default is 0)");
    struct arg_file *cmd_filename  = arg_file0("f",NULL,"<filename>","sequence to analyze");
    struct arg_lit  *cmd_hd    = arg_lit0(NULL,"hd","stores v[] on disk (in ./vfiles/)");
    struct arg_lit  *cmd_mkesa    = arg_lit0(NULL,"mkesa","use mkesa to build the suffix array and lcp (in ./mkesafiles/)");
    struct arg_int *cmd_printlength  = arg_int0("l","length",NULL,"set maximum read length (default is 16)");
    struct arg_lit  *cmd_nopairs    = arg_lit0(NULL,"nopairs","do not compute paired uniqueness, only single uniqueness");
    struct arg_lit  *cmd_verb    = arg_lit0("v","verbose","verbose");
    struct arg_end  *cmd_end     = arg_end(20);
    void *argtable[]={cmd_sigma,cmd_delta,cmd_filename,cmd_hd,cmd_mkesa,cmd_printlength,cmd_nopairs,cmd_verb,cmd_end};
    int nerrors;
    cmd_sigma->ival[0]=span;
    cmd_delta->ival[0]=delta;
    cmd_filename->filename[0]=filename;
    cmd_printlength->ival[0]=print_length;
    nerrors = arg_parse(argc,argv,argtable);
    if (nerrors>0)
    { 
     printf("Usage: pairs-v7");
     arg_print_syntax(stdout,argtable,"\n");
     arg_print_glossary(stdout,argtable,"  %-25s %s\n");
     exit(1);
    }
    span=cmd_sigma->ival[0];
    delta=cmd_delta->ival[0];
    filename=(char *)cmd_filename->filename[0];
    hd|=cmd_hd->count;
    verbose|=cmd_verb->count;
    mkesa|=cmd_mkesa->count;
    nopairs|=cmd_nopairs->count;
    print_length=cmd_printlength->ival[0];
    printf("File=%s (sigma,delta)=(%d,%d) maxlength=%d %s%s\n",filename,span,delta,print_length,\
            hd?"using vfiles":"not using vfiles",\
            mkesa?", using mkesa":", not using mkesa");

   // opens sequence
   if (!(len_sa=open_file())) 
      error("Can't open file");
   logging("file opened");
   printf("Sequence length (both strands): %d\n",len_sa);

   my_int i;
   no_r_below=0;
   v=NULL;

   // if --hd, check if we can use vfiles, if not, create everything
   if (hd==0 || check_vfiles()==0)
   {
     if (hd) logging("creating suffix array, lcp, etc.. for vfiles");

     if (!mkesa)
     {
       // alloc memory for suffix array
       sa = (my_int*)malloc(ALLOC_SA);
       lc2 = (my_lc*)malloc((len_sa+1)*sizeof(my_lc));
       if (sa==NULL || lc2==NULL)
	 error("Failed to allocate memory (already?)");

       // create suffix array and lcp
       bsarray(seq,sa,len_sa);
       logging("suffix array created");
       lcpa(sa,seq,lc2,len_sa+1);
       logging("lcp created");

       // replaces 64-bit lcp by a 32bits one
       lc = (my_sj*)alloc_sj(len_sa);
       for (i=0;i<len_sa;i++)
	 set_sj(lc,i,lc2[i+1]);
       free(lc2);

       // create the inverse of the suffix array
       invsa = (my_33bits*)alloc_33bits(len_sa);
       for (i=0;i<len_sa;i++)
	 set_33bits(invsa,sa[i+1],i);// sa[i+1] because pysarray shifts indices
       free(sa);
     }
     else
     {
       lc = (my_sj*)alloc_sj(len_sa);
       load_mkesa_lcp();
       logging("mkesa lcp loaded");
       invsa = (my_33bits*)alloc_33bits(len_sa);
       load_mkesa_invsa();
       logging("mkesa inverse suffix array loaded");
     }


     // allocate memory for the r array 
     r = (my_r*)malloc((len_sa+1)*sizeof(my_r));
     verif_alloc(r);

     logging("low memory arrays allocated");

     // allocate memory for h and create it
     build_h(); 

     logging("h array built");

     // allocate memory for the v array
     if (hd==1 && check_vfiles()==0)
     {

       v=(bool*)malloc(len_sa*sizeof(bool));
       verif_alloc(v);
       memset(v,0,len_sa);

       logging("creating v-files");
       create_vfiles();
       free_sj(lc);
       free_33bits(invsa);	
       free_h();
 
     }
   }
   
   if (check_vfiles() && hd) logging("using vfiles");
   if (v==NULL) 
   {
     v=(bool*)malloc(len_sa*sizeof(bool));
     verif_alloc(v);
     memset(v,0,len_sa);
   }
   
   // init p (for sieving paired) and rb (for right-mate dupes) array
   init_p(len_sa);
   init_rb();
   dupes=0;
   if (hd==0 && no_r_below==0) free(seq);
   if (hd && check_vfiles()) maxl=16;
   if (hd && no_r_below==0) error("r array needed, cannot use vfiles. Perhaps you should use the 64bit p array.");

   logging("beginning analysis");

   // begins processing
   //
   for (length=maxl;length>0;length--)
   {
      if (!hd && length==no_r_below) free(r);
      if (!hd) update_v(length);

      if (length>0 && length<=print_length)
      {
         
         if (nopairs) { printf("length %d dupes %d (%.2f uniques)\n",length,dupes,(100.0*(len_sa-length+1-dupes)/(len_sa-length+1)));continue;}

         if (hd) load_vfiles(length);

	 init_p_entries(len_sa);
	 paired_dupes=0;
	 if (!hd && length>no_r_below) build_r(length);

	 // pass 1: fill
	 logging("pass 1");
	 for (i=0;i<len_sa;i++)
	 {
	    if (vINVSA(i) && !((i>(len_seq-span+delta-length-length+1) && i<(len_seq+1))
		     || (i>(len_sa-span+delta-length-length))))
	    {
	       get_range(&left,&right,i,length);
	       my_int k;
	       for (k=left;k<=right;k++)
	       {

		  // p is too large when delta!=0..
	          if (delta==0 && !vINVSA(k)) break;
		  unsigned char rp=return_p(R(i),R(k));
		  if (rp!=2)
		     set_p(R(i),R(k),rp+1);
	       }
	    }
	 }
   	 logging("pass 2");
	 // pass 2: count
	 for (i=0;i<len_sa;i++)
	 {
	    if ((i>(len_seq-span+delta-length-length+1) && i<(len_seq+1))
		  || (i>(len_sa-span+delta-length-length)))
	       continue;
	    if (!vINVSA(i))
	    {
	       if (delta==0)
		  continue;
	       get_range(&left,&right,i,length);
	       my_int k;
	       rb_refresh();
	       for (k=left;k<=right;k++)
	       {
		  if (!vINVSA(k)) continue;
		  unsigned char rp=return_rb(R(k));
		  if (rp==0) set_rb(R(k));
		  else { 
		    paired_dupes+=1; 
		    if(sigma_delta_unicity)  break;
		  }
	       }
	      
	       continue;
	    }
	    get_range(&left,&right,i,length);
	    my_int k;
	    for (k=left;k<=right;k++)
	    {
	       if (delta==0 && !vINVSA(k)) break;
	       if (return_p(R(i),R(k))==2)
	       {
		  paired_dupes+=1;
		  if(sigma_delta_unicity)  break;
	       }
            }
	 }
         if (length%2==0) print_mem(length);
	 free_p_entries(len_sa);

        my_int len_range;
	len_range=len_seq-span+delta-length-length+1+(len_sa-span+delta-length-length+1)-(len_seq+1); 
        float uniques=100.0*(len_range*(1+(1-sigma_delta_unicity)*(2*delta))-paired_dupes)/(len_range*(1+(1-sigma_delta_unicity)*(2*delta)));

	printf("length %d\t paired_dupes: %d (%.2f %% uniques)\tdupes: %d (%.2f %% uniques)\n",length,paired_dupes, uniques, dupes,(100.0*(len_sa-length+1-dupes)/(len_sa-length+1)));
      }
   }

   if (no_r_below>0 || hd) free(seq);
   if (!hd)
    {
   free_sj(lc);
   free_33bits(invsa);	
   free_h();
   if (no_r_below==0) free(r);	
   }
   free_p();
   free(v);	
   arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
   return 0;
}
