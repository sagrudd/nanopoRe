#include <Rcpp.h>
using namespace Rcpp;
#include<Rinternals.h>
//#include <array>
#include <iostream>
#include <string.h>
#include <string>
using namespace std;
#include <zlib.h>
#include <stdio.h>
#include <time.h>
// #include <regex>


#include <vector>

static std::string fastq_filename;
static int isZipped = 1;
static gzFile gzfile;
static FILE *file;
static time_t start_time;
static int frameiterations = 0;
//static int allowedframeiterations = 25;
static long long int running_bases = 0;
static long int running_fastq = 0;
static time_t reference_time;


struct Fastq_tag {
  char header[1000000];
  char sequence[1000000];
  char delim[1000000];
  char quality[1000000];
};


static Fastq_tag fq;

static int LONGEST_SEQ = 1000000;


inline bool myfile_exists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}


inline int suffix_match(const char *query, const char *suffix) {
  return(strncmp(query + strlen(query) - strlen(suffix), suffix, strlen(suffix)));
}


inline int is_gzipped(std::string query)
{
vector<string> list;
  list.push_back(".gzip");
  list.push_back(".gz");
  for( vector<string>::const_iterator it = list.begin(); it != list.end(); ++it )
  {
    // Rcout << *it << endl;
    std::string suffix = *it;
    if (suffix_match(query.c_str(), suffix.c_str()) == 0)
    {
      return(1);
    }
  }
  return(0);
}



int has_next_fastq()
{
  if (isZipped == 1) {
    printf("hasNext(gzeof=%d)\n", gzeof(gzfile));
    if ((gzfile!=NULL) && (gzeof(gzfile)==0)) {
      return (1);
    }
  } else {
    if ((file!=NULL) && (feof(file)==0)) {
      return (1);
    }
  }
  return (0);
}






int realign_fastq() {
  return(0);
}


int validate_fastq()
{
  int validated_fastq = 1;  // innocent until found otherwise
  if (has_next_fastq()==0) return 0; // run over the end of the file ...


  // does header start with @ - if not there may be misalignment?
  char headdelim = fq.header[0];
  //printf("headdelim(c=%c)\n", headdelim);
  if (headdelim!='@')
  {
    return(realign_fastq());
  }

  // is delim == "+"? - * is used to indicate runover end of file ...
  // printf("fq.delim(%s)->%d\n", fq.delim, strcmp(fq.delim, "+"));
  if (strcmp(fq.delim, "+")!=0)
  {
    return(0); // this is not a valid fastq
  }
  // is sequence length > 0
  // does sequence length == quality length

  if (frameiterations>0) frameiterations=0;

  running_bases += strlen(fq.sequence);
  running_fastq += 1;

  time_t now;
  time(&now);

  //printf ( "%d - %d\n", now, reference_time);

  if ((now - reference_time) >= 1)
  {
    //
    //show_fastq_stats();
    //
    time ( &reference_time);
  }

  /**
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime);
  timeinfo = localtime(&rawtime);
  printf ( "%d - Current local time and date: %s", rawtime, asctime (timeinfo) );
  **/
  return (validated_fastq);
}







int get_next_fastq()
{
  int r1, r2, r3, r4;
  if (isZipped == 1)
  {
    r1 = gzgets(gzfile, fq.header, LONGEST_SEQ)==NULL;
    r2 = gzgets(gzfile, fq.sequence, LONGEST_SEQ)==NULL;
    r3 = gzgets(gzfile, fq.delim, LONGEST_SEQ)==NULL;
    r4 = gzgets(gzfile, fq.quality, LONGEST_SEQ)==NULL;
  }
  else
  {
    r1 = fgets(fq.header, LONGEST_SEQ, file)==NULL;
    r2 = fgets(fq.sequence, LONGEST_SEQ, file)==NULL;
    r3 = fgets(fq.delim, LONGEST_SEQ, file)==NULL;
    r4 = fgets(fq.quality, LONGEST_SEQ, file)==NULL;
  }
  // clip any trailing newlines
  fq.header[strcspn(fq.header, "\r\n")] = 0;
  fq.sequence[strcspn(fq.sequence, "\r\n")] = 0;
  fq.delim[strcspn(fq.delim, "\r\n")] = 0;
  fq.quality[strcspn(fq.quality, "\r\n")] = 0;
  return (validate_fastq());
}




//' parse a fastq file aiming to validate sequences
//'
//' @param x A fastq format DNA/RNA sequence file
//' @export
// [[Rcpp::export]]
std::string fastqValidator(std::string fastq) {

  Rcout << "set_fastq_file==" << fastq << std::endl;
  fastq_filename = fastq;

  // TEST (1) - DOES THE SPECIFIED FILE EXIST
  bool exists = myfile_exists(fastq_filename);
  if (!exists) {
    return("FastqFileNotFound");
  }

  // HOUSE-KEEPING - IS THE FILE COMPRESSED?
  if (is_gzipped(fastq_filename)==1)
  {
    Rcout << "provided with a gzip file - flying decompress initiated" << std::endl;
    isZipped = 1;
  }
  /**   else
  {
    isZipped = 0;
  }
**/


  if (isZipped == 1) {
    gzfile = gzopen(fastq_filename.c_str(), "r");
    Rcout << "gzopen" << std::endl;
  } else {
    file = fopen(fastq_filename.c_str(), "r");
  }
  time(&start_time);

  while (has_next_fastq() == 1)
  {
    if (get_next_fastq()==1)
    {

    }
  }



  return "Hello World";
}



