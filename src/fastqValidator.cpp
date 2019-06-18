#include <Rcpp.h>
using namespace Rcpp;
#include<Rinternals.h>
#include <iostream>
#include <string.h>
#include <string>
using namespace std;
#include <zlib.h>
#include <stdio.h>
#include <vector>

static std::string fastq_filename;
static int isZipped = 0;
static gzFile gzfile;
static FILE *file;
static time_t start_time;
static bool validFastq = true;
static long long int running_bases = 0;
static long int running_fastq = 0;
static int malformed_fastq_delim = 0;
static int fastq_plus_error = 0;
static int zeroLengthSequence = 0;
static int sequenceQualityLengthMismatch = 0;
static int linesSkipped = 0;

struct Fastq_tag {
  char header[10000000];
  char sequence[10000000];
  char delim[10000000];
  char quality[10000000];
};

static Fastq_tag fq;

static int LONGEST_SEQ = 10000000;

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
    // Rcout << "hasNext(gzeof==" << gzeof(gzfile) << ")" << std::endl;
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
  /**
   * the aim of this is to endeavour to resynchronise a fastq file that may contain an INS or DEL
   * resulting in the missing header and SEQ/QUAL separator
   */

  if (has_next_fastq())
  {

    linesSkipped ++;

    long filepos = 0;
    if (isZipped == 1)
    {
      filepos = gztell(gzfile);
    } else {
      filepos = ftell(file);
    }

  //Rcout << "trying to realign fastq ... " << filepos << std::endl;

  if (isZipped == 1)
  {
    gzgets(gzfile, fq.header, LONGEST_SEQ);
  } else {
    fgets(fq.header ,LONGEST_SEQ , file);
  }

  char headdelim = fq.header[0];
  if (headdelim=='@') {
    // we may be aligned ??
    // but @ is a valid quality score (phred=31) ... it therefore makes some sense just to check ???

    // KEEPING THIS FOR A FUTURE UPDATE ???
    // this does not break functionality, but a block selected on basis of a QUAL @ would sort itself out pretty quickly?

    if (isZipped == 1)
    {
      gzseek(gzfile, filepos, SEEK_SET);
    } else {
      fseek(file, filepos, SEEK_SET);
    }


  } else {
    return(realign_fastq());
  }

  }
  return(0);
}


int validate_fastq()
{
  int validated_fastq = 1;  // innocent until found otherwise
  if (has_next_fastq()==0) return 0; // run over the end of the file ...


  // does header start with @ - if not there may be misalignment?
  char headdelim = fq.header[0];
  if (headdelim!='@')
  {
    Rcout << "Malformed fastq entry delimitter " << headdelim << "!=@" << std::endl;
    malformed_fastq_delim ++;
    validFastq = false;
    return(realign_fastq());
  }

  // is delim == "+"? - * is used to indicate runover end of file ...
  if (strcmp(fq.delim, "+")!=0)
  {
    Rcout << "Malformed fastq entry ~ Line3[+] not present" << std::endl;
    fastq_plus_error ++;
    validFastq = false;
    return(0); // this is not a valid fastq
  }

  // is sequence length > 0
  size_t seqlength = strlen(fq.sequence);
  size_t qlength = strlen(fq.quality);
  if (seqlength == 0) {
    Rcout << "Malformed fastq entry ~ length(sequence)==0" << std::endl;
    zeroLengthSequence ++;
    validFastq = false;
    return(0); // this is not a valid fastq
  }

  // Rcout << "(s/q)Lengths==" << seqlength << "/" << qlength << std::endl;

  // does sequence length == quality length
  if (seqlength != qlength) {
    Rcout << "Malformed fastq entry ~ length(sequence)!=length(quality)" << std::endl;
    sequenceQualityLengthMismatch ++;
    validFastq = false;
    return(0);
  }

  running_bases += seqlength;
  running_fastq += 1;

  return (validated_fastq);
}







int get_next_fastq()
{
  if (isZipped == 1)
  {
    gzgets(gzfile, fq.header, LONGEST_SEQ);
    gzgets(gzfile, fq.sequence, LONGEST_SEQ);
    gzgets(gzfile, fq.delim, LONGEST_SEQ);
    gzgets(gzfile, fq.quality, LONGEST_SEQ);
  }
  else
  {
    fgets(fq.header, LONGEST_SEQ, file);
    fgets(fq.sequence, LONGEST_SEQ, file);
    fgets(fq.delim, LONGEST_SEQ, file);
    fgets(fq.quality, LONGEST_SEQ, file);
  }
  // clip any trailing newlines
  fq.header[strcspn(fq.header, "\r\n")] = 0;
  fq.sequence[strcspn(fq.sequence, "\r\n")] = 0;
  fq.delim[strcspn(fq.delim, "\r\n")] = 0;
  fq.quality[strcspn(fq.quality, "\r\n")] = 0;
  return (validate_fastq());
}





//' return the number of fastq entries previously parsed from provided Fastq file
//'
//' @return long integer of read fastq elements
//' @export
// [[Rcpp::export]]
long int getFastqCount()
{
  return(running_fastq);
}


//' return the number of fastq bases previously parsed from provided Fastq file
//'
//' @return long long integer of read fastq bases
//' @export
// [[Rcpp::export]]
long long int getFastqBases()
{
  return(running_bases);
}

//' count of fastq elements rejected due to malformed fastq header
//'
//' @return integer count of malformed fastq header entries
//' @export
// [[Rcpp::export]]
int getMalformedFastqHeaderCount()
{
  return(malformed_fastq_delim);
}

//' count of fastq elements rejected due to line 3 plus separator
//'
//' @return integer count of malformed fastq '+' separator entries
//' @export
// [[Rcpp::export]]
int getFastqPlusErrorCount()
{
 return(fastq_plus_error);
}

//' count of fastq elements rejected due to zero length sequence
//'
//' @return integer count of fastq entries with zero sequence length
//' @export
// [[Rcpp::export]]
int getZeroLengthSequenceCount()
{
  return(zeroLengthSequence);
}

//' count of fastq elements rejected due to mismatch between sequence and quality field lengths
//'
//' @return integer count of fastq entries with seq/qual length challenges
//' @export
// [[Rcpp::export]]
int getSequenceQualityMismatchCount()
{
  return(sequenceQualityLengthMismatch);
}

//' count of lines of fastq file skipped to enable fastq entry parsing
//'
//' @return integer count of lines skipped
//' @export
// [[Rcpp::export]]
int getSkippedLineCount()
{
  return(linesSkipped);
}

//' parse a fastq file aiming to validate sequences
//'
//' @param x A fastq format DNA/RNA sequence file
//' @export
// [[Rcpp::export]]
LogicalVector fastqValidator(std::string fastq) {

  // reset reported metrics for the given file ...
  //Rcout << "set_fastq_file==" << fastq << std::endl;
  fastq_filename = fastq;
  isZipped = 0;
  running_bases = 0;
  running_fastq = 0;
  malformed_fastq_delim = 0;
  fastq_plus_error = 0;
  zeroLengthSequence = 0;
  sequenceQualityLengthMismatch = 0;
  linesSkipped = 0;
  validFastq = true;

  // TEST (1) - DOES THE SPECIFIED FILE EXIST
  bool exists = myfile_exists(fastq_filename);
  if (!exists) {
    Rcout << "FastqFileNotFound" << std::endl;
    validFastq = false;
    return(LogicalVector(validFastq));
  }

  // HOUSE-KEEPING - IS THE FILE COMPRESSED?
  if (is_gzipped(fastq_filename)==1)
  {
    //Rcout << "provided with a gzip file - flying decompress initiated" << std::endl;
    isZipped = 1;
  }


  if (isZipped == 1) {
    gzfile = gzopen(fastq_filename.c_str(), "r");
  } else {
    file = fopen(fastq_filename.c_str(), "r");
  }
  time(&start_time);

  while (has_next_fastq() == 1)
  {
    get_next_fastq();
  }



  return(LogicalVector(validFastq));
}



