#include <Rcpp.h>
using namespace Rcpp;
#include<Rinternals.h>
#include <string.h>
using namespace std;
#include <zlib.h>
#include <stdio.h>

static std::string fastq_filename;
static std::string dest_filename;
static int isZipped = 0;
static int isDestZipped = 0;
static gzFile gzfile, gzDestfile;
static FILE *file, *destfile;
static bool validFastq = true;
static long int running_bases = 0;
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
    //Rcout << "Malformed fastq entry delimitter " << headdelim << "!=@" << std::endl;
    malformed_fastq_delim ++;
    validFastq = false;
    return(realign_fastq());
  }

  // is delim == "+"? - * is used to indicate runover end of file ...
  if (strcmp(fq.delim, "+")!=0)
  {
    //Rcout << "Malformed fastq entry ~ Line3[+] not present" << std::endl;
    fastq_plus_error ++;
    validFastq = false;
    return(0); // this is not a valid fastq
  }

  // is sequence length > 0
  size_t seqlength = strlen(fq.sequence);
  size_t qlength = strlen(fq.quality);
  if (seqlength == 0) {
    //Rcout << "Malformed fastq entry ~ length(sequence)==0" << std::endl;
    zeroLengthSequence ++;
    validFastq = false;
    return(0); // this is not a valid fastq
  }

  // Rcout << "(s/q)Lengths==" << seqlength << "/" << qlength << std::endl;

  // does sequence length == quality length
  if (seqlength != qlength) {
    //Rcout << "Malformed fastq entry ~ length(sequence)!=length(quality)" << std::endl;
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
long int getFastqBases()
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


void reset() {
  isZipped = 0;
  isDestZipped = 0;
  running_bases = 0;
  running_fastq = 0;
  malformed_fastq_delim = 0;
  fastq_plus_error = 0;
  zeroLengthSequence = 0;
  sequenceQualityLengthMismatch = 0;
  linesSkipped = 0;
  validFastq = true;
}



//' parse a fastq file aiming to validate sequences
//'
//' fastqValidator will parse the specified fastq (or fastq.gz) file looking for fastq
//' entry compliance. A boolean value of overall file compliance will be returned. Additional
//' summary counts describing the reason for rejection are also made available
//'
//' @seealso \code{\link{getFastqPlusErrorCount}}, \code{\link{getMalformedFastqHeaderCount}}, \code{\link{getZeroLengthSequenceCount}}, \code{\link{getSequenceQualityMismatchCount}} and \code{\link{getSkippedLineCount}}
//'
//' @param fastq A fastq format DNA/RNA sequence file
//' @return logical defining if fastq provided is indeed valid fastq
//' @export
// [[Rcpp::export]]
LogicalVector fastqValidator(std::string fastq) {

  // reset reported metrics for the given file ...
  //Rcout << "set_fastq_file==" << fastq << std::endl;
  fastq_filename = fastq;
  reset();

  // TEST (1) - DOES THE SPECIFIED FILE EXIST
  if (!myfile_exists(fastq_filename))
  {
    Rcout << "FastqFileNotFound" << std::endl;
    validFastq = false;
    return(LogicalVector(validFastq));
  }

  // HOUSE-KEEPING - IS THE FILE COMPRESSED?
  if (is_gzipped(fastq_filename)==1)
  {
    isZipped = 1;
    gzfile = gzopen(fastq_filename.c_str(), "r");
  } else {
    file = fopen(fastq_filename.c_str(), "r");
  }

  while (has_next_fastq() == 1)
  {
    get_next_fastq();
  }

  if (isZipped) {
    gzclose(gzDestfile);
  } else {
    fclose(file);
  }

  return(LogicalVector(validFastq));
}



char* getFastqEntry()
{
  size_t len1 = strlen(fq.header);
  size_t len2 = strlen(fq.sequence);
  size_t len3 = strlen(fq.delim);
  size_t len4 = strlen(fq.quality);
  //char* fastqE = malloc(len1 + len2 + len3 + len4 + 8);
  char* fastqE = new char[len1 + len2 + len3 + len4 + 8] ;

  fastqE = strcpy(fastqE, fq.header);
  fastqE = strcat(fastqE, "\n");
  fastqE = strcat(fastqE, fq.sequence);
  fastqE = strcat(fastqE, "\n");
  fastqE = strcat(fastqE, fq.delim);
  fastqE = strcat(fastqE, "\n");
  fastqE = strcat(fastqE, fq.quality);
  fastqE = strcat(fastqE, "\n");

  return(fastqE);
}



//' fix a corrupted fastq file (if fastq-like)
//'
//' fixFastq parses a fastq (or fastq.gz) file for compliant fastq entries and writes these to
//' the file specified in newfastq parameter. Any non-compliant reads are dropped
//'
//' @param fastq file location of fastq source
//' @param newfastq location of file to write content to
//' @return path to new fastq file (same as newfastq provided)
//'
//' @export
// [[Rcpp::export]]
std::string fixFastq(std::string fastq, std::string newfastq)
{
  fastq_filename = fastq;
  dest_filename = newfastq;

  reset();
  if (myfile_exists(fastq_filename))
  {
    if (is_gzipped(fastq_filename)==1)
    {
      isZipped = 1;
      gzfile = gzopen(fastq_filename.c_str(), "r");
    } else {
      file = fopen(fastq_filename.c_str(), "r");
    }

    if (is_gzipped(dest_filename)==1)
    {
      isDestZipped = 1;
      gzDestfile = gzopen(dest_filename.c_str(), "wb");
    } else {
      destfile = fopen(dest_filename.c_str(), "w");
    }


    while (has_next_fastq() == 1)
    {
      if (get_next_fastq())
      {
        if (isDestZipped)
        {
          gzputs (gzDestfile, getFastqEntry());
        }
        else
        {
          fputs(getFastqEntry(), destfile);
        }
      }
    }


    if (is_gzipped(dest_filename)==1) {
      Rcout << "closing gz dest" << std::endl;
      gzclose(gzDestfile);
    }
    else
    {
      fclose(destfile);
    }


    if (isZipped) {
      gzclose(gzfile);
    } else {
      fclose(file);
    }


  }
  return(dest_filename);
}
