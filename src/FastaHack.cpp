#include "Fasta.h"
#include <stdlib.h>
#include <getopt.h>
#include "disorder.h"
#include "Region.h"

void printSummary()
{
  std::cerr << "usage: fastahack [options] <fasta reference>" << std::endl
            << std::endl
            << "options:" << std::endl 
            << "    -i, --index          generate fasta index <fasta reference>.fai" << std::endl
            << "    -r, --region REGION  print the specified region" << std::endl
            << "    -c, --stdin          read a stream of line-delimited region specifiers on stdin" << std::endl
            << "                         and print the corresponding sequence for each on stdout" << std::endl
            << "    -e, --entropy        print the shannon entropy of the specified region" << std::endl
            << "    -d, --dump           print the fasta file in the form 'seq_name <tab> sequence'" << std::endl
            << std::endl
            << "REGION is of the form <seq>, <seq>:<start>[sep]<end>, <seq1>:<start>[sep]<seq2>:<end>" << std::endl
            << "where start and end are 1-based, and the region includes the end position." << std::endl
            << "[sep] is \"-\" or \"..\"" << std::endl
            << std::endl
            << "Specifying a sequence name alone will return the entire sequence, specifying" << std::endl
            << "range will return that range, and specifying a single coordinate pair, e.g." << std::endl
            << "<seq>:<start> will return just that base." << std::endl
            << std::endl
            << "author: Erik Garrison <erik.garrison@bc.edu>" << std::endl;
}


int main (int argc, char** argv)
{
  std::string command;
  std::string fastaFileName;
  std::string seqname;
  std::string longseqname;
  bool dump = false;
  bool buildIndex = false;  // flag to force index building
  bool printEntropy = false;  // entropy printing
  bool readRegionsFromStdin = false;
  std::string region;
  int c;

  while (true)
  {
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"help", no_argument, 0, 'h'},
        {"index",  no_argument, 0, 'i'},
        {"entropy", no_argument, 0, 'e'},
        {"region", required_argument, 0, 'r'},
        {"stdin", no_argument, 0, 'c'},
        {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;
    c = getopt_long (argc, argv, "hciedr:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
          break;
        printf ("option %s", long_options[option_index].name);
        if (optarg)
          printf (" with arg %s", optarg);
        printf ("\n");
        break;

      case 'e':
        printEntropy = true;
        break;

      case 'c':
        readRegionsFromStdin = true;
        break;

      case 'i':
        buildIndex = true;
        break;

      case 'r':
        region = optarg;
        break;

        case 'd':
            dump = true;
            break;

      case 'h':
        printSummary();
        exit(0);
        break;

      case '?':
        /* getopt_long already printed an error message. */
        printSummary();
        exit(1);
        break;

      default:
        abort ();
    }
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
  {
    //cerr << "fasta file: " << argv[optind] << std::endl;
    fastaFileName = argv[optind];
  }
  else
  {
    std::cerr << "Please specify a FASTA file." << std::endl;
    printSummary();
    exit(1);
  }

  if (buildIndex)
  {
    FastaIndex* fai = new FastaIndex();
    //cerr << "generating fasta index file for " << fastaFileName << std::endl;
    fai->indexReference(fastaFileName);
    fai->writeIndexFile((std::string) fastaFileName + fai->indexFileExtension());
  }
  
  std::string sequence;  // holds sequence so we can optionally process it

  FastaReference fr;
  fr.open(fastaFileName);

  if (dump)
  {
    for (vector<std::string>::iterator s = fr.index->sequenceNames.begin(); s != fr.index->sequenceNames.end(); ++s)
    {
      std::cout << *s << "\t" << fr.getSequence(*s) << std::endl;
    }

    return 0;
  }

  if (region != "")
  {
    FastaRegion target(region);
    sequence = fr.getTargetSubSequence(target);
  }

  if (readRegionsFromStdin)
  {
    std::string regionstr;

    while (getline(cin, regionstr))
    {
      FastaRegion target(regionstr);

      if (target.startPos == -1)
      {
        std::cout << fr.getSequence(target.startSeq) << std::endl;
      }
      else
      {
        std::cout << fr.getSubSequence(target.startSeq, target.startPos - 1, target.length()) << std::endl;
      }
    }
  }
  else
  {
    if (sequence != "")
    {
      if (printEntropy)
      {
        if (sequence.size() > 0)
        {
          std::cout << shannon_H((char*) sequence.c_str(), sequence.size()) << std::endl;
        }
        else
        {
          std::cerr << "please specify a region or sequence for which to calculate the shannon entropy" << std::endl;
        }
      }
      else
      {
        // if no statistical processing is requested, just print the sequence
        std::cout << sequence << std::endl;
      }
    }
  }

  return 0;
}
