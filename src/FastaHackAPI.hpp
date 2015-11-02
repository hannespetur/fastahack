/**
 *  This file provides a SeqAn API calls to fastahack.
 */

#ifndef __FASTA_HACK_API_H__
#define __FASTA_HACK_API_H__
#include "Fasta.h"
#include <stdlib.h>
#include "disorder.h"
#include "Region.h"

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <iostream>

class FastaHackAPI
{
 private:
  std::string reference_fasta_file_name;
  FastaReference fasta_reference;


 public:
  FastaHackAPI(seqan::String<char> reference_fasta_seqan);

  void index();

  seqan::String<seqan::Dna5> extract_region(std::string region);

  seqan::String<seqan::Dna5> extract_region(seqan::String<char> region);
};


#endif
