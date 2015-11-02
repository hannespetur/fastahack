#include "FastaHackAPI.hpp"

FastaHackAPI::FastaHackAPI(seqan::String<char> reference_fasta_seqan)
{
  reference_fasta_file_name = static_cast<std::string>(seqan::toCString(reference_fasta_seqan));
  fasta_reference.open(reference_fasta_file_name);
}

void FastaHackAPI::index()
{
  FastaIndex* fai = new FastaIndex();
  // std::cerr << "generating fasta index file for " << reference_fasta_file_name << std::endl;
  fai->indexReference(reference_fasta_file_name);
  fai->writeIndexFile((std::string) reference_fasta_file_name + fai->indexFileExtension());
}

seqan::String<seqan::Dna5> FastaHackAPI::extract_region(std::string region)
{
  FastaRegion target(region);
  std::string sequence = fasta_reference.getTargetSubSequence(target);
  char* sequence_cstr = strdup(sequence.c_str());
  seqan::String<seqan::Dna5> seqan_sequence = sequence_cstr;
  return seqan_sequence;
}

seqan::String<seqan::Dna5> FastaHackAPI::extract_region(seqan::String<char> region)
{
  std::string region_str(toCString(region));
  return extract_region(region_str);
}
