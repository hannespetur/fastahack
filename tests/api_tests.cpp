#include "../src/FastaHackAPI.hpp"

int main (int argc, char** argv)
{
  std::cout << "Hello world!" << std::endl;
  seqan::String<char> reference_file_name_1 = "/home/hannese/user/git/fastahack/tests/correct.fasta";
  FastaHackAPI reference1 = FastaHackAPI(reference_file_name_1);
  reference1.index();

  std::string region1 = "2:10-30";
  seqan::String<seqan::Dna5> reference_seqan1 = reference1.extract_region(region1);
  std::cout << reference_seqan1 << std::endl;

  seqan::String<char> reference_file_name_2 = "/home/hannese/user/git/fastahack/tests/correct_with_N.fasta";
  FastaHackAPI reference2 = FastaHackAPI(reference_file_name_2);
  reference2.index();

  std::string region2 = "2:8-120";
  seqan::String<seqan::Dna5> reference_seqan2 = reference2.extract_region(region2);
  std::cout << reference_seqan2 << std::endl;

}