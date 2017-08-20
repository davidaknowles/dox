// Compile me with
// g++  -std=c++11  pileup_filter_fast.cpp -o pileup_filter

#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>


template<typename KeyType, typename ValueType> 
std::pair<KeyType,ValueType> get_max( std::map<KeyType,ValueType>& x ) {
  using pairtype=std::pair<KeyType,ValueType>; 
  return *std::max_element(x.begin(), x.end(), [] (const pairtype & p1, const pairtype & p2) {
        return p1.second < p2.second;
  }); 
}

int main (int argc, char *argv[]) {
	int ref_filter=0, alt_filter=0, total_filter=0; 
   if (argc < 2) {
     std::cerr << "pileup_filter usage: samtools mpileup my.bam | pileup_filter <ref_count_threshold> <alt_count_theshold> > my.pileup" << std::endl; 
     return 0;
   }
   if (argc==3) {
	   ref_filter = atoi( argv[1] ); 
	   alt_filter = atoi( argv[2] ); 
	   total_filter = ref_filter + alt_filter;
   } else {
   	   total_filter = atoi( argv[1] ); 
   }

   int pos; 
   std::string pos_string, coverage_string; 
   std::string chr, consensus, reads, quality; 
   long line_count=0, output_lines=0, print_every=10000000;
   std::string ref_string=".,"; 
   std::string alt_string="atcgATCG";
   std::vector<char> bases={ 'A', 'T', 'G', 'C' };
   int ref_count;
   std::map< char, int > alt_counts; 
   while ( std::cin ) { 
     if (getline( std::cin, chr, '\t' )) {
        getline( std::cin, pos_string, '\t' );
        pos=atoi(pos_string.c_str());
        getline( std::cin, consensus, '\t' );
        getline( std::cin, coverage_string, '\t' );
        getline( std::cin, reads, '\t' );
        getline( std::cin, quality );
        if (reads.size() >= total_filter) {
          ref_count=0; 
          for (char base : bases) 
            alt_counts[base] = 0; 
          for(char& c : reads) 
            switch (c) {
            case '.' :
            case ',' : ref_count++; break; 
            case 'A' :
            case 'a' : alt_counts['A']++; break; 
            case 'T' :
            case 't' : alt_counts['T']++; break; 
            case 'G' :
            case 'g' : alt_counts['G']++; break; 
            case 'C' :
            case 'c' : alt_counts['C']++; break; 
          }
          auto max_kvp = get_max(alt_counts);
          int alt_count=max_kvp.second; 
          if (ref_count >= ref_filter & alt_count >= alt_filter & (alt_count+ref_count) >= total_filter) {
            fprintf( stdout, "%s\t%i\t%i\t%i\t%s\t%c\n", chr.c_str(), pos, ref_count, alt_count, consensus.c_str(), max_kvp.first ); 
            output_lines++; 
          }
        }
        line_count++; 
        if ( line_count % print_every == 0 ) 
          std::cerr << "Read " << line_count << " lines, output " << output_lines << std::endl;
      }
   }
   
   return 0; 
}

