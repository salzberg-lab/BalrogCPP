//
// Created by Markus on 9/13/2020.
//

#include <vector>
#include <string>

#ifndef BALROG_FASTAREADER_H
#define BALROG_FASTAREADER_H


class FastaReader {
public:
    FastaReader() = default;
    ~FastaReader() = default;

    void read_fasta(const std::string &filepath_in, std::vector<std::string> &seq_vec, std::vector<std::string> &contigname_vec);
};


#endif //BALROG_FASTAREADER_H
