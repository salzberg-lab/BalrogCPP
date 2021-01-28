//
// Created by Markus on 9/13/2020.
//

#include "FastaReader.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <gzip/compress.hpp>
#include <gzip/config.hpp>
#include <gzip/decompress.hpp>
#include <gzip/utils.hpp>
#include <gzip/version.hpp>

void FastaReader::read_fasta(const std::string &filepath_in, std::vector<std::string> &seq_vec,
                             std::vector<std::string> &contigname_vec) {
    // read input fasta
    std::ifstream ifile(filepath_in);
    if(ifile) {
        // read file buffer into string
        std::stringstream buffer;
        buffer << ifile.rdbuf();
        std::string filedata = buffer.str();

        // check if gzipped
        const char *compressed_pointer = filedata.data();
        std::size_t compressed_size = filedata.size();
        bool gzipped = gzip::is_compressed(compressed_pointer, compressed_size);
        if (gzipped) {
            // decompress and read all data into string
            std::string decompressed_data = gzip::decompress(compressed_pointer, compressed_size);

            // read compressed fasta data
            std::istringstream istring(decompressed_data);
            std::string line;
            while (std::getline(istring, line)) {
                if (line.back() == '\r'){
                    line.pop_back();
                }
                if (line[0] == '>') {
                    seq_vec.emplace_back(std::string());
                    contigname_vec.emplace_back(line);

                } else {
                    seq_vec.back() += line;
                }
            }
        } else {
            // read text fasta data
            std::istringstream istring(filedata);
            std::string line;
            while (std::getline(istring, line)) {
                if (line.back() == '\r'){
                    line.pop_back();
                }
                if (line[0] == '>') {
                    seq_vec.emplace_back(std::string());
                    contigname_vec.emplace_back(line);

                } else {
                    seq_vec.back() += line;
                }
            }
        }
    } else {
        std::cout << "Error opening file" << std::endl;
    }
}
