#include <iostream>
#include <cctype>
#include <cstdio>
#include "cxxopts.hpp"
#include "FastaReader.h"
#include "GeneFinder.h"
#include <cmrc/cmrc.hpp>

CMRC_DECLARE(cmakeresources);

int main(int argc, char* argv[]) {
    // parse command line options
    cxxopts::Options options("Balrog", "Balrog is a prokaryotic gene finder based on a temporal convolutional network");
    options.add_options()
            ("i, in", "Path to input fasta or gzipped fasta", cxxopts::value<std::string>())
            ("o, out", "Path to output annotation", cxxopts::value<std::string>())
            ("max-overlap", "Maximum allowable overlap between genes in nucleotides", cxxopts::value<int>()->default_value("60"))
            ("min-length", "Minimum allowable gene length in nucleotides", cxxopts::value<int>()->default_value("90"))
            ("table", "Nucleotide to amino acid translation table. 11 for most bacteria/archaea, 4 for Mycoplasma/Spiroplasma.",
                    cxxopts::value<int>()->default_value("11"))
            ("max-connections", "Maximum number of forward connections in the directed acyclic graph used to find a set of coherent genes in each genome.",
                    cxxopts::value<int>()->default_value("50"))
            ("gene-batch-size", "Batch size for the temporal convolutional network used to score genes.",
                    cxxopts::value<int>()->default_value("200"))
            ("TIS-batch-size", "Batch size for the temporal convolutional network used to score TIS.",
                    cxxopts::value<int>()->default_value("5000"))
            ("verbose", "Verbose output, set --verbose=false to suppress output text", cxxopts::value<bool>()->default_value("true"))
            ("mmseqs", "Use MMseqs2 to reduce false positive rate, set --mmseqs=false to run without mmseqs", cxxopts::value<bool>()->default_value("true"))
            ("clear-cache", "Force MMseqs2 to remake index, set --clear-cache=true to remake mmseqs reference database and index files", cxxopts::value<bool>()->default_value("false"))
            ("h, help", "Print Balrog usage")
            ;
    auto result = options.parse(argc, argv);

    // check validity and display help
    if (result.count("help") or not result.count("in")){
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // check translation table
    int table = result["table"].as<int>();
    if (table != 11 and table != 4){
        std::cout << "Only translation tables 11 and 4 are currently implemented. Please open a GitHub issue if you need another." << std::endl;
        exit(0);
    }


    // PREPARE MODELS AND DATA
    // open embedded filesystem to load models and data
    auto fs = cmrc::cmakeresources::get_filesystem();

    // get path to temp directory
    std::string tmp_dir;
    char* tmp = std::tmpnam(nullptr);
    int lastslash = ((std::string) tmp).rfind('/');
    tmp_dir = ((std::string) tmp).substr(0, lastslash);
    if (result["verbose"].as<bool>()){
        std::cout << "Saving temp files to " + tmp_dir << std::endl;
    }


    // load LibTorch jit traced gene model
    if (result["verbose"].as<bool>()){
        std::cout << "Importing gene model..." << std::endl;
    }
    // write gene model from embedded filesystem to tmp file
    std::ofstream stream;
    char* tmp_genemodel_path = std::tmpnam(nullptr);
    auto rc = fs.open("gene_model_v1.0.pt");
    stream.open(tmp_genemodel_path);
    auto it = rc.begin();
    while (it != rc.end()){
        stream << *it;
        it ++;
    }
    stream.close();
    // load gene model using jit
    torch::jit::script::Module gene_model;
    try {
        gene_model = torch::jit::load(tmp_genemodel_path);
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the LibTorch gene model\n";
        remove(tmp_genemodel_path);
        return -1;
    }
    // clean up temp file
    remove(tmp_genemodel_path);

    // load LibTorch jit traced TIS model
    if (result["verbose"].as<bool>()){
        std::cout << "Importing TIS model..." << std::endl;
    }
    char* tmp_TISmodel_path = std::tmpnam(nullptr);
    rc = fs.open("TIS_model_v1.0.pt");
    stream.open(tmp_TISmodel_path);
    it = rc.begin();
    while (it != rc.end()){
        stream << *it;
        it ++;
    }
    stream.close();
    // load TIS model using jit
    torch::jit::script::Module TIS_model;
    try {
        TIS_model = torch::jit::load(tmp_TISmodel_path);
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the LibTorch TIS model\n";
        remove(tmp_TISmodel_path);
        return -1;
    }
    // clean up temp file
    remove(tmp_TISmodel_path);

    // load reference gene sequence
    if (result["mmseqs"].as<bool>()){
        std::string ref_fasta_path = tmp_dir + "/reference_genes.fasta";
        std::string ref_db_path = tmp_dir + "/reference_genes.db";
        std::string ref_index_path = tmp_dir + "/balrog_mmseqs_index";

        // use precomputed mmseqs database to save time on multiple runs
        std::ifstream infile(ref_index_path);
        if (result["verbose"].as<bool>() and (not result["clear-cache"].as<bool>())) {
            std::cout << "Found MMseqs2 index at " + ref_index_path << std::endl;
        }

        // create new mmseqs index if needed
        if ((not infile.good()) or (result["clear-cache"].as<bool>())) {
            if (result["verbose"].as<bool>()) {
                std::cout << "Loading reference genes..." << std::endl;
            }
            char *tmp_reference_path = std::tmpnam(nullptr);
            rc = fs.open("reference_genes.tar.gz");
            stream.open(tmp_reference_path);
            it = rc.begin();
            while (it != rc.end()) {
                stream << *it;
                it++;
            }
            stream.close();
            // unzip reference genes
            std::string command = "tar -xf " + (std::string) tmp_reference_path + " -C " + tmp_dir;
            int status = std::system(command.c_str());
            if (status != 0) {
                std::cerr << "error loading reference genes\n";
                return -1;
            }
            // clean temp file
            remove(tmp_reference_path);

            // create mmseqs reference database and index
            command = "mmseqs createdb " + ref_fasta_path + " " + ref_db_path;
            if (not result["verbose"].as<bool>()){
                command += " -v 0";
            }
            status = std::system(command.c_str());
            if (status != 0) {
                std::cerr << "error creating mmseqs database\n";
                return -1;
            }
            command = "mmseqs createindex " + ref_db_path + " " + ref_index_path;
            if (not result["verbose"].as<bool>()){
                command += " -v 0";
            }
            status = std::system(command.c_str());
            if (status != 0) {
                std::cerr << "error creating mmseqs index\n";
                return -1;
            }
        }
    }


    // PREDICT GENES
    // read fasta
    if (result["verbose"].as<bool>()){
        std::cout << "Reading fasta..." << std::endl;
    }
    std::vector<std::string> seq_vec;
    std::vector<std::string> contigname_vec;

    std::string in_path = result["in"].as<std::string>();
    FastaReader io;
    io.read_fasta(in_path, seq_vec, contigname_vec);

    // capitalize all nucleotides
    for (auto &seq : seq_vec){
        for(auto &c: seq){
            c = toupper(c);
        }
    }

    // find genes on all contigs
    std::vector<std::vector<std::pair<int, int>>> gene_coord_vec;
    std::vector<std::vector<bool>> gene_strand_vec;
    std::vector<std::vector<std::string>> gene_nucseq_vec;
    std::vector<std::vector<std::string>> gene_protseq_vec;
    std::vector<std::vector<double>> gene_score_vec;

    int i = 0;
    for (std::string seq : seq_vec){
        ++i;

        GeneFinder gf(gene_model, TIS_model);
        if (result["verbose"].as<bool>()) {
            std::cout << std::endl << "contig " << i << " of " << seq_vec.size() << " : length " << seq.length() << std::endl;
        }

        std::vector<std::pair<int, int>> gene_coord;
        std::vector<bool> gene_strand;
        std::vector<std::string> gene_nucseq;
        std::vector<std::string> gene_protseq;
        std::vector<double> gene_score;

        gf.find_genes(seq,
                      gene_coord,
                      gene_strand,
                      gene_nucseq,
                      gene_protseq,
                      gene_score,
                      table,
                      result["min-length"].as<int>(),
                      result["max-overlap"].as<int>(),
                      result["verbose"].as<bool>(),
                      result["gene-batch-size"].as<int>(),
                      result["TIS-batch-size"].as<int>(),
                      result["mmseqs"].as<bool>());

        gene_coord_vec.emplace_back(gene_coord);
        gene_strand_vec.emplace_back(gene_strand);
        gene_nucseq_vec.emplace_back(gene_nucseq);
        gene_protseq_vec.emplace_back(gene_protseq);
        gene_score_vec.emplace_back(gene_score);
    }


    // OUTPUT PREDICTIONS
    // export genes to gff annotation file
    if (result["verbose"].as<bool>()) {
        std::cout << "Writing predicted genes to " << result["out"].as<std::string>() << std::endl;
    }
    std::ofstream out;
    out.open(result["out"].as<std::string>());
    out << "##gff-version 3" << std::endl;
    std::string contigname;
    for (int j=0; j < contigname_vec.size(); ++j) {
        // extract name up to first space
        contigname = contigname_vec[j].substr(0, contigname_vec[j].find(' '));
        contigname.erase(std::remove(contigname.begin(), contigname.end(), '>'), contigname.end());
        // write all 1-indexed sequence region names and coords
        out << "##sequence-region " << contigname << " " << 1 << " " << seq_vec[j].length() << std::endl;
    }
    for (int j=0; j < contigname_vec.size(); ++j) {
        contigname = contigname_vec[j].substr(0, contigname_vec[j].find(' '));
        contigname.erase(std::remove(contigname.begin(), contigname.end(), '>'), contigname.end());
        // write CDS features
        for (int k=0; k < gene_coord_vec[j].size(); ++k){
            int coord0;
            int coord1;
            if (gene_strand_vec[j][k]){
                coord0 = gene_coord_vec[j][k].first + 1;
                coord1 = gene_coord_vec[j][k].second + 3;
            } else{
                coord0 = gene_coord_vec[j][k].second + 1;
                coord1 = gene_coord_vec[j][k].first + 3;
            }
            out << contigname << "\tbalrog\tCDS\t" << coord0 << "\t" << coord1 << "\t" << "." << "\t";
            if (gene_strand_vec[j][k]){
                out << "+\t";
            } else{
                out << "-\t";
            }
            out << "0\tinference=ab initio prediction:Balrog;product=hypothetical protein" << std::endl;
        }
    }
    out.close();

    if (result["verbose"].as<bool>()) {
        std::cout << "Done...\n" << std::endl;
    }

    return 0;
}
