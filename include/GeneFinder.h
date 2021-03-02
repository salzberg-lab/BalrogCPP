//
// Created by Markus on 9/13/2020.
//

#ifndef BALROG_GENEFINDER_H
#define BALROG_GENEFINDER_H

#include <utility>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <random>

class GeneFinder {
public:
    explicit GeneFinder(std::string path){
        // temporary directory for all tmp files
        tmp_dir = std::move(path);
        if (tmp_dir.back() != '/'){
            tmp_dir += "/";
        }

        // ORF score parameters
        float max_start_codon = 0;
        if (weight_ATG > max_start_codon){
            max_start_codon = weight_ATG;
        }
        if (weight_TTG > max_start_codon){
            max_start_codon = weight_TTG;
        }
        if (weight_GTG > max_start_codon){
            max_start_codon = weight_GTG;
        }
        float maxprob = weight_gene_prob + weight_TIS_prob + max_start_codon;
        probthresh = score_threshold * maxprob;

        // translation table encoding, 5' to 3' on DNA sense strand
        // NOTE: lazily converted from LibTorch tensor pair encoding, TODO translation can be much faster if needed
        trantab_standard.insert({"TTT", 'F'});
        trantab_standard.insert({"TTC", 'F'});
        trantab_standard.insert({"TTA", 'L'});
        trantab_standard.insert({"TTG", 'L'});
        trantab_standard.insert({"CTT", 'L'});
        trantab_standard.insert({"CTC", 'L'});
        trantab_standard.insert({"CTA", 'L'});
        trantab_standard.insert({"CTG", 'L'});
        trantab_standard.insert({"ATT", 'I'});
        trantab_standard.insert({"ATC", 'I'});
        trantab_standard.insert({"ATA", 'I'});
        trantab_standard.insert({"ATG", 'M'});
        trantab_standard.insert({"GTT", 'V'});
        trantab_standard.insert({"GTC", 'V'});
        trantab_standard.insert({"GTA", 'V'});
        trantab_standard.insert({"GTG", 'V'});
        trantab_standard.insert({"TCT", 'S'});
        trantab_standard.insert({"TCC", 'S'});
        trantab_standard.insert({"TCA", 'S'});
        trantab_standard.insert({"TCG", 'S'});
        trantab_standard.insert({"CCT", 'P'});
        trantab_standard.insert({"CCC", 'P'});
        trantab_standard.insert({"CCA", 'P'});
        trantab_standard.insert({"CCG", 'P'});
        trantab_standard.insert({"ACT", 'T'});
        trantab_standard.insert({"ACC", 'T'});
        trantab_standard.insert({"ACA", 'T'});
        trantab_standard.insert({"ACG", 'T'});
        trantab_standard.insert({"GCT", 'A'});
        trantab_standard.insert({"GCC", 'A'});
        trantab_standard.insert({"GCA", 'A'});
        trantab_standard.insert({"GCG", 'A'});
        trantab_standard.insert({"TAT", 'Y'});
        trantab_standard.insert({"TAC", 'Y'});
        trantab_standard.insert({"CAT", 'H'});
        trantab_standard.insert({"CAC", 'H'});
        trantab_standard.insert({"CAA", 'Q'});
        trantab_standard.insert({"CAG", 'Q'});
        trantab_standard.insert({"AAT", 'N'});
        trantab_standard.insert({"AAC", 'N'});
        trantab_standard.insert({"AAA", 'K'});
        trantab_standard.insert({"AAG", 'K'});
        trantab_standard.insert({"GAT", 'D'});
        trantab_standard.insert({"GAC", 'D'});
        trantab_standard.insert({"GAA", 'E'});
        trantab_standard.insert({"GAG", 'E'});
        trantab_standard.insert({"TGT", 'C'});
        trantab_standard.insert({"TGC", 'C'});
        trantab_standard.insert({"TGG", 'W'});
        trantab_standard.insert({"CGT", 'R'});
        trantab_standard.insert({"CGC", 'R'});
        trantab_standard.insert({"CGA", 'R'});
        trantab_standard.insert({"CGG", 'R'});
        trantab_standard.insert({"AGT", 'S'});
        trantab_standard.insert({"AGC", 'S'});
        trantab_standard.insert({"AGA", 'R'});
        trantab_standard.insert({"AGG", 'R'});
        trantab_standard.insert({"GGT", 'G'});
        trantab_standard.insert({"GGC", 'G'});
        trantab_standard.insert({"GGA", 'G'});
        trantab_standard.insert({"GGG", 'G'});
    }
    ~GeneFinder() = default;

    void find_genes(std::string &contig_seq,
            std::vector<std::pair<int, int>> &gene_coord,
            std::vector<bool> &gene_strand,
            std::vector<std::string> &gene_nucseq,
            std::vector<std::string> &gene_protseq,
            std::vector<double> &gene_score,
            int table,
            int min_ORF_length,
            int max_ORF_overlap,
            bool verbose,
            int gene_batch_size,
            int TIS_batch_size,
            bool mmseqs);
private:
    // general parameters
    std::string tmp_dir;
    bool verbosity;
    const int k_seengene = 10;
    const int multimer_threshold = 2;
    const int max_forward_connections = 1000;
    const int cutoff = 100;

    // previously optimized parameters
    const float weight_gene_prob = 0.9746869839852076;
    const float weight_TIS_prob = 0.25380288790532707;
    const float weight_ATG = 0.84249804151264;
    const float weight_TTG = 0.7512400826652517;
    const float weight_GTG = 0.7083689705744909;
    const float score_threshold = 0.47256101519707244;
    float probthresh;
    const float unidirectional_penalty_per_base = 3.895921717182765; // 3' 5' overlap
    const float convergent_penalty_per_base = 4.603432608883688; // 3' 3' overlap
    const float divergent_penalty_per_base = 3.3830814940689975; // 5' 5' overlap

    // prelims
    std::string seq;
    int trantab;
    int minimum_length;
    int max_overlap;
    int gene_batch;
    int TIS_batch;

    // functions
    void complement(std::string &s, std::string &comp);

    void find_ORFs();
    void translate_seq();

    void score_ORFs();

    void get_gene_model_scores();
    void get_TIS_model_scores();
    void get_start_codon_scores();
    void combine_scores();
    void maximize_coherent_score();

    void get_kmer_filter_scores();
    void run_mmseqs();

    int mt_rand();

    // gene info
    std::vector<std::pair<int, int>> *_gene_coord;
    std::vector<bool> *_gene_strand;
    std::vector<std::string> *_gene_nucseq;
    std::vector<std::string> *_gene_protseq;
    std::vector<double> *_gene_score;

    // seq information
    std::string seq_translated_f0;
    std::string seq_translated_f1;
    std::string seq_translated_f2;
    std::string seq_translated_r0;
    std::string seq_translated_r1;
    std::string seq_translated_r2;

    std::vector<std::pair<int, int>> ORF_coords;
    std::vector<int> ORF_lengths;
    std::vector<std::vector<char>> ORF_protein_seq_3to5_nostart_nostop;
    std::vector<std::string> ORF_start_codon;

    // calculated values
    std::vector<double> ORF_scores;
    std::vector<double> gene_model_scores;
    std::vector<double> TIS_model_scores;
    std::vector<double> start_codon_scores;
    std::vector<int> length_scores;
    std::vector<bool> kmer_filter_scores;
    std::vector<double> MMSeqs2_scores;

    // translation table
    std::unordered_map<std::string, char> trantab_standard;

};


#endif //BALROG_GENEFINDER_H
