//
// Created by Markus on 9/13/2020.
//

#ifndef BALROG_GENEFINDER_H
#define BALROG_GENEFINDER_H

#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <torch/script.h>

class GeneFinder {
public:
    GeneFinder(torch::jit::script::Module &gene_model, torch::jit::script::Module &TIS_model){
        // models
        _gene_model = gene_model;
        _TIS_model = TIS_model;

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

        // amino acid encoding
        torch::Tensor tensor_L = torch::zeros({1, 21, 1});
        tensor_L[0][1] += 1;
        torch::Tensor tensor_V = torch::zeros({1, 21, 1});
        tensor_V[0][2] += 1;
        torch::Tensor tensor_I = torch::zeros({1, 21, 1});
        tensor_I[0][3] += 1;
        torch::Tensor tensor_M = torch::zeros({1, 21, 1});
        tensor_M[0][4] += 1;
        torch::Tensor tensor_C = torch::zeros({1, 21, 1});
        tensor_C[0][5] += 1;
        torch::Tensor tensor_A = torch::zeros({1, 21, 1});
        tensor_A[0][6] += 1;
        torch::Tensor tensor_G = torch::zeros({1, 21, 1});
        tensor_G[0][7] += 1;
        torch::Tensor tensor_S = torch::zeros({1, 21, 1});
        tensor_S[0][8] += 1;
        torch::Tensor tensor_T = torch::zeros({1, 21, 1});
        tensor_T[0][9] += 1;
        torch::Tensor tensor_P = torch::zeros({1, 21, 1});
        tensor_P[0][10] += 1;
        torch::Tensor tensor_F = torch::zeros({1, 21, 1});
        tensor_F[0][11] += 1;
        torch::Tensor tensor_Y = torch::zeros({1, 21, 1});
        tensor_Y[0][12] += 1;
        torch::Tensor tensor_W = torch::zeros({1, 21, 1});
        tensor_W[0][13] += 1;
        torch::Tensor tensor_E = torch::zeros({1, 21, 1});
        tensor_E[0][14] += 1;
        torch::Tensor tensor_D = torch::zeros({1, 21, 1});
        tensor_D[0][15] += 1;
        torch::Tensor tensor_N = torch::zeros({1, 21, 1});
        tensor_N[0][16] += 1;
        torch::Tensor tensor_Q = torch::zeros({1, 21, 1});
        tensor_Q[0][17] += 1;
        torch::Tensor tensor_K = torch::zeros({1, 21, 1});
        tensor_K[0][18] += 1;
        torch::Tensor tensor_R = torch::zeros({1, 21, 1});
        tensor_R[0][19] += 1;
        torch::Tensor tensor_H = torch::zeros({1, 21, 1});
        tensor_H[0][20] += 1;
        torch::Tensor tensor_null = torch::zeros({1, 21, 1});

        // translation table encoding, 5' to 3' on DNA sense strand
        trantab_standard.insert({"TTT", std::make_pair('F', tensor_F)});
        trantab_standard.insert({"TTC", std::make_pair('F', tensor_F)});
        trantab_standard.insert({"TTA", std::make_pair('L', tensor_L)});
        trantab_standard.insert({"TTG", std::make_pair('L', tensor_L)});
        trantab_standard.insert({"CTT", std::make_pair('L', tensor_L)});
        trantab_standard.insert({"CTC", std::make_pair('L', tensor_L)});
        trantab_standard.insert({"CTA", std::make_pair('L', tensor_L)});
        trantab_standard.insert({"CTG", std::make_pair('L', tensor_L)});
        trantab_standard.insert({"ATT", std::make_pair('I', tensor_I)});
        trantab_standard.insert({"ATC", std::make_pair('I', tensor_I)});
        trantab_standard.insert({"ATA", std::make_pair('I', tensor_I)});
        trantab_standard.insert({"ATG", std::make_pair('M', tensor_M)});
        trantab_standard.insert({"GTT", std::make_pair('V', tensor_V)});
        trantab_standard.insert({"GTC", std::make_pair('V', tensor_V)});
        trantab_standard.insert({"GTA", std::make_pair('V', tensor_V)});
        trantab_standard.insert({"GTG", std::make_pair('V', tensor_V)});
        trantab_standard.insert({"TCT", std::make_pair('S', tensor_S)});
        trantab_standard.insert({"TCC", std::make_pair('S', tensor_S)});
        trantab_standard.insert({"TCA", std::make_pair('S', tensor_S)});
        trantab_standard.insert({"TCG", std::make_pair('S', tensor_S)});
        trantab_standard.insert({"CCT", std::make_pair('P', tensor_P)});
        trantab_standard.insert({"CCC", std::make_pair('P', tensor_P)});
        trantab_standard.insert({"CCA", std::make_pair('P', tensor_P)});
        trantab_standard.insert({"CCG", std::make_pair('P', tensor_P)});
        trantab_standard.insert({"ACT", std::make_pair('T', tensor_T)});
        trantab_standard.insert({"ACC", std::make_pair('T', tensor_T)});
        trantab_standard.insert({"ACA", std::make_pair('T', tensor_T)});
        trantab_standard.insert({"ACG", std::make_pair('T', tensor_T)});
        trantab_standard.insert({"GCT", std::make_pair('A', tensor_A)});
        trantab_standard.insert({"GCC", std::make_pair('A', tensor_A)});
        trantab_standard.insert({"GCA", std::make_pair('A', tensor_A)});
        trantab_standard.insert({"GCG", std::make_pair('A', tensor_A)});
        trantab_standard.insert({"TAT", std::make_pair('Y', tensor_Y)});
        trantab_standard.insert({"TAC", std::make_pair('Y', tensor_Y)});
        trantab_standard.insert({"CAT", std::make_pair('H', tensor_H)});
        trantab_standard.insert({"CAC", std::make_pair('H', tensor_H)});
        trantab_standard.insert({"CAA", std::make_pair('Q', tensor_Q)});
        trantab_standard.insert({"CAG", std::make_pair('Q', tensor_Q)});
        trantab_standard.insert({"AAT", std::make_pair('N', tensor_N)});
        trantab_standard.insert({"AAC", std::make_pair('N', tensor_N)});
        trantab_standard.insert({"AAA", std::make_pair('K', tensor_K)});
        trantab_standard.insert({"AAG", std::make_pair('K', tensor_K)});
        trantab_standard.insert({"GAT", std::make_pair('D', tensor_D)});
        trantab_standard.insert({"GAC", std::make_pair('D', tensor_D)});
        trantab_standard.insert({"GAA", std::make_pair('E', tensor_E)});
        trantab_standard.insert({"GAG", std::make_pair('E', tensor_E)});
        trantab_standard.insert({"TGT", std::make_pair('C', tensor_C)});
        trantab_standard.insert({"TGC", std::make_pair('C', tensor_C)});
        trantab_standard.insert({"TGG", std::make_pair('W', tensor_W)});
        trantab_standard.insert({"CGT", std::make_pair('R', tensor_R)});
        trantab_standard.insert({"CGC", std::make_pair('R', tensor_R)});
        trantab_standard.insert({"CGA", std::make_pair('R', tensor_R)});
        trantab_standard.insert({"CGG", std::make_pair('R', tensor_R)});
        trantab_standard.insert({"AGT", std::make_pair('S', tensor_S)});
        trantab_standard.insert({"AGC", std::make_pair('S', tensor_S)});
        trantab_standard.insert({"AGA", std::make_pair('R', tensor_R)});
        trantab_standard.insert({"AGG", std::make_pair('R', tensor_R)});
        trantab_standard.insert({"GGT", std::make_pair('G', tensor_G)});
        trantab_standard.insert({"GGC", std::make_pair('G', tensor_G)});
        trantab_standard.insert({"GGA", std::make_pair('G', tensor_G)});
        trantab_standard.insert({"GGG", std::make_pair('G', tensor_G)});

        // nucleotide encoding
        torch::Tensor tensor_nucA = torch::zeros({1, 4, 1}); // TODO ensure having a single 0-hot category in a onehot encoded TCN model doesnt impact performance, may help to retrain TIS model with full one-hot encoding
        torch::Tensor tensor_nucT = torch::zeros({1, 4, 1});
        tensor_nucT[0][1] += 1;
        torch::Tensor tensor_nucG = torch::zeros({1, 4, 1});
        tensor_nucG[0][2] += 1;
        torch::Tensor tensor_nucC = torch::zeros({1, 4, 1});
        tensor_nucC[0][3] += 1;

        // nucleotide table, on DNA sense strand
        nuctab.insert({'A', tensor_nucA});
        nuctab.insert({'T', tensor_nucT});
        nuctab.insert({'G', tensor_nucG});
        nuctab.insert({'C', tensor_nucC});

        // hardcoded linear model weights and bias from trained gene model
        gene_model_linear_weights = torch::empty(32);
        int i = 0;
        for (double x : {0.3667307496070862,
                        0.644554615020752,
                        -0.8310428261756897,
                        0.6387593150138855,
                        -0.8136463761329651,
                        0.5082396864891052,
                        0.49443182349205017,
                        -0.5467190146446228,
                        0.6815416216850281,
                        -0.6435322761535645,
                        -0.5412205457687378,
                        0.48874253034591675,
                        0.41401681303977966,
                        0.44753748178482056,
                        0.4844626188278198,
                        0.5019670128822327,
                        0.5385213494300842,
                        -0.45349496603012085,
                        0.16833889484405518,
                        -0.520717442035675,
                        -0.47433432936668396,
                        -0.5133400559425354,
                        -0.9626453518867493,
                        0.5767906308174133,
                        -0.4973398447036743,
                        0.4504794478416443,
                        0.4339512884616852,
                        -0.7323018908500671,
                        0.5616238117218018,
                        0.4348728656768799,
                        0.4719271659851074,
                        -0.4932844340801239}){
            gene_model_linear_weights[i] = x;
            ++i;
        }
        gene_model_linear_bias = 0.19441358745098114;

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
    torch::jit::script::Module _gene_model;
    torch::jit::script::Module _TIS_model;
    torch::Tensor gene_model_linear_weights;
    double gene_model_linear_bias;
    std::string seq;
    int trantab;
    int minimum_length;
    int max_overlap;
    int gene_batch;
    int TIS_batch;

    // functions
    void complement(std::string &s, std::string &comp);

    void find_ORFs();
    void translate_and_tensorize_seq();

    void score_ORFs();

    void get_gene_model_scores();
    void get_TIS_model_scores();
    void get_start_codon_scores();
    void combine_scores();
    void maximize_coherent_score();

    void get_kmer_filter_scores();
    void run_mmseqs();

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

    std::vector<torch::Tensor> seq_tensorized_f0;
    std::vector<torch::Tensor> seq_tensorized_f1;
    std::vector<torch::Tensor> seq_tensorized_f2;
    std::vector<torch::Tensor> seq_tensorized_r0;
    std::vector<torch::Tensor> seq_tensorized_r1;
    std::vector<torch::Tensor> seq_tensorized_r2;

    std::vector<std::pair<int, int>> ORF_coords;
    std::vector<int> ORF_lengths;
    std::vector<std::string> ORF_protein_seq_3to5_nostart_nostop;
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
    std::unordered_map<std::string, std::pair<char, torch::Tensor>> trantab_standard;

    // nucelotide encoding table
    std::unordered_map<char, torch::Tensor> nuctab;

};


#endif //BALROG_GENEFINDER_H
