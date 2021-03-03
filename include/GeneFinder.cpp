//
// Created by Markus on 9/13/2020.
//

#include <iostream>
#include <algorithm>
#include "GeneFinder.h"
#include "tqdm.h"
#include <set>
#include <fstream>
#include <iterator>
#include <sstream>

void GeneFinder::find_genes(std::string &contig_seq,
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
        bool mmseqs) {

    _gene_coord = &gene_coord;
    _gene_strand = &gene_strand;
    _gene_nucseq = &gene_nucseq;
    _gene_protseq = &gene_protseq;
    _gene_score = &gene_score;


    verbosity = verbose;
    trantab = table;
    minimum_length = min_ORF_length;
    seq = contig_seq;
    gene_batch = gene_batch_size;
    TIS_batch = TIS_batch_size;
    max_overlap = max_ORF_overlap;

    if (verbosity){
        std::cout << "Finding open reading frames..." << std::endl;
    }
    find_ORFs();
    if (verbosity){
        std::cout << "Translating open reading frames..." << std::endl;
    }
    translate_seq();
    score_ORFs();
    if (verbosity){
        std::cout << "Solving for geneiest path through contig..." << std::endl;
    }
    maximize_coherent_score();

    if (mmseqs){
        if (verbosity) {
            std::cout << "Checking low-scoring ORFs against reference genes with MMseqs2..." << std::endl;
        }
        run_mmseqs();
    }

//    for (int j=0; j < ORF_scores.size(); ++j){
////        int i = ORF_scores.size() - j - 1;
//        int i = j;
//        std::cout <<
//                  ORF_coords[i].first << " " <<
//                  ORF_coords[i].second << " " <<
//                  (ORF_coords[i].first < ORF_coords[i].second) << " " <<
//                  ORF_lengths[i] << " " <<
//                  ORF_scores[i] << " " <<
//                  gene_model_scores[i] << " ";
//        for (char k : ORF_protein_seq_3to5_nostart_nostop[i]){
//            std::cout << k;
//        }
//        std::cout << std::endl;
//    }
//    std::cout << ORF_protein_seq_3to5_nostart_nostop.size() << std::endl;
//    std::cout << ORF_scores.size() << std::endl;

//    for (auto & i : ORF_protein_seq_3to5_nostart_nostop){
//        for (char c : i){
//            std::cout << c;
//        }
//        std::cout << std::endl;
//    }
//    for (int i=0; i < ORF_protein_seq_3to5_nostart_nostop.size(); i++){
//        int j = ORF_protein_seq_3to5_nostart_nostop.size() - i - 1;
//        std::cout << ORF_protein_seq_3to5_nostart_nostop[j].size();
//        std::cout << std::endl;
//    }
}

void GeneFinder::run_mmseqs(){
    // locate tmp dir
    std::string tmp_dir;
    char* tmp = std::tmpnam(nullptr);
    int lastslash = ((std::string) tmp).rfind('/');
    tmp_dir = ((std::string) tmp).substr(0, lastslash);

    // find low-scoring ORFs
    std::vector<int> low_score_idx;
    for (int i=0; i < (*_gene_score).size(); ++i){
        if ((*_gene_score)[i] < cutoff){
            low_score_idx.emplace_back(i);
        }
    }

    // stop if no low-scorers
    if (low_score_idx.empty()){
        return;
    }

    // write low-scoring ORFs to temp file
    std::ofstream stream;
    char* tmp_lowscore_path = std::tmpnam(nullptr);
    stream.open(tmp_lowscore_path);
    for (int i=0; i < low_score_idx.size(); ++i){
        stream << ">" << i << "\n";
        stream << (*_gene_protseq)[low_score_idx[i]] << "\n";
    }
    stream.close();


    // run mmseqs
    char tmp_alignment_path[L_tmpnam];
    std::tmpnam(tmp_alignment_path);

    std::string ref_db_path = tmp_dir + "/reference_genes.db"; // TODO link these to main definition
    std::string ref_index_path = tmp_dir + "/balrog_mmseqs_index";

    std::string mmseqs_sensitivity = "7.0";
    std::string command = "mmseqs easy-search -s " + mmseqs_sensitivity + " " \
    + (std::string) tmp_lowscore_path + " " \
    + (std::string) ref_db_path + " " \
    + (std::string) tmp_alignment_path + " " \
    + (std::string) ref_index_path;

    if (verbosity){
        command += " -v 1";
    }
    else {
        command += " -v 0";
    }

    int status = std::system(command.c_str());
    if (status != 0) {
        std::cerr << "error running mmseqs, try --clear-cache\n";
        exit(status);
    }

    // read mmseqs output
    std::ifstream instream;
    instream.open(tmp_alignment_path);
    std::set<int> hitset;
    for(std::string line; getline(instream, line);){
        int tab = ((std::string) line).find('\t');
        int hitidx = stoi(((std::string) line).substr(0, tab));
        hitset.insert(hitidx);
    }
    instream.close();

    // remove all low-scoring ORFs without mmseqs hits
    // remove hits backward from highest idx
    std::vector<int> hitvec;
    hitvec.reserve(hitset.size());
    for (int h : hitset){
        hitvec.push_back(h);
    }
    std::reverse(hitvec.begin(), hitvec.end()); // TODO this can be replaced by reverse iterating over the sorted set

    for (int idx : hitvec){
        low_score_idx.erase(low_score_idx.begin() + idx);
    }

    // remove non-hit low-scores backward from highest idx
    std::reverse(low_score_idx.begin(), low_score_idx.end());
    for (int idx : low_score_idx){
        _gene_coord->erase(_gene_coord->begin() + idx);
        _gene_strand->erase(_gene_strand->begin() + idx);
        _gene_nucseq->erase(_gene_nucseq->begin() + idx);
        _gene_protseq->erase(_gene_protseq->begin() + idx);
        _gene_score->erase(_gene_score->begin() + idx);
    }

}

void GeneFinder::score_ORFs() {
    if (verbosity){
        std::cout << "Scoring longest open reading frames with gene model..." << std::endl;
    }
    get_gene_model_scores();
    if (verbosity){
        std::cout << "Scoring potential start sites with TIS model..." << std::endl;
    }
    get_TIS_model_scores();
    get_start_codon_scores();
    combine_scores();
}

bool toposort_ORFs(const std::pair<std::pair<int, int>, int>& a,
                   const std::pair<std::pair<int, int>, int>& b){
    return (a.first.first < b.first.first);
}

bool threeprimefiveprimesort_ORFs(const std::pair<std::pair<int, int>, int>& a,
                                  const std::pair<std::pair<int, int>, int>& b){
    // returns sorting by stop codon, with longest ORF first if sharing stop codon
    // assumes a and b are on same strand
    // sort by stop codon first
    if (a.first.second != b.first.second){
        return (a.first.second < b.first.second);
    }
    // if stop codon is the same, sort by start codon depending on strand
    else if (a.first.first < a.first.second){
        // forward strand
        return (a.first.first < b.first.first);
    } else{
        // reverse strand
        return (a.first.first > b.first.first);
    }
}

bool ORFsort_vector(const std::pair<int, std::vector<float>>& a,
                    const std::pair<int, std::vector<float>>& b){
    return (a.first < b.first);
}

//bool ORFsort_tensor(const std::pair<int, torch::Tensor>& a,
//                    const std::pair<int, torch::Tensor>& b){
//    return (a.first < b.first);
//}

// LIBTORCH RELIANT
//void GeneFinder::get_gene_model_scores() {
//    // get lengths for later scoring
//    for (auto & ORF_coord : ORF_coords) {
//        int l = std::abs(ORF_coord.first - ORF_coord.second);
//        ORF_lengths.emplace_back(l);
//    }
//
//    // separate strands
//    std::vector<std::pair<std::pair<int, int>, int>> ORFs_forward;
//    std::vector<std::pair<std::pair<int, int>, int>> ORFs_reverse;
//
//    for (int i=0; i < ORF_coords.size(); ++i) {
//        std::pair<std::pair<int, int>, int> coord_idx = std::make_pair(ORF_coords[i], i);
//        if (coord_idx.first.first < coord_idx.first.second){
//            ORFs_forward.emplace_back(coord_idx);
//        } else{
//            ORFs_reverse.emplace_back(coord_idx);
//        }
//    }
//
//    // sort by stop codon, with ORFs sorted by length within shared stop codon
//    std::sort(ORFs_forward.begin(), ORFs_forward.end(), threeprimefiveprimesort_ORFs);
//    std::sort(ORFs_reverse.begin(), ORFs_reverse.end(), threeprimefiveprimesort_ORFs);
//
//    // extract longest ORF per stop codon
//    std::vector<std::pair<int, int>> ORF_coords_longest;
//
//    // forward
//    int scoredstop = -1;
//    for (auto & O : ORFs_forward){
//        auto coords = O.first;
//        int stop = coords.second;
//        if (scoredstop == stop){
//            continue;
//        }
//        scoredstop = stop;
//
//        ORF_coords_longest.emplace_back(coords);
//    }
//    // reverse
//    scoredstop = -1;
//    for (auto & O : ORFs_reverse){
//        auto coords = O.first;
//        int stop = coords.second;
//        if (scoredstop == stop){
//            continue;
//        }
//        scoredstop = stop;
//
//        ORF_coords_longest.emplace_back(coords);
//    }
//
//
//    // sort by length, keep track of original index
//    std::vector<std::tuple<int, std::pair<int, int>, int>> ORF_lengths_coords_sort;
//    for (int i=0; i < ORF_coords_longest.size(); ++i){
//        int l = std::abs(ORF_coords_longest[i].first-ORF_coords_longest[i].second);
//        ORF_lengths_coords_sort.emplace_back(std::make_tuple(l, ORF_coords_longest[i], i));
//    }
//    std::sort(ORF_lengths_coords_sort.begin(), ORF_lengths_coords_sort.end());
//
//    // score in length-sorted blocks
//    std::vector<std::pair<int, torch::Tensor>> ORF_coords_longest_tensorscore;
//    tqdm progressbar;
//    progressbar.set_theme_braille_spin();
//    for (int i=0; i < ORF_lengths_coords_sort.size(); i += gene_batch){
//        std::vector<torch::Tensor> ORF_block;
//        int block_len;
//        if (i + gene_batch <= ORF_lengths_coords_sort.size()){
//            block_len = std::get<0>(ORF_lengths_coords_sort[i+gene_batch]) / 3;
//        } else{
//            // last block may be of different size
//            block_len = std::get<0>(ORF_lengths_coords_sort.back()) / 3;
//        }
//        // create length-sorted tensor blocks
//        unsigned long block_end = std::min((unsigned long)(i + gene_batch), ORF_lengths_coords_sort.size());
//
//        for (int j=i; j < block_end; ++j){
//            if (verbosity) {
//                progressbar.progress(j, (int) ORF_lengths_coords_sort.size());
//            }
//            // extract amino acid sequence in correct reading frame and add it to tensor block
//            std::pair<int, int> coords = std::get<1>(ORF_lengths_coords_sort[j]);
//            // forward strand
//            if (coords.first < coords.second) {
//                if (coords.second%3 == 0) {
//                    auto begin = seq_tensorized_f0.begin() + coords.first / 3 + 1;
//                    auto end = begin + (coords.second - coords.first) / 3 - 1;
//                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
//                    std::reverse(vec.begin(), vec.end()); // reverse so sequence is 3' to 5'
//                    unsigned long pad = block_len - vec.size(); // pad so all tensors are the same size in each block
//                    for (int k=0; k < pad; ++k){
//                        vec.emplace_back(torch::zeros({1, 21,1}));
//                    }
//                    ORF_block.emplace_back(torch::cat(vec, -1));
//                } else if (coords.second%3 == 1) {
//                    auto begin = seq_tensorized_f1.begin() + coords.first / 3 + 1;
//                    auto end = begin + (coords.second - coords.first) / 3 - 1;
//                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
//                    std::reverse(vec.begin(), vec.end());
//                    unsigned long pad = block_len - vec.size();
//                    for (int k=0; k < pad; ++k){
//                        vec.emplace_back(torch::zeros({1, 21,1}));
//                    }
//                    ORF_block.emplace_back(torch::cat(vec, -1));
//                } else {
//                    auto begin = seq_tensorized_f2.begin() + coords.first / 3 + 1;
//                    auto end = begin + (coords.second - coords.first) / 3 - 1;
//                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
//                    std::reverse(vec.begin(), vec.end());
//                    unsigned long pad = block_len - vec.size();
//                    for (int k=0; k < pad; ++k){
//                        vec.emplace_back(torch::zeros({1, 21,1}));
//                    }
//                    ORF_block.emplace_back(torch::cat(vec, -1));
//                }
//            }
//            // reverse strand
//            else {
//                if (coords.first%3 == 0) {
//                    auto begin = seq_tensorized_r0.begin() + coords.second / 3 + 1;
//                    auto end = begin + (coords.first - coords.second) / 3 - 1;
//                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
//                    unsigned long pad = block_len - vec.size();
//                    for (int k=0; k < pad; ++k){
//                        vec.emplace_back(torch::zeros({1, 21,1}));
//                    }
//                    ORF_block.emplace_back(torch::cat(vec, -1));
//                } else if (coords.first%3 == 1){
//                    auto begin = seq_tensorized_r1.begin() + coords.second / 3 + 1;
//                    auto end = begin + (coords.first - coords.second) / 3 - 1;
//                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
//                    unsigned long pad = block_len - vec.size();
//                    for (int k=0; k < pad; ++k){
//                        vec.emplace_back(torch::zeros({1, 21,1}));
//                    }
//                    ORF_block.emplace_back(torch::cat(vec, -1));
//                } else {
//                    auto begin = seq_tensorized_r2.begin() + coords.second / 3 + 1;
//                    auto end = begin + (coords.first - coords.second) / 3 - 1;
//                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
//                    unsigned long pad = block_len - vec.size();
//                    for (int k=0; k < pad; ++k){
//                        vec.emplace_back(torch::zeros({1, 21,1}));
//                    }
//                    ORF_block.emplace_back(torch::cat(vec, -1));
//                }
//            }
//        }
//        // score each batch block
//        if (not ORF_block.empty()) {
//            // get values for last hidden unit layer {batch_size, hidden_units, length}
//            torch::Tensor ORF_batch = torch::cat({ORF_block}, 0);
//            torch::NoGradGuard nograd;
//            _gene_model.eval();
//            auto score = _gene_model.forward({ORF_batch});
//            torch::Tensor score_tensor = score.toTensor();
//
//            // apply linear function across all causal hidden units (equivalent to getting score at each character in the sequence)
//            std::vector<torch::Tensor> linear_layer_score_vec;
//            for (int j=0; j < score_tensor.sizes()[2]; ++j){
//                torch::Tensor linear_layer_score = torch::matmul(score_tensor.index({torch::indexing::Slice(), torch::indexing::Slice(), j}), gene_model_linear_weights) + gene_model_linear_bias;
//                linear_layer_score_vec.emplace_back(linear_layer_score);
//            }
//
//            // combine post-linear-layer scores into {batch_size, length} matrix
//            torch::Tensor linear_layer_matrix = torch::stack({linear_layer_score_vec}, 1);
//
//            int k = 0;
//            for (int j=i; j < block_end; ++j) {
//                torch::Tensor logit_all = linear_layer_matrix.index({k, torch::indexing::Slice(torch::indexing::None, torch::indexing::None)});
//                int original_idx = std::get<2>(ORF_lengths_coords_sort[j]);
//                ORF_coords_longest_tensorscore.emplace_back(std::make_pair(original_idx, logit_all));
//                k += 1;
//            }
//
//        } else{
//            // TODO process contigs with 0 ORFs
//        }
//    }
//    if (verbosity) {
//        progressbar.finish();
//    }
//
//    // resort by longest_ORF idx
//    std::sort(ORF_coords_longest_tensorscore.begin(), ORF_coords_longest_tensorscore.end(), ORFsort_tensor);
//
//
//    // calculate scores for all ORFs by indexing into longest ORF scores
//    if (verbosity){
//        std::cout << "Calculating geneiness for all open reading frames..." << std::endl;
//    }
//    std::vector<std::pair<int, double>> ORF_score_sort;
//
//    int longorf_idx = 0;
//    int longorf_threeprime = ORF_coords_longest[0].second;
//    torch::Tensor * longorf_logit = &ORF_coords_longest_tensorscore[0].second;
//    int orf_threeprime;
//    // forward
//    for (auto & O : ORFs_forward){
//        int ORF_len = (O.first.second - O.first.first) / 3;
//        orf_threeprime = O.first.second;
//
//        // update index and logit scores for new corresponding longest ORF
//        if (orf_threeprime != longorf_threeprime){
//            longorf_idx ++;
//            longorf_threeprime = ORF_coords_longest[longorf_idx].second;
//            longorf_logit = &ORF_coords_longest_tensorscore[longorf_idx].second;
//        }
//
//        // average logits up to actual length in 3' to 5' direction
//        torch::Tensor logit_mean = torch::mean(longorf_logit->index({torch::indexing::Slice(torch::indexing::None, ORF_len)}));
//
//        // expit to get back to probability space
//        auto logit_score = (1 / (1 + torch::exp(-logit_mean))).item<double>();
//
//        ORF_score_sort.emplace_back(std::make_pair(O.second, logit_score));
//    }
//    // reverse
//    for (auto & O : ORFs_reverse){
//        int ORF_len = (O.first.first - O.first.second) / 3;
//        orf_threeprime = O.first.second;
//
//        // update index and logit scores for new corresponding longest ORF
//        if (orf_threeprime != longorf_threeprime){
//            longorf_idx ++;
//            longorf_threeprime = ORF_coords_longest[longorf_idx].second;
//            longorf_logit = &ORF_coords_longest_tensorscore[longorf_idx].second;
//        }
//
//        // average logits up to actual length in 3' to 5' direction
//        torch::Tensor logit_mean = torch::mean(longorf_logit->index({torch::indexing::Slice(torch::indexing::None, ORF_len)}));
//
//        // expit to get back to probability space
//        auto logit_score = (1 / (1 + torch::exp(-logit_mean))).item<double>();
//
//        ORF_score_sort.emplace_back(std::make_pair(O.second, logit_score));
//    }
//
//    // save scores
//    std::sort(ORF_score_sort.begin(), ORF_score_sort.end());
//    for (auto x: ORF_score_sort){
//        gene_model_scores.emplace_back(x.second);
//    }
//}

void GeneFinder::get_gene_model_scores() {
    // get lengths for later scoring
    for (auto & ORF_coord : ORF_coords) {
        int l = std::abs(ORF_coord.first - ORF_coord.second);
        ORF_lengths.emplace_back(l);
    }

    // separate strands
    std::vector<std::pair<std::pair<int, int>, int>> ORFs_forward;
    std::vector<std::pair<std::pair<int, int>, int>> ORFs_reverse;

    for (int i=0; i < ORF_coords.size(); ++i) {
        std::pair<std::pair<int, int>, int> coord_idx = std::make_pair(ORF_coords[i], i);
        if (coord_idx.first.first < coord_idx.first.second){
            ORFs_forward.emplace_back(coord_idx);
        } else{
            ORFs_reverse.emplace_back(coord_idx);
        }
    }

    // sort by stop codon, with ORFs sorted by length within shared stop codon
    std::sort(ORFs_forward.begin(), ORFs_forward.end(), threeprimefiveprimesort_ORFs);
    std::sort(ORFs_reverse.begin(), ORFs_reverse.end(), threeprimefiveprimesort_ORFs);

    // extract longest ORF per stop codon
    std::vector<std::pair<int, int>> ORF_coords_longest;

    // forward
    int scoredstop = -1;
    for (auto & O : ORFs_forward){
        auto coords = O.first;
        int stop = coords.second;
        if (scoredstop == stop){
            continue;
        }
        scoredstop = stop;

        ORF_coords_longest.emplace_back(coords);
    }
    // reverse
    scoredstop = -1;
    for (auto & O : ORFs_reverse){
        auto coords = O.first;
        int stop = coords.second;
        if (scoredstop == stop){
            continue;
        }
        scoredstop = stop;

        ORF_coords_longest.emplace_back(coords);
    }

    // sort by length, keep track of original index
    std::vector<std::tuple<int, std::pair<int, int>, int>> ORF_lengths_coords_sort;
    for (int i=0; i < ORF_coords_longest.size(); ++i){
        int l = std::abs(ORF_coords_longest[i].first-ORF_coords_longest[i].second);
        ORF_lengths_coords_sort.emplace_back(std::make_tuple(l, ORF_coords_longest[i], i));
    }
    std::sort(ORF_lengths_coords_sort.begin(), ORF_lengths_coords_sort.end());




    // extract amino acid sequence of longest ORFs for pytorch scoring
    std::pair<int, int> coords;
    std::vector<char> aa_seqveq_tmp;

    for (int i=0; i < ORF_lengths_coords_sort.size(); i++){
        // index into translated character vector based on coordinates
        coords = std::get<1>(ORF_lengths_coords_sort[i]);
        if (coords.first < coords.second) { // forward strand
            if (coords.second%3 == 0) {
                auto begin = seq_translated_f0.begin() + coords.first / 3 + 1;
                auto end = begin + (coords.second - coords.first) / 3 - 1;
                aa_seqveq_tmp = std::vector<char>(begin, end);
            } else if (coords.second%3 == 1) {
                auto begin = seq_translated_f1.begin() + coords.first / 3 + 1;
                auto end = begin + (coords.second - coords.first) / 3 - 1;
                aa_seqveq_tmp = std::vector<char>(begin, end);
            } else {
                auto begin = seq_translated_f2.begin() + coords.first / 3 + 1;
                auto end = begin + (coords.second - coords.first) / 3 - 1;
                aa_seqveq_tmp = std::vector<char>(begin, end);
            }
            std::reverse(aa_seqveq_tmp.begin(), aa_seqveq_tmp.end()); // reverse forward strand so sequence is 3' to 5'
        } else { // reverse strand
            if (coords.first%3 == 0) {
                auto begin = seq_translated_r0.begin() + coords.second / 3 + 1;
                auto end = begin + (coords.first - coords.second) / 3 - 1;
                aa_seqveq_tmp = std::vector<char>(begin, end);
            } else if (coords.first%3 == 1){
                auto begin = seq_translated_r1.begin() + coords.second / 3 + 1;
                auto end = begin + (coords.first - coords.second) / 3 - 1;
                aa_seqveq_tmp = std::vector<char>(begin, end);
            } else {
                auto begin = seq_translated_r2.begin() + coords.second / 3 + 1;
                auto end = begin + (coords.first - coords.second) / 3 - 1;
                aa_seqveq_tmp = std::vector<char>(begin, end);
            }
        }
        ORF_protein_seq_3to5_nostart_nostop.emplace_back(aa_seqveq_tmp);
    }

    std::string tmp_ORF_path = tmp_dir + "balrog.tmp" + std::to_string(mt_rand());
    std::ofstream stream;
    stream.open(tmp_ORF_path);
    for (const std::vector<char>& s : ORF_protein_seq_3to5_nostart_nostop){
        for (char c : s){
            stream << c;
        }
        stream << "\n";
    }
    stream.close();

    // run pytorch gene model
    // this is a mega janky way to do this, but still less janky than LibTorch itself
    std::vector<std::string> py_score_ORFs;
    py_score_ORFs.emplace_back("import torch");
    py_score_ORFs.emplace_back("import os");
    py_score_ORFs.emplace_back("import sys");
    py_score_ORFs.emplace_back("import time");
    py_score_ORFs.emplace_back("import numpy as np");
    py_score_ORFs.emplace_back("import torch.nn as nn");
    py_score_ORFs.emplace_back("import torch.nn.functional as F");
    py_score_ORFs.emplace_back("from torch.nn.utils import weight_norm");
    py_score_ORFs.emplace_back("from scipy.special import expit");
    py_score_ORFs.emplace_back("from scipy.special import logit");
    py_score_ORFs.emplace_back("");
    py_score_ORFs.emplace_back("aa_ORF_seq_path = \"" + tmp_ORF_path + "\"");
    py_score_ORFs.emplace_back("with open(aa_ORF_seq_path, \"rt\") as f:\n"
                               "    aa_ORF_seq = f.readlines()");
    py_score_ORFs.emplace_back(R"(aa_ORF_seq = [x.rstrip("\n") for x in aa_ORF_seq])");
    py_score_ORFs.emplace_back("");
    py_score_ORFs.emplace_back("gene_batch_size = " + std::to_string(gene_batch));
    py_score_ORFs.emplace_back("");
    py_score_ORFs.emplace_back("# load gene model");
    py_score_ORFs.emplace_back("repo = \"Markusjsommer/balrog_models\"");
    py_score_ORFs.emplace_back("torch.hub.set_dir(\"" + tmp_dir + "\")");

    py_score_ORFs.emplace_back("if torch.cuda.device_count() > 0:\n"
                               "    print(\"GPU detected...\")\n"
                               "    print(\"Loading gene model...\")\n"
                               "    model = torch.hub.load(repo, \"geneTCN\", force_reload=False).cuda()\n"
                               "    time.sleep(0.5)\n"
                               "else:\n"
                               "    print(\"No GPU detected, using CPU...\")\n"
                               "    print(\"Loading gene model...\")\n"
                               "    model = torch.hub.load(repo, \"geneTCN\", force_reload=False)\n"
                               "    time.sleep(0.5)\n");
    py_score_ORFs.emplace_back("");
    py_score_ORFs.emplace_back("def tokenize_aa_seq(aa_seq):\n"
                               "    \"\"\" Convert amino acid letters to integers.\"\"\"\n"
                               "    table = {\"L\": 1,\n"
                               "             \"V\": 2,\n"
                               "             \"I\": 3,\n"
                               "             \"M\": 4,\n"
                               "             \"C\": 5,\n"
                               "             \"A\": 6,\n"
                               "             \"G\": 7,\n"
                               "             \"S\": 8,\n"
                               "             \"T\": 9,\n"
                               "             \"P\": 10,\n"
                               "             \"F\": 11,\n"
                               "             \"Y\": 12,\n"
                               "             \"W\": 13,\n"
                               "             \"E\": 14,\n"
                               "             \"D\": 15,\n"
                               "             \"N\": 16,\n"
                               "             \"Q\": 17,\n"
                               "             \"K\": 18,\n"
                               "             \"R\": 19,\n"
                               "             \"H\": 20,\n"
                               "             \"*\": 0,\n"
                               "             \"X\": 0}\n"
                               "    tokenized = torch.tensor([table[aa] for aa in aa_seq])\n"
                               "    return tokenized");
    py_score_ORFs.emplace_back("");
    py_score_ORFs.emplace_back("def predict(X):\n"
                               "    model.eval()\n"
                               "    with torch.no_grad():\n"
                               "        if torch.cuda.device_count() > 0:\n"
                               "            X_enc = F.one_hot(X, 21).permute(0,2,1).float().cuda()\n"
                               "            probs = expit(model(X_enc).cpu())\n"
                               "            del X_enc\n"
                               "            torch.cuda.empty_cache()\n"
                               "        else:\n"
                               "            X_enc = F.one_hot(X, 21).permute(0,2,1).float()\n"
                               "            probs = expit(model(X_enc).cpu())\n"
                               "    return probs");
    py_score_ORFs.emplace_back("");
    py_score_ORFs.emplace_back("# encode amino acids as integers");
    py_score_ORFs.emplace_back("ORF_seq_enc = [tokenize_aa_seq(x) for x in aa_ORF_seq]");
    py_score_ORFs.emplace_back("");
    py_score_ORFs.emplace_back("# sort by length to minimize impact of batch padding");
    py_score_ORFs.emplace_back("ORF_lengths = np.asarray([len(x) for x in ORF_seq_enc])");
    py_score_ORFs.emplace_back("length_idx = np.argsort(ORF_lengths)");
    py_score_ORFs.emplace_back("ORF_seq_sorted = [ORF_seq_enc[i] for i in length_idx]");
    py_score_ORFs.emplace_back("");
    py_score_ORFs.emplace_back("# pad to allow creation of batch matrix\n"
                               "print(\"Scoring ORFs with gene model...\\n\")\n"
                               "prob_list = []\n"
                               "for i in range(0, len(ORF_seq_sorted), gene_batch_size):\n"
                               "    if (i%gene_batch_size==0):\n"
                               "        print (\"\\033[A                             \\033[A\")\n"
                               "        print(i, \" of \", len(ORF_seq_sorted))\n"
                               "    batch = ORF_seq_sorted[i:i+gene_batch_size]\n"
                               "    seq_lengths = torch.LongTensor(list(map(len, batch)))\n"
                               "    seq_tensor = torch.zeros((len(batch), seq_lengths.max())).long()\n"
                               "\n"
                               "    for idx, (seq, seqlen) in enumerate(zip(batch, seq_lengths)):\n"
                               "        seq_tensor[idx, :seqlen] = torch.LongTensor(seq)\n"
                               "\n"
                               "    pred_all = predict(seq_tensor)\n"
                               "\n"
                               "    pred = []\n"
                               "    for j, length in enumerate(seq_lengths):\n"
                               "        subseq = pred_all[j, 0, 0:int(length)]\n"
                               "        predprob = float(expit(torch.mean(logit(subseq))))\n"
//                               "        pred.append(predprob)\n"
                               "        subseq = logit(pred_all[j, 0, 0:int(length)])\n"
                               "        pred.append(subseq)\n" // need all logits
                               "    prob_list.extend(pred)\n"

                               "print (\"\\033[A                             \\033[A\")\n"
                               "print(len(ORF_seq_sorted), \" of \", len(ORF_seq_sorted))\n");
    py_score_ORFs.emplace_back("# unsort\n"
                               "unsort_idx = np.argsort(length_idx)\n"
//                               "ORF_prob = prob_arr[unsort_idx]");
                               "ORF_prob = []\n"
                               "for idx in unsort_idx:\n"
                               "    ORF_prob.append(prob_list[idx])");

    // save to temp file
    std::string tmp_gene_logit_path = tmp_dir + "balrog.tmp" + std::to_string(mt_rand());

    py_score_ORFs.emplace_back("with open(\"" + tmp_gene_logit_path + "\", \"wt\") as f :");
//    py_score_ORFs.emplace_back("    for tensor in prob_list:\n"
    py_score_ORFs.emplace_back("    for tensor in ORF_prob:\n"
                               "        for x in tensor:\n"
                               "            f.write(str(float(x)) + \" \")\n"
                               "        f.write(\"\\n\")");
    py_score_ORFs.emplace_back("");

    // write temp python file
    std::string tmp_score_ORFs_py_path = tmp_dir + "balrog.tmp" + std::to_string(mt_rand());
    stream.open(tmp_score_ORFs_py_path);
    for (const std::string& l : py_score_ORFs){
        stream << l << "\n";
    }
    stream.close();

    std::string command = "python " + tmp_score_ORFs_py_path;

    // run gene scoring script
    int status = std::system(command.c_str());
    if (status != 0) {
        std::cerr << "error running PyTorch, balrog requires python with PyTorch >= 1.5.0 to be callable from the environment where balrog is run\n";
        exit(status);
    }

    // read scores from temp file
    std::vector<std::vector<float>> ORF_longest_logit;
    std::ifstream ORF_logit_file (tmp_gene_logit_path);
    std::string line;
    if (ORF_logit_file.is_open()){
        while ( getline(ORF_logit_file,line) ){
            std::vector<float> v;
            std::istringstream iss(line);
            std::copy(std::istream_iterator<float>(iss),
                      std::istream_iterator<float>(),
                      std::back_inserter(v));
            ORF_longest_logit.emplace_back(v);
        }
        ORF_logit_file.close();
    }

    // clean tmp files
    std::remove(tmp_ORF_path.c_str());
    std::remove(tmp_gene_logit_path.c_str());
    std::remove(tmp_score_ORFs_py_path.c_str());

    // resort by longest_ORF idx
    std::vector<std::pair<int, std::vector<float>>> ORF_coords_longest_tensorscore;
    std::pair<int, std::vector<float>> ORF_idx_scorevec;
    int idx = 0;
    for (auto x : ORF_lengths_coords_sort){
        ORF_idx_scorevec = std::make_pair(std::get<2>(x), ORF_longest_logit[idx]);
        ORF_coords_longest_tensorscore.emplace_back(ORF_idx_scorevec);
        idx ++;
    }
    std::sort(ORF_coords_longest_tensorscore.begin(), ORF_coords_longest_tensorscore.end(), ORFsort_vector);

//
//    std::vector<int> foo_lens; // THIS CODE SHOWS A DIFFERENCE IN ORF LENGTH IN LIBTORCH RELIANT VERSION
//    int l;
//    for (auto x : ORF_coords_longest){
//        l = std::abs(x.first - x.second)/3;
//        foo_lens.emplace_back(l);
//    }
//
//
//    for (int i = 0; i < 20; i++){
//        std::cout << ORF_coords_longest_tensorscore[i].second.size() << std::endl;
//        std::cout << foo_lens[i] << std::endl << std::endl;
//    }


    // calculate scores for all ORFs by indexing into longest ORF scores
    if (verbosity){
        std::cout << "Calculating geneiness for all open reading frames..." << std::endl;
    }
    std::vector<std::pair<int, double>> ORF_score_sort;

    int longorf_idx = 0;
    int longorf_threeprime = ORF_coords_longest[0].second;
    std::vector<float> longorf_logit = ORF_coords_longest_tensorscore[0].second;
    int orf_threeprime;
    // forward
    for (auto & O : ORFs_forward){
        int ORF_len = (O.first.second - O.first.first) / 3;
        orf_threeprime = O.first.second;

        // update index and logit scores for new corresponding longest ORF
        if (orf_threeprime != longorf_threeprime){
            longorf_idx ++;
            longorf_threeprime = ORF_coords_longest[longorf_idx].second;
            longorf_logit = ORF_coords_longest_tensorscore[longorf_idx].second;
        }

        // average logits up to actual length in 3' to 5' direction
        float logit_sum = 0;
        for (int i=0; i < ORF_len; ++i){
            logit_sum += longorf_logit[i];
        }
        float logit_mean;
        logit_mean = logit_sum / (float)ORF_len;

        // expit to get back to probability space
        auto expit_score = (1 / (1 + std::exp(-logit_mean)));
        ORF_score_sort.emplace_back(std::make_pair(O.second, expit_score));
    }
    // reverse
    for (auto & O : ORFs_reverse){
        int ORF_len = (O.first.first - O.first.second) / 3;
        orf_threeprime = O.first.second;

        // update index and logit scores for new corresponding longest ORF
        if (orf_threeprime != longorf_threeprime){
            longorf_idx ++;
            longorf_threeprime = ORF_coords_longest[longorf_idx].second;
            longorf_logit = ORF_coords_longest_tensorscore[longorf_idx].second;
        }

        // average logits up to actual length in 3' to 5' direction
        float logit_sum = 0;
        for (int i=0; i < ORF_len; ++i){
            logit_sum += longorf_logit[i];
        }
        float logit_mean;
        logit_mean = logit_sum / (float)ORF_len;

        // expit to get back to probability space
        auto logit_score = (1 / (1 + std::exp(-logit_mean)));

        ORF_score_sort.emplace_back(std::make_pair(O.second, logit_score));
    }

    // save scores
    std::sort(ORF_score_sort.begin(), ORF_score_sort.end());
    for (auto x: ORF_score_sort){
        gene_model_scores.emplace_back(x.second);
    }
}

void GeneFinder::get_TIS_model_scores(){
    // extract nucleotide sequence flanking ORF start sites
    std::vector<std::string> ORF_TIS_seq;

    int seqlen = seq.size();
    for (int i = 0; i < ORF_coords.size(); i += 1) {
        std::pair<int, int> coords = ORF_coords[i];
        std::string TIS;

        // forward strand
        if (coords.first < coords.second) {
            if (coords.first >= 16 and seqlen - coords.first >= 19) {
                std::string window = seq.substr(coords.first - 16,
                                                35); // flanking 16 bp up and down stream, excluding start codon
                std::string start_codon = window.substr(16, 3);
                ORF_start_codon.emplace_back(start_codon); // save start codona
                window.erase(16, 3); // remove start codon from window
                for (int k = 31; k >= 0; --k) { // TIS is scored in 3' to 5' direction
                    char nuc = window.at(k);
                    TIS += nuc;
                }
            } else {
                // TIS fragment
                for (int k = 0; k < 32; ++k) {
                    TIS += 'A';
                }
            }
            ORF_TIS_seq.emplace_back(TIS);
        }
        else {
            if (coords.first >= 16 and seqlen - coords.first >= 19) {
                std::string window = seq.substr(coords.first - 16, 35); // flanking 16 bp up and down stream, excluding start codon
                complement(window, window);
                std::string start_codon = window.substr(16, 3);
                std::swap(start_codon[0], start_codon[2]);
                ORF_start_codon.emplace_back(start_codon); // save start codon
                window.erase(16, 3); // remove start codon from window
                for (int k = 0; k < 32; ++k){ // TIS is scored in 3' to 5' direction
                    char nuc = window.at(k);
                    TIS += nuc;
                }
            } else {
                // TIS fragment
                for (int k = 0; k < 32; ++k){
                    TIS += 'A';

                }
            }
            ORF_TIS_seq.emplace_back(TIS);
        }
    }

    // save TIS seq to temp file
    std::string tmp_TIS_path = tmp_dir + "balrog.tmp" + std::to_string(mt_rand());
    std::ofstream stream;
    stream.open(tmp_TIS_path);
    for (const std::string& s : ORF_TIS_seq){
        for (char c : s){
            stream << c;
        }
        stream << "\n";
    }
    stream.close();

    // score TIS
    // this is a mega janky way to do this, but still less janky than LibTorch itself
    std::vector<std::string> py_score_TIS;
    py_score_TIS.emplace_back("import torch");
    py_score_TIS.emplace_back("import os");
    py_score_TIS.emplace_back("import sys");
    py_score_TIS.emplace_back("import time");
    py_score_TIS.emplace_back("import numpy as np");
    py_score_TIS.emplace_back("import torch.nn as nn");
    py_score_TIS.emplace_back("import torch.nn.functional as F");
    py_score_TIS.emplace_back("from torch.nn.utils import weight_norm");
    py_score_TIS.emplace_back("from scipy.special import expit");
    py_score_TIS.emplace_back("from scipy.special import logit");
    py_score_TIS.emplace_back("");
    py_score_TIS.emplace_back("");
    py_score_TIS.emplace_back("TIS_seq_path = \"" + tmp_TIS_path + "\"");
    py_score_TIS.emplace_back("with open(TIS_seq_path, \"rt\") as f:\n"
                               "    TIS_seq = f.readlines()");
    py_score_TIS.emplace_back(R"(TIS_seq = [x.rstrip("\n") for x in TIS_seq])");
    py_score_TIS.emplace_back("");
    py_score_TIS.emplace_back("TIS_batch_size = " + std::to_string(TIS_batch));
    py_score_TIS.emplace_back("");
    py_score_TIS.emplace_back("# load gene model");
    py_score_TIS.emplace_back("repo = \"Markusjsommer/balrog_models\"");
    py_score_TIS.emplace_back("torch.hub.set_dir(\"" + tmp_dir + "\")");

    py_score_TIS.emplace_back("if torch.cuda.device_count() > 0:\n"
                               "    print(\"GPU detected...\")\n"
                               "    print(\"Loading TIS model...\")\n"
                               "    model_tis = torch.hub.load(repo, \"tisTCN\", force_reload=False).cuda()\n"
                               "    time.sleep(0.5)\n"
                               "else:\n"
                               "    print(\"No GPU detected, using CPU...\")\n"
                               "    print(\"Loading TIS model...\")\n"
                               "    model_tis = torch.hub.load(repo, \"tisTCN\", force_reload=False)\n"
                               "    time.sleep(0.5)\n");
    py_score_TIS.emplace_back("");
    py_score_TIS.emplace_back("def predict_tis(X):\n"
                              "    model_tis.eval()\n"
                              "    with torch.no_grad():\n"
                              "        if torch.cuda.device_count() > 0:\n"
                              "            X_enc = F.one_hot(X, 4).permute(0,2,1).float().cuda()\n"
                              "        else:\n"
                              "            X_enc = F.one_hot(X, 4).permute(0,2,1).float()\n"
                              "        probs = expit(model_tis(X_enc).cpu())\n"
                              "    return probs");
    py_score_TIS.emplace_back("");
    py_score_TIS.emplace_back("nuc_encode = {\"A\":0,\n"
                              "              \"T\":1,\n"
                              "              \"G\":2,\n"
                              "              \"C\":3,\n"
                              "              \"N\":0,\n"
                              "              \"M\":0,\n"
                              "              \"R\":0,\n"
                              "              \"Y\":0,\n"
                              "              \"W\":0,\n"
                              "              \"K\":0}\n"
                              "              ");
    py_score_TIS.emplace_back("");
    py_score_TIS.emplace_back("TIS_seq_enc = [torch.tensor([nuc_encode[c] for c in TIS], dtype=int) for TIS in TIS_seq]");
    py_score_TIS.emplace_back("");
    py_score_TIS.emplace_back("# batch score TIS\n"
                              "TIS_prob_list = []\n"
                              "print()\n"
                              "for i in range(0, len(TIS_seq), TIS_batch_size):\n"
                              "    if (i%TIS_batch_size==0):\n"
                              "        print (\"\\033[A                             \\033[A\")\n"
                              "        print(i, \" of \", len(TIS_seq))\n"
                              "    batch = TIS_seq_enc[i:i+TIS_batch_size]\n"
                              "    TIS_stacked = torch.stack(batch)\n"
                              "    pred = predict_tis(TIS_stacked)\n"
                              "    TIS_prob_list.extend(pred)\n"
                              "y_pred_TIS = np.asarray(TIS_prob_list, dtype=float)\n"
                              "print (\"\\033[A                             \\033[A\")\n"
                              "print(len(TIS_seq), \" of \", len(TIS_seq))\n");

    // save to temp file
    std::string tmp_TIS_prob_path = tmp_dir + "balrog.tmp" + std::to_string(mt_rand());

    py_score_TIS.emplace_back("with open(\"" + tmp_TIS_prob_path + "\", \"wt\") as f :");
    py_score_TIS.emplace_back("    for prob in y_pred_TIS:\n"
                              "        f.write(str(float(prob)) + \"\\n\")");

    // write temp python file
    std::string tmp_score_TIS_py_path = tmp_dir + "balrog.tmp" + std::to_string(mt_rand());
    stream.open(tmp_score_TIS_py_path);
    for (const std::string& l : py_score_TIS){
        stream << l << "\n";
    }
    stream.close();

    std::string command = "python " + tmp_score_TIS_py_path;

    // run gene scoring script
    int status = std::system(command.c_str());
    if (status != 0) {
        std::cerr << "error running PyTorch, balrog requires python with PyTorch >= 1.5.0 to be callable from the environment where balrog is run\n";
        exit(status);
    }

    // read TIS scoresfrom temp file
    std::ifstream TIS_prob_file (tmp_TIS_prob_path);
    std::string line;
    if (TIS_prob_file.is_open()){
        while (getline(TIS_prob_file,line) ){
            double prob;
            prob = std::atof(line.c_str());
            TIS_model_scores.emplace_back(prob); // save scores
        }
        TIS_prob_file.close();
    }

    // clean temp files
    std::remove(tmp_TIS_path.c_str());
    std::remove(tmp_score_TIS_py_path.c_str());
    std::remove(tmp_TIS_prob_path.c_str());


}

void GeneFinder::translate_seq() {
    // TODO might be faster to use C++17 string indexing rather than a bunch of copies
    // NOTE: this was lazily converted from libtorch tensor pair encoding
    std::string codon;
    char res;
    char aa;

    // reserve memory for translated and tensorized sequence // TODO does this actually save time?
    seq_translated_f0.reserve(seq.length()/3);
    seq_translated_f1.reserve(seq.length()/3);
    seq_translated_f2.reserve(seq.length()/3);
    seq_translated_r0.reserve(seq.length()/3);
    seq_translated_r1.reserve(seq.length()/3);
    seq_translated_r2.reserve(seq.length()/3);

    for (int idx=0; idx < seq.length() - 3; idx += 3){ // TODO ignoring a few codons at the end
        // frame f0
        codon = seq.substr(idx, 3);
        res = trantab_standard[codon];
        aa = res;
        if (aa){
            seq_translated_f0 += aa;
        } else{
            seq_translated_f0 += 'X'; // catch non-standard characters
        }
        // frame r0
        std::swap(codon[0], codon[2]); // reverse codon order by swapping first and last
        complement(codon, codon);
        res = trantab_standard[codon];
        aa = res;
        if (aa){
            seq_translated_r0 += aa;
        } else{
            seq_translated_r0 += 'X'; // catch non-standard characters
        }
        // frame f1
        codon = seq.substr(idx + 1, 3);
        res = trantab_standard[codon];
        aa = res;
        if (aa){
            seq_translated_f1 += aa;
        } else{
            seq_translated_f1 += 'X'; // catch non-standard characters
        }
        // frame r1
        std::swap(codon[0], codon[2]); // reverse codon order by swapping first and last
        complement(codon, codon);
        res = trantab_standard[codon];
        aa = res;
        if (aa){
            seq_translated_r1 += aa;
        } else{
            seq_translated_r1 += 'X'; // catch non-standard characters
        }
        // frame f2
        codon = seq.substr(idx + 2, 3);
        res = trantab_standard[codon];
        aa = res;
        if (aa){
            seq_translated_f2 += aa;
        } else{
            seq_translated_f2 += 'X'; // catch non-standard characters
        }
        // frame r2
        std::swap(codon[0], codon[2]); // reverse codon order by swapping first and last
        complement(codon, codon);
        res = trantab_standard[codon];
        aa = res;
        if (aa){
            seq_translated_r2 += aa;
        } else{
            seq_translated_r2 += 'X'; // catch non-standard characters
        }
    }
}

void GeneFinder::find_ORFs() {
    // define possible start and stop codons
    std::vector<std::string> starts;
    starts.emplace_back("ATG");
    starts.emplace_back("GTG");
    starts.emplace_back("TTG");

    std::vector<std::string> starts_r;
    starts_r.emplace_back("CAT");
    starts_r.emplace_back("CAC");
    starts_r.emplace_back("CAA");

    std::vector<std::string> stops;
    stops.emplace_back("TAA");
    stops.emplace_back("TAG");
    if (trantab == 11){
        stops.emplace_back("TGA");
    }

    std::vector<std::string> stops_r;
    stops_r.emplace_back("TTA");
    stops_r.emplace_back("CTA");
    if (trantab == 11){
        stops_r.emplace_back("TCA");
    }

    // scan for ORFs in each reading frame, allowing starts and stops off contig ends
    // forward 0
    std::vector<int> temp_starts;
    temp_starts.emplace_back(0);
    for (int i=0; i < seq.size()-3; i+=3){
        std::string codon = seq.substr(i, 3);
        if (std::find(stops.begin(), stops.end(), codon) != stops.end()){ // TODO check if string_view is better
            for (int s : temp_starts){
                uint64_t ORF_length = i - s + 1;
                if (ORF_length >= minimum_length){
                    ORF_coords.emplace_back(std::make_pair(s, i)); // NOTE stop i is index of first nucleotide in stop codon, s is index of first nucelotide in start codon
                }
            }
            temp_starts.clear();
        }
        else if (std::find(starts.begin(), starts.end(), codon) != starts.end()){
            temp_starts.emplace_back(i);
        }
    }
    for (int s : temp_starts){
        uint64_t i = seq.size() - 1; // NOTE using index of last nucelotide for all ORFs extending past edge regardless of frame
        uint64_t ORF_length = i - s + 1;
        if (ORF_length >= minimum_length){
            ORF_coords.emplace_back(std::make_pair(s, i));
        }
    }
    // forward 1
    temp_starts.clear();
    temp_starts.emplace_back(1);
    for (int i=1; i < seq.size()-3; i+=3){ // NOTE using 0 for all ORFs extending past edge, regardless of frame
        std::string codon = seq.substr(i, 3);
        if (std::find(stops.begin(), stops.end(), codon) != stops.end()){
            for (int s : temp_starts){
                uint64_t ORF_length = i - s + 1;
                if (ORF_length >= minimum_length){
                    ORF_coords.emplace_back(std::make_pair(s, i));
                }
            }
            temp_starts.clear();
        }
        else if (std::find(starts.begin(), starts.end(), codon) != starts.end()){
            temp_starts.emplace_back(i);
        }
    }
    for (int s : temp_starts){
        uint64_t i = seq.size() - 1;
        uint64_t ORF_length = i - s + 1;
        if (ORF_length >= minimum_length){
            ORF_coords.emplace_back(std::make_pair(s, i));
        }
    }
    // forward 2
    temp_starts.clear();
    temp_starts.emplace_back(0);
    for (int i=2; i < seq.size()-3; i+=3){
        std::string codon = seq.substr(i, 3);
        if (std::find(stops.begin(), stops.end(), codon) != stops.end()){
            for (int s : temp_starts){
                uint64_t ORF_length = i - s + 1;
                if (ORF_length >= minimum_length){
                    ORF_coords.emplace_back(std::make_pair(s, i));
                }
            }
            temp_starts.clear();
        }
        else if (std::find(starts.begin(), starts.end(), codon) != starts.end()){
            temp_starts.emplace_back(i);
        }
    }
    for (int s : temp_starts){
        uint64_t i = seq.size() - 1;
        uint64_t ORF_length = i - s + 1;
        if (ORF_length >= minimum_length){
            ORF_coords.emplace_back(std::make_pair(s, i));
        }
    }

    // reverse ending at 0
    temp_starts.clear();
    int stop_idx = 0; // use temporary stop idx to keep track of stop for future starts
    for (int i=0; i < seq.size()-3; i+=3){
        std::string codon_r = seq.substr(i, 3);
        if (std::find(stops_r.begin(), stops_r.end(), codon_r) != stops_r.end()) {
            for (int s : temp_starts) {
                uint64_t ORF_length = s - stop_idx + 1;
                if (ORF_length >= minimum_length) {
                    if (stop_idx == 0){ // check for end of strand to be consistent with coords
                        ORF_coords.emplace_back(std::make_pair(s, 0));
                    } else {
                        ORF_coords.emplace_back(std::make_pair(s, stop_idx)); // NOTE add 2 to account for strand switch
                    }
                }
            }
            temp_starts.clear();
            stop_idx = i; // NOTE offsetting stop idx for strand direction here
        }
        else if (std::find(starts_r.begin(), starts_r.end(), codon_r) != starts_r.end()){
            temp_starts.emplace_back(i);
        }
    }
    temp_starts.emplace_back(seq.size() - 1); // allow start past end
    for (int s : temp_starts){
        uint64_t ORF_length = s - stop_idx + 1;
        if (ORF_length >= minimum_length){
            ORF_coords.emplace_back(std::make_pair(s, stop_idx));
        }
    }
    // reverse ending at 1
    temp_starts.clear();
    stop_idx = 0; // use temporary stop idx to keep track of stop for future starts
    for (int i=1; i < seq.size()-3; i+=3){
        std::string codon_r = seq.substr(i, 3);
        if (std::find(stops_r.begin(), stops_r.end(), codon_r) != stops_r.end()) {
            for (int s : temp_starts) {
                uint64_t ORF_length = s - stop_idx + 1;
                if (ORF_length >= minimum_length) {
                    if (stop_idx == 0){ // check for end of strand to be consistent with coords
                        ORF_coords.emplace_back(std::make_pair(s, 0));
                    } else {
                        ORF_coords.emplace_back(std::make_pair(s, stop_idx)); // NOTE add 2 to account for strand switch
                    }
                }
            }
            temp_starts.clear();
            stop_idx = i; // NOTE offsetting stop idx for strand direction here
        }
        else if (std::find(starts_r.begin(), starts_r.end(), codon_r) != starts_r.end()){
            temp_starts.emplace_back(i);
        }
    }
    temp_starts.emplace_back(seq.size() - 1); // allow start past end
    for (int s : temp_starts){
        uint64_t ORF_length = s - stop_idx + 1;
        if (ORF_length >= minimum_length){
            ORF_coords.emplace_back(std::make_pair(s, stop_idx));
        }
    }
    // reverse ending at 2
    temp_starts.clear();
    stop_idx = 0; // use temporary stop idx to keep track of stop for future starts
    for (int i=2; i < seq.size()-3; i+=3){
        std::string codon_r = seq.substr(i, 3);
        if (std::find(stops_r.begin(), stops_r.end(), codon_r) != stops_r.end()) {
            for (int s : temp_starts) {
                uint64_t ORF_length = s - stop_idx + 1;
                if (ORF_length >= minimum_length) {
                    if (stop_idx == 0){ // check for end of strand to be consistent with coords
                        ORF_coords.emplace_back(std::make_pair(s, 0));
                    } else {
                        ORF_coords.emplace_back(std::make_pair(s, stop_idx)); // NOTE add 2 to account for strand switch
                    }
                }
            }
            temp_starts.clear();
            stop_idx = i; // NOTE offsetting stop idx for strand direction here
        }
        else if (std::find(starts_r.begin(), starts_r.end(), codon_r) != starts_r.end()){
            temp_starts.emplace_back(i);
        }
    }
    temp_starts.emplace_back(seq.size() - 1); // allow start past end
    for (int s : temp_starts){
        uint64_t ORF_length = s - stop_idx + 1;
        if (ORF_length >= minimum_length){
            ORF_coords.emplace_back(std::make_pair(s, stop_idx));
        }
    }
}

void GeneFinder::get_start_codon_scores() {
    for (std::string s : ORF_start_codon){
        if (s == "ATG"){
            start_codon_scores.emplace_back(weight_ATG);
        } else if (s == "TTG"){
            start_codon_scores.emplace_back(weight_TTG);
        } else if (s == "GTG"){
            start_codon_scores.emplace_back(weight_GTG);
        } else {
            start_codon_scores.emplace_back(weight_GTG);
        }
    }
}

void GeneFinder::combine_scores() {
    for (int i=0; i < ORF_coords.size(); ++i){
        double score = 0;
        score += weight_gene_prob * gene_model_scores[i];
        score += weight_TIS_prob * TIS_model_scores[i];
        score += start_codon_scores[i];
        score -= probthresh;
        score *= ORF_lengths[i];
        ORF_scores.emplace_back(score);
    }
}


void GeneFinder::complement(std::string &s, std::string &comp) {
    std::string c;
    for (char ch : s){
        if (ch == 'A'){
            c += 'T';
        } else if (ch == 'T'){
            c += 'A';
        } else if (ch == 'G'){
            c += 'C';
        } else if (ch == 'C'){
            c += 'G';
        } else{
            c += ch; // just preserve non-ATGC nucleotides, they don't really matter
        }
    }
    comp = c;
}

void GeneFinder::maximize_coherent_score() {
    // topological sort by minimum coordinate
    std::vector<std::pair<std::pair<int, int>, int>> ORF_coords_topsort;
    ORF_coords_topsort.reserve(ORF_coords.size());
    for (int i=0; i < ORF_coords.size(); ++i) {
        int mincoord = std::min(ORF_coords[i].first, ORF_coords[i].second);
        int maxcoord = std::max(ORF_coords[i].first, ORF_coords[i].second);
        std::pair<int, int> coords_ordered = std::make_pair(mincoord, maxcoord);
        ORF_coords_topsort.emplace_back(std::make_pair(coords_ordered, i));
    }
    std::sort(ORF_coords_topsort.begin(), ORF_coords_topsort.end(), toposort_ORFs);

    // forward sweep, keeping track of max score at each node and parent of that max score
    std::vector<double> max_path_score(ORF_coords_topsort.size(), 0.0);
    std::vector<int> max_path_parent(ORF_coords_topsort.size(), 0);
    for (int i=0; i < ORF_coords_topsort.size(); ++i){
        int threeprime_parent = ORF_coords[ORF_coords_topsort[i].second].second;
        int fiveprime_parent = ORF_coords[ORF_coords_topsort[i].second].first;
        double parent_max_score = max_path_score[i];

        int n_forward_connections = 0;
        int child = i;
        while (n_forward_connections < max_forward_connections and child+1 < ORF_coords_topsort.size()){
            child ++;
            int threeprime_child = ORF_coords[ORF_coords_topsort[child].second].second;
            int fiveprime_child = ORF_coords[ORF_coords_topsort[child].second].first;

            // ensure 3' sites are different if on same strand
            if ((threeprime_parent < fiveprime_parent) == (threeprime_child < fiveprime_child)){
                if (threeprime_child == threeprime_parent){
                    continue;
                }
            }

            // check total overlap
            int overlap = ORF_coords_topsort[i].first.second - ORF_coords_topsort[child].first.first + 1;
            if (overlap > max_overlap){
                continue;
            }

            // penalize by overlap category // TODO double check the directionality here
            double overlap_penalty = 0;
            if (overlap > 0){
                // unidirectional
                if ((threeprime_child < fiveprime_child) and (threeprime_parent < fiveprime_parent)){
                    overlap_penalty += (double)overlap * unidirectional_penalty_per_base;
                }
                // convergent
                else if (((fiveprime_child < threeprime_parent and (threeprime_parent <= threeprime_child))) or ((fiveprime_parent < threeprime_child) and (threeprime_child <= threeprime_parent))){
                    overlap_penalty += (double)overlap * convergent_penalty_per_base;
                }
                // divergent
                else {
                    overlap_penalty += (double)overlap * divergent_penalty_per_base;
                }
            }

            // update maximum attainable score and parent of that score
            double edge_score = ORF_scores[ORF_coords_topsort[child].second] - overlap_penalty;
            double child_score = parent_max_score + edge_score;
            if (child_score > max_path_score[child]){
                max_path_score[child] = child_score;
                max_path_parent[child] = i;
                n_forward_connections ++;
            }
        }
    }

    // backtrack from max score to 0 to find highest scoring path
    std::vector<int> gene_ORF_idx_topsort;
    int graph_idx = std::distance(max_path_score.begin(), std::max_element(max_path_score.begin(), max_path_score.end()));
    gene_ORF_idx_topsort.emplace_back(graph_idx);
    double s = max_path_score[graph_idx];
    while (s != 0){
//    while (s > 0){
        graph_idx = max_path_parent[graph_idx];
        gene_ORF_idx_topsort.emplace_back(graph_idx);
        s = max_path_score[graph_idx];
    }

    // save predicted genes
    (*_gene_coord).reserve(gene_ORF_idx_topsort.size());
    (*_gene_strand).reserve(gene_ORF_idx_topsort.size());
    (*_gene_nucseq).reserve(gene_ORF_idx_topsort.size());
    (*_gene_protseq).reserve(gene_ORF_idx_topsort.size());
    (*_gene_score).reserve(gene_ORF_idx_topsort.size());

    for (int i = (int)gene_ORF_idx_topsort.size()-1; i >= 0; --i){ // flip order of backtracked graph
        // coords
//        std::pair<int, int> coords = ORF_coords_topsort[gene_ORF_idx_topsort[i]].first;
        std::pair<int, int> coords = ORF_coords[ORF_coords_topsort[gene_ORF_idx_topsort[i]].second];
        (*_gene_coord).emplace_back(coords);

        // strand
        (*_gene_strand).emplace_back(coords.first < coords.second);

        // nucleotide sequence 5' to 3'
        std::string nucseq;
        if ((*_gene_strand).back()){
            nucseq = seq.substr(coords.first,coords.second - coords.first);
            (*_gene_nucseq).emplace_back(nucseq);
        } else{
            int stopcoord = coords.second + 3;
            nucseq = seq.substr(stopcoord,coords.first - coords.second);

            std::reverse(nucseq.begin(), nucseq.end());
            std::string nucseq_rt;
            complement(nucseq, nucseq_rt);
            nucseq = nucseq_rt;
        }
        (*_gene_nucseq).emplace_back(nucseq);


        // protein sequence 5' to 3'
        std::string protseq;
        for (int j=0; j < nucseq.size(); j +=3){
            protseq += trantab_standard[nucseq.substr(j, 3)];
        }
        (*_gene_protseq).emplace_back(protseq);

        // get scores
        double score = ORF_scores[ORF_coords_topsort[gene_ORF_idx_topsort[i]].second];
        (*_gene_score).emplace_back(score);
    }
}

int GeneFinder::mt_rand() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, INT32_MAX);
    return distrib(gen);
}
