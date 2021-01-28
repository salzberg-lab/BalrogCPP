//
// Created by Markus on 9/13/2020.
//

#include <iostream>
#include <algorithm>
#include "GeneFinder.h"
#include <torch/script.h>
#include "tqdm.h"
#include <omp.h>

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
        std::cout << "Tensorizing open reading frames..." << std::endl;
    }
    translate_and_tensorize_seq();
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

bool ORFsort_tensor(const std::pair<int, torch::Tensor>& a,
                    const std::pair<int, torch::Tensor>& b){
    return (a.first < b.first);
}


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

    // score in length-sorted blocks
    std::vector<std::pair<int, torch::Tensor>> ORF_coords_longest_tensorscore;
    tqdm progressbar;
    progressbar.set_theme_braille_spin();
    for (int i=0; i < ORF_lengths_coords_sort.size(); i += gene_batch){
        std::vector<torch::Tensor> ORF_block;
        int block_len;
        if (i + gene_batch <= ORF_lengths_coords_sort.size()){
            block_len = std::get<0>(ORF_lengths_coords_sort[i+gene_batch]) / 3;
        } else{
            // last block may be of different size
            block_len = std::get<0>(ORF_lengths_coords_sort.back()) / 3;
        }
        // create length-sorted tensor blocks
        unsigned long block_end = std::min((unsigned long)(i + gene_batch), ORF_lengths_coords_sort.size());

        for (int j=i; j < block_end; ++j){
            if (verbosity) {
                progressbar.progress(j, (int) ORF_lengths_coords_sort.size());
            }
            // extract amino acid sequence in correct reading frame and add it to tensor block
            std::pair<int, int> coords = std::get<1>(ORF_lengths_coords_sort[j]);
            // forward strand
            if (coords.first < coords.second) {
                if (coords.second%3 == 0) {
                    auto begin = seq_tensorized_f0.begin() + coords.first / 3 + 1;
                    auto end = begin + (coords.second - coords.first) / 3 - 1;
                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
                    std::reverse(vec.begin(), vec.end()); // reverse so sequence is 3' to 5'
                    unsigned long pad = block_len - vec.size(); // pad so all tensors are the same size in each block
                    for (int k=0; k < pad; ++k){
                        vec.emplace_back(torch::zeros({1, 21,1}));
                    }
                    ORF_block.emplace_back(torch::cat(vec, -1));
                } else if (coords.second%3 == 1) {
                    auto begin = seq_tensorized_f1.begin() + coords.first / 3 + 1;
                    auto end = begin + (coords.second - coords.first) / 3 - 1;
                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
                    std::reverse(vec.begin(), vec.end());
                    unsigned long pad = block_len - vec.size();
                    for (int k=0; k < pad; ++k){
                        vec.emplace_back(torch::zeros({1, 21,1}));
                    }
                    ORF_block.emplace_back(torch::cat(vec, -1));
                } else {
                    auto begin = seq_tensorized_f2.begin() + coords.first / 3 + 1;
                    auto end = begin + (coords.second - coords.first) / 3 - 1;
                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
                    std::reverse(vec.begin(), vec.end());
                    unsigned long pad = block_len - vec.size();
                    for (int k=0; k < pad; ++k){
                        vec.emplace_back(torch::zeros({1, 21,1}));
                    }
                    ORF_block.emplace_back(torch::cat(vec, -1));
                }
            }
            // reverse strand
            else {
                if (coords.first%3 == 0) {
                    auto begin = seq_tensorized_r0.begin() + coords.second / 3 + 1;
                    auto end = begin + (coords.first - coords.second) / 3 - 1;
                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
                    unsigned long pad = block_len - vec.size();
                    for (int k=0; k < pad; ++k){
                        vec.emplace_back(torch::zeros({1, 21,1}));
                    }
                    ORF_block.emplace_back(torch::cat(vec, -1));
                } else if (coords.first%3 == 1){
                    auto begin = seq_tensorized_r1.begin() + coords.second / 3 + 1;
                    auto end = begin + (coords.first - coords.second) / 3 - 1;
                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
                    unsigned long pad = block_len - vec.size();
                    for (int k=0; k < pad; ++k){
                        vec.emplace_back(torch::zeros({1, 21,1}));
                    }
                    ORF_block.emplace_back(torch::cat(vec, -1));
                } else {
                    auto begin = seq_tensorized_r2.begin() + coords.second / 3 + 1;
                    auto end = begin + (coords.first - coords.second) / 3 - 1;
                    std::vector<torch::Tensor> vec = std::vector<torch::Tensor>(begin, end);
                    unsigned long pad = block_len - vec.size();
                    for (int k=0; k < pad; ++k){
                        vec.emplace_back(torch::zeros({1, 21,1}));
                    }
                    ORF_block.emplace_back(torch::cat(vec, -1));
                }
            }
        }
        // score each batch block
        if (not ORF_block.empty()) {
            // get values for last hidden unit layer {batch_size, hidden_units, length}
            torch::Tensor ORF_batch = torch::cat({ORF_block}, 0);
            torch::NoGradGuard nograd;
            _gene_model.eval();
            auto score = _gene_model.forward({ORF_batch});
            torch::Tensor score_tensor = score.toTensor();

            // apply linear function across all causal hidden units (equivalent to getting score at each character in the sequence)
            std::vector<torch::Tensor> linear_layer_score_vec;
            for (int j=0; j < score_tensor.sizes()[2]; ++j){
                torch::Tensor linear_layer_score = torch::matmul(score_tensor.index({torch::indexing::Slice(), torch::indexing::Slice(), j}), gene_model_linear_weights) + gene_model_linear_bias;
                linear_layer_score_vec.emplace_back(linear_layer_score);
            }

            // combine post-linear-layer scores into {batch_size, length} matrix
            torch::Tensor linear_layer_matrix = torch::stack({linear_layer_score_vec}, 1);

            int k = 0;
            for (int j=i; j < block_end; ++j) {
                torch::Tensor logit_all = linear_layer_matrix.index({k, torch::indexing::Slice(torch::indexing::None, torch::indexing::None)});
                int original_idx = std::get<2>(ORF_lengths_coords_sort[j]);
                ORF_coords_longest_tensorscore.emplace_back(std::make_pair(original_idx, logit_all));
                k += 1;
            }

        } else{
            // TODO process contigs with 0 ORFs
        }
    }
    if (verbosity) {
        progressbar.finish();
    }

    // resort by longest_ORF idx
    std::sort(ORF_coords_longest_tensorscore.begin(), ORF_coords_longest_tensorscore.end(), ORFsort_tensor);


    // calculate scores for all ORFs by indexing into longest ORF scores
    if (verbosity){
        std::cout << "Calculating geneiness for all open reading frames..." << std::endl;
    }
    std::vector<std::pair<int, double>> ORF_score_sort;

    int longorf_idx = 0;
    int longorf_threeprime = ORF_coords_longest[0].second;
    torch::Tensor * longorf_logit = &ORF_coords_longest_tensorscore[0].second;
    int orf_threeprime;
    // forward
    for (auto & O : ORFs_forward){
        int ORF_len = (O.first.second - O.first.first) / 3;
        orf_threeprime = O.first.second;

        // update index and logit scores for new corresponding longest ORF
        if (orf_threeprime != longorf_threeprime){
            longorf_idx ++;
            longorf_threeprime = ORF_coords_longest[longorf_idx].second;
            longorf_logit = &ORF_coords_longest_tensorscore[longorf_idx].second;
        }

        // average logits up to actual length in 3' to 5' direction
        torch::Tensor logit_mean = torch::mean(longorf_logit->index({torch::indexing::Slice(torch::indexing::None, ORF_len)}));

        // expit to get back to probability space
        auto logit_score = (1 / (1 + torch::exp(-logit_mean))).item<double>();

        ORF_score_sort.emplace_back(std::make_pair(O.second, logit_score));
    }
    // reverse
    for (auto & O : ORFs_reverse){
        int ORF_len = (O.first.first - O.first.second) / 3;
        orf_threeprime = O.first.second;

        // update index and logit scores for new corresponding longest ORF
        if (orf_threeprime != longorf_threeprime){
            longorf_idx ++;
            longorf_threeprime = ORF_coords_longest[longorf_idx].second;
            longorf_logit = &ORF_coords_longest_tensorscore[longorf_idx].second;
        }

        // average logits up to actual length in 3' to 5' direction
        torch::Tensor logit_mean = torch::mean(longorf_logit->index({torch::indexing::Slice(torch::indexing::None, ORF_len)}));

        // expit to get back to probability space
        auto logit_score = (1 / (1 + torch::exp(-logit_mean))).item<double>();

        ORF_score_sort.emplace_back(std::make_pair(O.second, logit_score));
    }

    // save scores
    std::sort(ORF_score_sort.begin(), ORF_score_sort.end());
    for (auto x: ORF_score_sort){
        gene_model_scores.emplace_back(x.second);
    }
}

// SCORES ALL ORFS RATHER THAN LONGEST ORFS, ~5X SLOWER
//void GeneFinder::get_gene_model_scores() {
//    // sort by length, keep track of original index
//    std::vector<std::tuple<int, std::pair<int, int>, int>> ORF_lengths_coords_sort;
//    for (int i=0; i < ORF_coords.size(); ++i){
//        int l = std::abs(ORF_coords[i].first-ORF_coords[i].second);
//        ORF_lengths.emplace_back(l);
//        ORF_lengths_coords_sort.emplace_back(std::make_tuple(l, ORF_coords[i], i));
//    }
//    std::sort(ORF_lengths_coords_sort.begin(), ORF_lengths_coords_sort.end());
//
//    // score in length-sorted blocks
//    std::vector<std::pair<int, double>> ORF_score_sort;
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
//            // logit-average up to the true index of the last character (end of matrix may be padding; exclude padding from score)
//            // NOTE: network already outputs logits, no need to logit again
//            int k = 0;
//            for (int j=i; j < block_end; ++j) {
//                int true_len = std::get<0>(ORF_lengths_coords_sort[j]) / 3;
//                torch::Tensor logit_mean = torch::mean(linear_layer_matrix.index({k, torch::indexing::Slice(torch::indexing::None, true_len)}));
//                auto logit_score = (1 / (1 + torch::exp(-logit_mean))).item<double>();
//                int original_idx = std::get<2>(ORF_lengths_coords_sort[j]);
//                ORF_score_sort.emplace_back(std::make_pair(original_idx, logit_score));
//                k += 1;
//            }
//        } else{
//            // TODO process contigs with 0 ORFs
//        }
//    }
//    if (verbosity) {
//        progressbar.finish();
//    }
//    // resort by original idx
//    std::sort(ORF_score_sort.begin(), ORF_score_sort.end());
//
//    // save scores
//    for (auto x: ORF_score_sort){
//        gene_model_scores.emplace_back(x.second);
//    }
//
//}

void GeneFinder::get_TIS_model_scores() {
    int seqlen = seq.size();

    tqdm progressbar;
    progressbar.set_theme_braille_spin();
    for (int i = 0; i < ORF_coords.size(); i += TIS_batch) {
        unsigned long block_end = std::min((unsigned long)(i + TIS_batch), ORF_coords.size());
        std::vector<torch::Tensor> ORF_block;

        for (int j = i; j < block_end; ++j){
            if (verbosity){
                progressbar.progress(j, (int)ORF_coords.size());
            }
            std::pair<int, int> coords = ORF_coords[j];

            // forward strand
            if (coords.first < coords.second) {
                std::vector<torch::Tensor> TIS_tensor_vec;
                if (coords.first >= 16 and seqlen - coords.first >= 19) {
                    std::string window = seq.substr(coords.first - 16, 35); // flanking 16 bp up and down stream, excluding start codon
                    std::string start_codon = window.substr(16, 3);
                    ORF_start_codon.emplace_back(start_codon); // save start codon
                    window.erase(16, 3); // remove start codon from window
                    for (int k = 31; k >= 0; --k){ // TIS is scored in 3' to 5' direction
                        char nuc = window.at(k);
                        torch::Tensor enc;
                        if (nuctab.count(nuc)){
                            enc = nuctab[nuc];
                        } else{ // non-ATGC nucleotides
                            enc = nuctab['A'];
                        }
                        TIS_tensor_vec.emplace_back(enc);
                    }
                } else {
                    // TIS fragment
                    for (int k = 0; k < 32; ++k){
                        TIS_tensor_vec.emplace_back(nuctab['A']);
                    }
                }
                ORF_block.emplace_back(torch::cat(TIS_tensor_vec, -1));
            }

            // reverse strand
            else {
                std::vector<torch::Tensor> TIS_tensor_vec;
                if (coords.first >= 16 and seqlen - coords.first >= 19) {
                    std::string window = seq.substr(coords.first - 16, 35); // flanking 16 bp up and down stream, excluding start codon
                    complement(window, window);
                    std::string start_codon = window.substr(16, 3);
                    std::swap(start_codon[0], start_codon[2]);
                    ORF_start_codon.emplace_back(start_codon); // save start codon
                    window.erase(16, 3); // remove start codon from window
                    for (int k = 0; k < 32; ++k){ // TIS is scored in 3' to 5' direction
                        char nuc = window.at(k);
                        torch::Tensor enc;
                        if (nuctab.count(nuc)){
                            enc = nuctab[nuc];
                        } else{ // non-ATGC nucleotides
                            enc = nuctab['A'];
                        }
                        TIS_tensor_vec.emplace_back(enc);
                    }
                } else {
                    // TIS fragment
                    for (int k = 0; k < 32; ++k){
                        TIS_tensor_vec.emplace_back(nuctab['A']);
                    }
                }
                ORF_block.emplace_back(torch::cat({TIS_tensor_vec}, -1));
            }

        }

        // score each TIS batch
        torch::Tensor ORF_TIS_batch = torch::cat({ORF_block}, 0);
        torch::NoGradGuard nograd;
        _TIS_model.eval();
        auto TIS_score = _TIS_model.forward({ORF_TIS_batch});
        torch::Tensor score_tensor = TIS_score.toTensor();

        // save scores
        for (int j=0; j<score_tensor.sizes()[0]; ++j){
            auto expit_score = (1 / (1 + torch::exp(-score_tensor[j]))).item<double>();
            TIS_model_scores.emplace_back(expit_score);
        }
    }

    if (verbosity) {
        progressbar.progress((int)ORF_coords.size(), (int)ORF_coords.size());
        progressbar.finish();
    }
}


void GeneFinder::translate_and_tensorize_seq() {
    // TODO might be faster to use C++17 string indexing rather than a bunch of copies
    std::string codon;
    std::pair<char, torch::Tensor> res;
    char aa;

    // reserve memory for translated and tensorized sequence // TODO does this actually help?
    seq_translated_f0.reserve(seq.length()/3);
    seq_translated_f1.reserve(seq.length()/3);
    seq_translated_f2.reserve(seq.length()/3);
    seq_translated_r0.reserve(seq.length()/3);
    seq_translated_r1.reserve(seq.length()/3);
    seq_translated_r2.reserve(seq.length()/3);

    seq_tensorized_f0.reserve(seq.length()/3);
    seq_tensorized_f1.reserve(seq.length()/3);
    seq_tensorized_f2.reserve(seq.length()/3);
    seq_tensorized_r0.reserve(seq.length()/3);
    seq_tensorized_r1.reserve(seq.length()/3);
    seq_tensorized_r2.reserve(seq.length()/3);

    for (int idx=0; idx < seq.length() - 3; idx += 3){ // TODO ignoring a few codons at the end
        // frame f0
        codon = seq.substr(idx, 3);
        res = trantab_standard[codon];
        aa = res.first;
        if (aa){
            seq_translated_f0 += aa;
            seq_tensorized_f0.emplace_back(res.second);
        } else{
            seq_translated_f0 += 'X'; // catch non-standard characters
            seq_tensorized_f0.emplace_back(torch::zeros({1, 21, 1}));
        }
        // frame r0
        std::swap(codon[0], codon[2]); // reverse codon order by swapping first and last
        complement(codon, codon);
        res = trantab_standard[codon];
        aa = res.first;
        if (aa){
            seq_translated_r0 += aa;
            seq_tensorized_r0.emplace_back(res.second);
        } else{
            seq_translated_r0 += 'X'; // catch non-standard characters
            seq_tensorized_r0.emplace_back(torch::zeros({1, 21, 1}));
        }
        // frame f1
        codon = seq.substr(idx + 1, 3);
        res = trantab_standard[codon];
        aa = res.first;
        if (aa){
            seq_translated_f1 += aa;
            seq_tensorized_f1.emplace_back(res.second);
        } else{
            seq_translated_f1 += 'X'; // catch non-standard characters
            seq_tensorized_f1.emplace_back(torch::zeros({1, 21, 1}));
        }
        // frame r1
        std::swap(codon[0], codon[2]); // reverse codon order by swapping first and last
        complement(codon, codon);
        res = trantab_standard[codon];
        aa = res.first;
        if (aa){
            seq_translated_r1 += aa;
            seq_tensorized_r1.emplace_back(res.second);
        } else{
            seq_translated_r1 += 'X'; // catch non-standard characters
            seq_tensorized_r1.emplace_back(torch::zeros({1, 21, 1}));
        }
        // frame f2
        codon = seq.substr(idx + 2, 3);
        res = trantab_standard[codon];
        aa = res.first;
        if (aa){
            seq_translated_f2 += aa;
            seq_tensorized_f2.emplace_back(res.second);
        } else{
            seq_translated_f2 += 'X'; // catch non-standard characters
            seq_tensorized_f2.emplace_back(torch::zeros({1, 21, 1}));
        }
        // frame r2
        std::swap(codon[0], codon[2]); // reverse codon order by swapping first and last
        complement(codon, codon);
        res = trantab_standard[codon];
        aa = res.first;
        if (aa){
            seq_translated_r2 += aa;
            seq_tensorized_r2.emplace_back(res.second);
        } else{
            seq_translated_r2 += 'X'; // catch non-standard characters
            seq_tensorized_r2.emplace_back(torch::zeros({1, 21, 1}));
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
            protseq += trantab_standard[nucseq.substr(j, 3)].first;
        }
        (*_gene_protseq).emplace_back(protseq);

        // get scores
        double score = ORF_scores[ORF_coords_topsort[gene_ORF_idx_topsort[i]].second];
        (*_gene_score).emplace_back(score);
    }
}
