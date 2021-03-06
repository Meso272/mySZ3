#ifndef _SZ_SZ_GENERAL_HPP
#define _SZ_SZ_GENERAL_HPP

#include "predictor/Predictor.hpp"
#include "predictor/LorenzoPredictor.hpp"
#include "quantizer/Quantizer.hpp"
#include "encoder/Encoder.hpp"
#include "lossless/Lossless.hpp"
#include "utils/Iterator.hpp"
#include "utils/MemoryOps.hpp"
#include "utils/Config.hpp"
#include "utils/fileUtil.h"
#include "def.hpp"
#include <cstring>

namespace SZ {
    template<class T, size_t N, class Predictor, class Quantizer, class Encoder, class Lossless>
    class SZ_General_Compressor {
    public:


        SZ_General_Compressor(const Config<T, N> &conf,
                              Predictor predictor, Quantizer quantizer, Encoder encoder, Lossless lossless) :
                fallback_predictor(LorenzoPredictor<T, N, 1>(conf.eb)),
                predictor(predictor), quantizer(quantizer), encoder(encoder), lossless(lossless),
                block_size(conf.block_size), stride(conf.stride),
                global_dimensions(conf.dims), num_elements(conf.num) {
            static_assert(std::is_base_of_v<concepts::PredictorInterface<T, N>, Predictor>,
                          "must implement the predictor interface");
            static_assert(std::is_base_of_v<concepts::QuantizerInterface<T>, Quantizer>, "must implement the quatizer interface");
            static_assert(std::is_base_of_v<concepts::EncoderInterface<int>, Encoder>, "must implement the encoder interface");
            static_assert(std::is_base_of_v<concepts::LosslessInterface, Lossless>, "must implement the lossless interface");
        }

        // compress given the error bound
        uchar *compress(T *data, size_t &compressed_size) {
            // TODO: new quantizer if eb does not match

            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), stride, 0);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1, 0);
            std::array<size_t, N> intra_block_dims;
            std::vector<int> quant_inds(num_elements);
            predictor.precompress_data(inter_block_range->begin());
            quantizer.precompress_data();
            size_t quant_count = 0;
            struct timespec start, end;
            clock_gettime(CLOCK_REALTIME, &start);
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; ++block) {

                    // std::cout << *block << " " << lp.predict(block) << std::endl;
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1 && global_dimensions[i] - cur_index * stride < block_size) ?
                                              global_dimensions[i] - cur_index * stride : block_size;
                    }

                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());
                    concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                    if (!predictor.precompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }
                    predictor_withfallback->precompress_block_commit();
//                    quantizer.precompress_block();
                    {
                        auto intra_begin = intra_block_range->begin();
                        auto intra_end = intra_block_range->end();
                        for (auto element = intra_begin; element != intra_end; ++element) {
                            quant_inds[quant_count++] = quantizer.quantize_and_overwrite(
                                    *element, predictor_withfallback->predict(element));
                        }
                    }
                }
            }

            clock_gettime(CLOCK_REALTIME, &end);
            std::cout << "Predition & Quantization time = "
                      << (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / (double) 1000000000
                      << "s" << std::endl;

            predictor.postcompress_data(inter_block_range->begin());
            quantizer.postcompress_data();
            uchar *compressed_data;
            if (stride > block_size) {
                std::cout << "Sampling Compress Mode is ON" << std::endl;
                quant_inds.resize(quant_count);
                compressed_data = new uchar[3 * quant_count * sizeof(T)];
            } else {
                compressed_data = new uchar[2 * num_elements * sizeof(T)];
            }

            uchar *compressed_data_pos = compressed_data;
            write(global_dimensions.data(), N, compressed_data_pos);
            write(block_size, compressed_data_pos);
            predictor.save(compressed_data_pos);
            quantizer.save(compressed_data_pos);

            encoder.preprocess_encode(quant_inds, 4 * quantizer.get_radius());
            encoder.save(compressed_data_pos);
            encoder.encode(quant_inds, compressed_data_pos);
            encoder.postprocess_encode();

            uchar *lossless_data = lossless.compress(compressed_data,
                                                     compressed_data_pos - compressed_data,
                                                     compressed_size);
            lossless.postcompress_data(compressed_data);
            return lossless_data;
        }


        T *decompress(uchar const *lossless_compressed_data, const size_t length) {
            auto compressed_data = lossless.decompress(lossless_compressed_data, length);
            uchar const *compressed_data_pos = compressed_data;
            size_t remaining_length = length;
            read(global_dimensions.data(), N, compressed_data_pos, remaining_length);
            num_elements = 1;
            for (const auto &d : global_dimensions) {
                num_elements *= d;
                std::cout << d << " ";
            }
            std::cout << std::endl;
            uint block_size = 0;
            read(block_size, compressed_data_pos, remaining_length);
            predictor.load(compressed_data_pos, remaining_length);
            // std::cout << "load predictor done\n";fflush(stdout);
            quantizer.load(compressed_data_pos, remaining_length);
            // std::cout << "load quantizer done\n";fflush(stdout);
            encoder.load(compressed_data_pos, remaining_length);
            // size_t unpred_data_size = 0;
            // read(unpred_data_size, compressed_data_pos, remaining_length);
            // T const * unpred_data_pos = (T const *) compressed_data_pos;
            // compressed_data_pos += unpred_data_size * sizeof(T);
            auto quant_inds = encoder.decode(compressed_data_pos, num_elements);
            encoder.postprocess_decode();
            lossless.postdecompress_data(compressed_data);

            // std::cout << "load encoder done\n";fflush(stdout);
//    	std::cout << quant_inds[157684267] << std::endl;
            int const *quant_inds_pos = (int const *) quant_inds.data();
            std::array<size_t, N> intra_block_dims;
            auto dec_data = std::make_unique<T[]>(num_elements);
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.get(),
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), block_size,
                                                                                         0);

            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(dec_data.get(),
                                                                                         std::begin(global_dimensions),
                                                                                         std::end(global_dimensions), 1, 0);

            predictor.predecompress_data(inter_block_range->begin());
            quantizer.predecompress_data();

            std::cout << "start decompression" << std::endl;
            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                for (auto block = inter_begin; block != inter_end; block++) {
                    for (int i = 0; i < intra_block_dims.size(); i++) {
                        size_t cur_index = block.get_local_index(i);
                        size_t dims = inter_block_range->get_dimensions(i);
                        intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size
                                                                      : block_size;
                    }
                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
                    intra_block_range->set_offsets(block.get_offset());
                    intra_block_range->set_starting_position(block.get_local_index());

                    concepts::PredictorInterface<T, N> *predictor_withfallback = &predictor;
                    if (!predictor.predecompress_block(intra_block_range)) {
                        predictor_withfallback = &fallback_predictor;
                    }


//                    quantizer.predecompress_block();
                    // std::cout << "dimensions: " << intra_block_range->get_dimensions(0) << " " << intra_block_range->get_dimensions(1) << " " << intra_block_range->get_dimensions(2) << std::endl;
                    // std::cout << "index: " << block.get_local_index(0) << " " << block.get_local_index(1) << " " << block.get_local_index(2) << std::endl;
                    {
                        auto intra_begin = intra_block_range->begin();
                        auto intra_end = intra_block_range->end();
                        for (auto element = intra_begin; element != intra_end; ++element) {
                            *element = quantizer.recover(predictor_withfallback->predict(element), *(quant_inds_pos++));
                        }
                    }
                }
            }
            predictor.postdecompress_data(inter_block_range->begin());
            quantizer.postdecompress_data();
            return dec_data.release();
        }


    private:
        Predictor predictor;
        LorenzoPredictor<T, N, 1> fallback_predictor;
        Quantizer quantizer;
        Encoder encoder;
        Lossless lossless;
        uint block_size;
        uint stride;
        size_t num_elements;
        std::array<size_t, N> global_dimensions;
    };

    template<class T, uint N, class Predictor, class Quantizer, class Encoder, class Lossless>
    SZ_General_Compressor<T, N, Predictor, Quantizer, Encoder, Lossless>
    make_sz_general_compressor(const Config<T, N> &conf, Predictor predictor, Quantizer quantizer, Encoder encoder,
                               Lossless lossless) {
        return SZ_General_Compressor<T, N, Predictor, Quantizer, Encoder, Lossless>(conf, predictor, quantizer, encoder,
                                                                                    lossless);
    }
}
#endif

