#include <def.hpp>
#include <quantizer/IntegerQuantizer.hpp>
#include <lossless/Lossless_zstd.hpp>
#include <SZ.hpp>
#include <cstdio>
#include <iostream>
#include <memory>
#include<cmath>
#include<cstdlib>
#define radius 32768

template<class T>
inline int quantize_api(T &data, T pred, double error_bound) {

    T diff = data - pred;
    int quant_index = (int) (fabs(diff) / error_bound) + 1;
    if (quant_index < radius * 2) {
        quant_index >>= 1;
        int half_index = quant_index;
        quant_index <<= 1;
        int quant_index_shifted;
        if (diff < 0) {
            quant_index = -quant_index;
            quant_index_shifted = radius - half_index;
        } else {
            quant_index_shifted = radius + half_index;
        }
        T decompressed_data = pred + quant_index * error_bound;
        if (fabs(decompressed_data - data) > error_bound) {
            return 0;
        } else {
            data = decompressed_data;
            return quant_index_shifted;
        }
    } else {
        return 0;
    }
}

int main(int argc, char **argv) {

    std::cout << "Usage .\\sz_backend latent_vecs.dat latentformatsize quantization_bins.dat unpred_data.dat" << std::endl;

    //latent part

    size_t num_elements=0;
    auto elements=SZ::readfile<int>(argv[1], num_elements);
    std::vector<int> latent_elements(elements.get(), elements.get() + num_elements);

    int bitnum=atoi(argv[2]);
    long long radius1=(long long)(pow(2,bitnum-1));
    SZ::uchar *compressed_data_0 = new SZ::uchar[3 * num_elements * sizeof(int)];
    SZ::uchar *compressed_data_pos_0 = compressed_data;
    auto encoder_0 = SZ::HuffmanEncoder<int>();
    encoder_0.preprocess_encode(latent_elements, 4 * radius1);
    encoder_0.save(compressed_data_pos_0);
    encoder_0.encode(latent_elements, compressed_data_pos_0);
    encoder_0.postprocess_encode();
    auto lossless_0 = SZ::Lossless_zstd();
    size_t compressed_size_0;
    SZ::uchar *lossless_data_0 = lossless_0.compress(compressed_data_0,
                                                 compressed_data_pos_0 - compressed_data, compressed_size_0);
    lossless_0.postcompress_data(compressed_data_0);

    

    delete[] lossless_data_0;


    //quant part
    num_elements = 0;
    size_t num_unpred = 0;
    auto quan = SZ::readfile<int>(argv[3], num_elements);
    std::cout << "Read " << num_elements << " quantization_bins\n";
    std::vector<int> quant_inds(quan.get(), quan.get() + num_elements);

    auto unpred = SZ::readfile<int>(argv[4], num_unpred);
    std::cout << "Read " << num_unpred << " unpred_data\n";


    SZ::uchar *compressed_data = new SZ::uchar[3 * num_elements * sizeof(int)];
    SZ::uchar *compressed_data_pos = compressed_data;

    // store unpred data
    *reinterpret_cast<size_t *>(compressed_data_pos) = num_unpred;
    compressed_data_pos += sizeof(size_t);
    memcpy(compressed_data_pos, unpred.get(), num_unpred * sizeof(float));
    compressed_data_pos += num_unpred * sizeof(float);

    // huffman encoder
    auto encoder = SZ::HuffmanEncoder<int>();
    encoder.preprocess_encode(quant_inds, 4 * radius);
    encoder.save(compressed_data_pos);
    encoder.encode(quant_inds, compressed_data_pos);
    encoder.postprocess_encode();

    // lossless compressor
    auto lossless = SZ::Lossless_zstd();
    size_t compressed_size;
    SZ::uchar *lossless_data = lossless.compress(compressed_data,
                                                 compressed_data_pos - compressed_data, compressed_size);
    lossless.postcompress_data(compressed_data);

    std::cout << "Compressed size = " << (compressed_size+compressed_size_0) << std::endl;
    float ratio = 1.0 * num_elements * sizeof(float) / (compressed_size+compressed_size_0);
    std::cout << "Compression ratio = " << ratio << std::endl;

    delete[] lossless_data;

    return 0;
}
