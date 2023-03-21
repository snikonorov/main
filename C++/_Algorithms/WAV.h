#ifndef WAV_H_INCLUDED
#define WAV_H_INCLUDED

#include <vector>
#include <string>

namespace
{
    constexpr unsigned long as_unsigned(const void* s, const unsigned n = 4)
    {
        ///return *(unsigned*)s;

        /**return (unsigned char)s[0] +
              ((unsigned char)s[1] << 8) +
              ((unsigned char)s[2] << 16) +
              ((unsigned char)s[3] << 24);*/

        unsigned long res = 0;

        for(int k = int(n) - 1; k > 0; k--)
        {
            res = (res + ((unsigned char*)s)[k]) << 8;
        }

        return res + ((unsigned char*)s)[0];
    }

    long as_integer(const void* s, const unsigned n = 4)
    {
        switch(n)
        {
            case 1:
                return *(char*)s;

            case 2:
                return *(short*)s;

            case 4:
                return *(long*)s;

            case 8:
                return *(long long*)s;
        }

        return 0;
    }

    /*inline void print(char (&s)[4])
    {
        for(unsigned k = 0; k < 4; k++)
        {
            putchar(s[k]);
        }
    }*/
}

enum class WAVFormat
{
    Signed16,
    Signed32,
    Float32,
    Float64
};

/// multi-channel sampled signal
///
struct TSampledSignal: public std::vector<std::vector<double>>
{
    unsigned sample_rate;

    TSampledSignal(unsigned N = 0): std::vector<std::vector<double>>(N)
    {
        sample_rate = 0;
    }

    TSampledSignal(const std::vector<std::vector<double>>& data, unsigned sample_rate)
                  : std::vector<std::vector<double>>(data), sample_rate(sample_rate)
    {
    }

    double GetDuration() const
    {
        return size() ? 1.*(this->at(0).size() - 1)/sample_rate : 0.0;
    }
};

///-----------------------------------------------------------------------------------------

std::map<std::string, std::map<std::string, unsigned>> GetWAVMetadata(const char* path)
{
    struct TChunkHeader
    {
        char id[4];
        unsigned size;
    };

    struct TWAVHeader
    {
        char RIFF[4];   /// "RIFF"
        unsigned n;     /// <chunk size>
        char WAVE[4];   /// "WAVE"
    };

    ///--------------------------------------------------------

    std::map<std::string, std::map<std::string, unsigned>> res;

    FILE* f = fopen(path, "rb");

    if(f)
    {
        TChunkHeader header;
        std::vector<unsigned char> chunk_data;

        TWAVHeader WAVHeader;

        fread(&WAVHeader, sizeof(TWAVHeader), 1, f);

        res["RIFF"]["size"] = WAVHeader.n;

        if
        (
            as_unsigned(WAVHeader.RIFF) == as_unsigned("RIFF") &&
            as_unsigned(WAVHeader.WAVE) == as_unsigned("WAVE")
        )
        {
            while(fread(&header, sizeof(TChunkHeader), 1, f))
            {
                auto ID = as_unsigned(header.id);

                if(ID == as_unsigned("fmt "))
                {
                    auto& fmt = res["fmt "];

                    fmt["size"] = header.size;

                    chunk_data.resize(header.size);
                    fread((unsigned char*) chunk_data.data(), sizeof(unsigned char), header.size, f);

                    fmt["format_code"]     = as_unsigned(chunk_data.data(),      2);
                    fmt["channel_number"]  = as_unsigned(chunk_data.data() + 2,  2);
                    fmt["sample_rate"]     = as_unsigned(chunk_data.data() + 4,  4);
                    fmt["data_rate"]       = as_unsigned(chunk_data.data() + 8,  4);
                    fmt["block_align"]     = as_unsigned(chunk_data.data() + 12, 2);
                    fmt["bits_per_sample"] = as_unsigned(chunk_data.data() + 14, 2);

                    if(header.size > 16)
                    {
                        auto extension_size = as_unsigned(chunk_data.data() + 16, 2);

                        fmt["extension_size"] = extension_size;

                        if(extension_size)
                        {
                            fmt["ex_format_code"] = as_unsigned(chunk_data.data() + 24, 2);
                        }
                    }
                }
                else
                {
                    std::string s;
                    s.resize(4);

                    for(unsigned k = 0; k < s.size(); k++)
                    {
                        s[k] = header.id[k];
                    }

                    res[s]["size"] = header.size;

                    fseek(f, header.size, SEEK_CUR);
                }
            }
        }

        fclose(f);
    }

    return res;
}

TSampledSignal ReadWAV(FILE* f, bool normalize = true, bool close_handle = false)
{
    struct TChunkHeader
    {
        char id[4];
        unsigned size;
    };

    struct TWAVHeader
    {
        char RIFF[4];   /// "RIFF"
        unsigned n;     /// <chunk size>
        char WAVE[4];   /// "WAVE"
    };

    enum struct FormatCode
    {
        None = 0,
        PCM = 1,
        Float = 3
    };

    ///--------------------------------------------------------

    TSampledSignal res;

    FormatCode format = FormatCode::None;
    unsigned bytes_per_sample = 0;

    if(f)
    {
        TChunkHeader header;
        std::vector<unsigned char> chunk_data;

        TWAVHeader WAVHeader;

        fread(&WAVHeader, sizeof(TWAVHeader), 1, f);

        ///println(WAVHeader.n);

        if
        (
            as_unsigned(WAVHeader.RIFF) == as_unsigned("RIFF") &&
            as_unsigned(WAVHeader.WAVE) == as_unsigned("WAVE")
        )
        {
            while(fread(&header, sizeof(TChunkHeader), 1, f))
            {
                auto ID = as_unsigned(header.id);

                if(ID == as_unsigned("fmt "))
                {
                    if(header.size < 16)
                    {
                        println("Invalid [fmt] chunk size (", header.size, " bytes)");
                        break;
                    }

                    chunk_data.resize(header.size);
                    fread((unsigned char*) chunk_data.data(), sizeof(unsigned char), header.size, f);

                    auto format_code     = as_unsigned(chunk_data.data(),      2);
                    auto channel_number  = as_unsigned(chunk_data.data() + 2,  2);
                    auto sample_rate     = as_unsigned(chunk_data.data() + 4,  4);
                    //auto data_rate     = as_unsigned(chunk_data.data() + 8,  4);
                    //auto block_align   = as_unsigned(chunk_data.data() + 12,  2);
                    auto bits_per_sample = as_unsigned(chunk_data.data() + 14, 2);

                    if(header.size > 16)
                    {
                        auto extension_size = as_unsigned(chunk_data.data() + 16, 2);

                        if(extension_size)
                        {
                            format_code = as_unsigned(chunk_data.data() + 24, 2);
                        }
                    }

                    switch(format_code)
                    {
                        case 1:
                            format = FormatCode::PCM;
                            break;

                        case 3:
                            format = FormatCode::Float;
                            break;
                    }

                    res.resize(channel_number);
                    res.sample_rate = sample_rate;

                    bytes_per_sample = bits_per_sample/8;

                    //println("format_code:     ", format_code);
                    //println("channel_number:  ", channel_number);
                    //println("sample_rate:     ", sample_rate);
                    //println("bits_per_sample: ", bits_per_sample);
                }
                else if(ID == as_unsigned("data"))
                {
                    chunk_data.resize(header.size);
                    fread((unsigned char*) chunk_data.data(), sizeof(unsigned char), header.size, f);

                    /// number of samples
                    unsigned N = chunk_data.size()/res.size()/bytes_per_sample;

                    for(unsigned k = 0; k < res.size(); k++)
                    {
                        res[k].resize(N);
                    }

                    unsigned k = 0;
                    long long B = (1 << (8*bytes_per_sample - 2));

                    long long _min = -B*2;
                    //long long _max = (B-1)*2 + 1;

                    unsigned long long A = 2*(2*B - 1) + 1;

                    for(unsigned n = 0; n < N; n++)
                    {
                        for(unsigned channel = 0; channel < res.size(); channel++)
                        {
                            switch(format)
                            {
                                case FormatCode::PCM:
                                    {
                                        res[channel][n] = as_integer(chunk_data.data() + k, bytes_per_sample);

                                        if(normalize)
                                            res[channel][n] = 2.*(res[channel][n] - _min)/A - 1;
                                    }
                                    break;

                                case FormatCode::Float:

                                    if(bytes_per_sample == 4)
                                    {
                                        res[channel][n] = *(float*)(chunk_data.data() + k);
                                    }
                                    else if(bytes_per_sample == 8)
                                    {
                                        res[channel][n] = *(double*)(chunk_data.data() + k);
                                    }

                                    break;
                            }

                            k += bytes_per_sample;
                        }
                    }

                    break;
                }
                else
                {
                    //print("[chunk: ");
                    //print(header.id);
                    //println(" (size: ", header.size, ")]");

                    fseek(f, header.size, SEEK_CUR);
                }
            }
        }

        if(close_handle)
            fclose(f);
    }

    return res;
}

inline TSampledSignal ReadWAV(const char* path, bool normalize = true)
{
    return ReadWAV(fopen(path, "rb"), normalize, true);
}

/*inline TSampledSignal ReadWAV(const wchar_t* path, bool normalize = true)
{
    return ReadWAV(_wfopen(path, L"rb"), normalize, true);
}*/

bool WriteWAV
(
    FILE* f,
    const TSampledSignal& data,
    WAVFormat format = WAVFormat::Signed16,
    bool denormalize = true,
    bool close_handle = false
)
{
    auto write_value = [](auto a, FILE* f)
    {
        fwrite(&a, sizeof(a), 1, f);
    };

    auto clamp = [](double _min = -1.0, double _max = 1.0)
    {
        return [_min, _max](double x)
        {
            return x < _min ? _min : x > _max ? _max : x;
        };
    };

    if(data.size() && data[0].size())
    {
        if(f)
        {
            ///short bits_per_sample = 0;
            short bytes_per_sample = 0;

            switch(format)
            {
                case WAVFormat::Signed16:
                    ///bits_per_sample = 16;
                    bytes_per_sample = 2;
                    break;

                case WAVFormat::Signed32:
                case WAVFormat::Float32:
                    ///bits_per_sample = 32;
                    bytes_per_sample = 4;
                    break;

                case WAVFormat::Float64:
                    ///bits_per_sample = 64;
                    bytes_per_sample = 8;
                    break;
            }

            unsigned data_size = data.size()*data[0].size()*bytes_per_sample;
            bool format_extensible = (format == WAVFormat::Signed32) || (data.size() > 2);

            ///-------------------------------------------------------

            fwrite("RIFF", sizeof(char), 4, f);

            if(format_extensible)
                write_value(60u + data_size, f);         /// file_size - 8
            else
                write_value(38u + data_size, f);         /// file_size - 8

            fwrite("WAVE", sizeof(char), 4, f);

            ///-------------------------------------------------------

            short block_align = bytes_per_sample*data.size();
            short format_code = (format == WAVFormat::Signed16 || format == WAVFormat::Signed32) ? 1 : 3;

            fwrite("fmt ", sizeof(char), 4, f);

            if(format_extensible)
                write_value(40u, f);
            else
                write_value(18u, f);

            if(format_extensible)
                write_value(short(0xFFFE), f);
            else
                write_value(format_code, f);

            write_value(short(data.size()), f);                 /// number of channels
            write_value(data.sample_rate, f);                   /// sampling rate
            write_value(data.sample_rate*block_align, f);       /// byte rate
            write_value(block_align, f);                        /// block align
            write_value(short(8*bytes_per_sample), f);          /// bits per sample

            if(format_extensible)
            {
                write_value(short(22), f);

                write_value(short(8*bytes_per_sample), f);      /// valid bits per sample
                write_value(0u, f);                             /// channel mask
                write_value(format_code, f);                    /// format code

                fwrite("\x00\x00\x00\x00\x10\x00\x80\x00\x00\xAA\x00\x38\x9B\x71", sizeof(char), 14, f);
            }
            else
                write_value(short(0), f);

            ///-------------------------------------------------------

            fwrite("data", sizeof(char), 4, f);
            write_value(data_size, f);

            /// number of samples
            const unsigned& N = data[0].size();

            unsigned k = 0;
            long long B = (1 << (8*bytes_per_sample - 2));

            long long _min = -B*2;
            ///long long _max = (B-1)*2 + 1;

            unsigned long long A = 2*(2*B - 1) + 1;

            auto C = clamp();

            for(unsigned n = 0; n < N; n++)
            {
                for(unsigned channel = 0; channel < data.size(); channel++)
                {
                    double res = C(data[channel][n]);
                    ///double res = data[channel][n];

                    switch(format)
                    {
                        case WAVFormat::Signed16:
                        case WAVFormat::Signed32:

                            if(denormalize)
                                res = round(A*(res + 1)/2 + _min);

                            if(bytes_per_sample == 2)
                            {
                                write_value(short(res), f);
                            }
                            else //if(bytes_per_sample == 4)
                            {
                                write_value(long(res), f);
                            }

                            break;

                        case WAVFormat::Float32:
                        case WAVFormat::Float64:

                            if(bytes_per_sample == 4)
                            {
                                write_value(float(res), f);
                            }
                            else //if(bytes_per_sample == 8)
                            {
                                write_value(res, f);
                            }

                            break;
                    }
                }
            }

            if(close_handle)
                fclose(f);
        }
        else
            return false;
    }
    else
        return false;

    return true;
}

bool WriteWAV
(
    const char* path,
    const TSampledSignal& data,
    WAVFormat format = WAVFormat::Signed16,
    bool denormalize = true
)
{
    return WriteWAV(fopen(path, "wb"), data, format, denormalize, true);
}

/*bool WriteWAV
(
    const wchar_t* path,
    const TSampledSignal& data,
    WAVFormat format = WAVFormat::Signed16,
    bool denormalize = true
)
{
    return WriteWAV(_wfopen(path, L"wb"), data, format, denormalize, true);
}*/

#endif // WAV_H_INCLUDED
