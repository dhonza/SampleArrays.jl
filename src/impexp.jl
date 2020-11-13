using FLAC
using WAV

export normalize, normalize!, readwav, readflac, readwav_rfft, writewav

function normalize!(a::SampleArray)
    a ./= maximum(abs.(a))
    a
end

normalize(a::SampleArray) = normalize!(similar(a))

function readwav(fname; format="double")
    y, Fs, nbits, opt = wavread(fname; format=format)
    SampleArray(y, Fs)
end

function readflac(fname)
    y, Fs = load(fname)
    SampleArray(y, Fs)
end

function readwav_rfft(fname; normalized=false, format="double")
    x = readwav(fname; format=format)
    if normalized
        normalize!(x)
    end
    x, rfft(x)
end

function readflac_rfft(fname; normalized=false)
    x = readflac(fname)
    if normalized
        normalize!(x)
    end
    x, rfft(x)
end

writewav(x::SampleArray, fname::String; kwargs...) = wavwrite(data(x), fname; Fs=rate(x), kwargs...)
