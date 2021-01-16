using CSVFiles
using DataFrames
using FileIO
using FLAC
using WAV

export normalize, normalize!, readwav, readflac, readwav_rfft, writewav, readfrd

function normalize!(a::SampleArray)
    a ./= maximum(abs.(a))
    a
end

normalize(a::SampleArray) = normalize!(similar(a))

function readwav(fname; format="double", names=nothing)
    y, Fs, nbits, opt = wavread(fname; format=format)
    SampleArray(y, Fs, names)
end

function readflac(fname; names=nothing)
    y, Fs = load(fname)
    SampleArray(y, Fs, names)
end

function readwav_rfft(fname; normalized=false, format="double", names=nothing)
    x = readwav(fname; format=format, names=names)
    if normalized
        normalize!(x)
    end
    x, rfft(x)
end

function readflac_rfft(fname; normalized=false, names=nothing)
    x = readflac(fname; names=names)
    if normalized
        normalize!(x)
    end
    x, rfft(x)
end

writewav(x::SampleArray, fname::String; kwargs...) = wavwrite(data(x), fname; Fs=rate(x), kwargs...)

function readfrd(fname, rate; name=nothing, magf=db2amp, phasef=deg2rad)
    df = DataFrame(load(File(format"TSV", fname); header_exists=false, colnames=[:freq, :mag, :phase], spacedelim=true))
    # df = DataFrames.readtable(fname, header=false, colnames=[:freq, :mag, :phase])
    sort!(df, [:freq])
    df = df[df.freq .<= toHz(rate)/2, :]
    SpectrumArray(Complex.(MagPhase.(magf.(df.mag), phasef.(df.phase))), rate, df.freq, isnothing(name) ? nothing : [name])
    # SpectrumArray(MagPhase.(magf.(df.mag), phasef.(df.phase)), rate, df.freq, isnothing(name) ? nothing : [name])
end