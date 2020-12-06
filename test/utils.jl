using Unitful

@testset "utils" begin
    r = toHz(12.3)
    @test typeof(r) == Float64
    r = toHz(12.3u"kHz")
    @test typeof(r) == Float64
    @test r == 12300

    r = tos(12.3)
    @test typeof(r) == Float64
    r = tos(12.3u"s")
    @test typeof(r) == Float64
    @test r == 12.3
end