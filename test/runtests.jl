using Fite2dmt
using Test

@testset "Fite2dmt.jl" begin
    @test Fite2dmt.InversionRhoaPhaseTipper("Model01_data_RhoaPhaseTipper.txt","Model01_m0.txt","Model01_topo.txt","Model01_setting_RhoaPhaseTipper.txt") == nothing
    @test Fite2dmt.InversionRhoaPhase("Model01_data_RhoaPhase.txt","Model01_m0.txt","Model01_topo.txt","Model01_setting_RhoaPhase.txt") == nothing
    @test Fite2dmt.InversionTipper("Model01_data_Tipper.txt","Model01_m0.txt","Model01_topo.txt","Model01_setting_Tipper.txt") == nothing
end
