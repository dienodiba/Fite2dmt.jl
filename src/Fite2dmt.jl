# Dieno Diba 2024
# FITE2DMT
# Main module

module Fite2dmt

using MKL
using DelimitedFiles
using Statistics
using SparseArrays
using LinearAlgebra
using LinearInterpolations
using Printf
using Dates
using KLU
using StaticArrays
using Base.Threads

export InversionRhoaPhaseTipper
export InversionRhoaPhase
export InversionTipper

include("CalculateFieldAtSideBoundaryTE.jl")
include("CalculateFieldAtSideBoundaryTM.jl")
include("CalculateModelUpdateDataSpace.jl")
include("CalculateModelUpdateModelSpace.jl")
include("CalculateObjectiveFunctionRhoaPhase.jl")
include("CalculateObjectiveFunctionRhoaPhaseTipper.jl")
include("CalculateObjectiveFunctionTipper.jl")
include("CalculateResponseFunctionsTERhoaPhase.jl")
include("CalculateResponseFunctionsTERhoaPhaseTipper.jl")
include("CalculateResponseFunctionsTETipper.jl")
include("CalculateResponseFunctionsTMRhoaPhase.jl")
include("ConstructJacobianTERhoaPhase.jl")
include("ConstructJacobianTERhoaPhaseTipper.jl")
include("ConstructJacobianTETipper.jl")
include("ConstructJacobianTMRhoaPhase.jl")
include("ConstructRoughnessMatrix.jl")
include("ConstructStiffnessMatrixTE.jl")
include("ConstructStiffnessMatrixTM.jl")
include("ExportFinalModelResistivity.jl")
include("ExportFinalModelRhoaPhase.jl")
include("ExportFinalModelRhoaPhaseTipper.jl")
include("ExportFinalModelTipper.jl")
include("ExportLogFileRhoaPhase.jl")
include("ExportLogFileRhoaPhaseTipper.jl")
include("ExportLogFileTipper.jl")
include("FindBoundaryElementsTE.jl")
include("FindBoundaryElementsTM.jl")
include("FindBoundaryNodesTE.jl")
include("FindBoundaryNodesTM.jl")
include("FindNodesBeneathStationsTE.jl")
include("FindNodesBeneathStationsTM.jl")
include("InversionRhoaPhase.jl")
include("InversionRhoaPhaseTipper.jl")
include("InversionTipper.jl")
include("MergeJacobianRhoaPhase.jl")
include("MergeJacobianRhoaPhaseTipper.jl")
include("MergeJacobianTipper.jl")
include("ReadInitialModel.jl")
include("ReadInversionSettingParametersRhoaPhase.jl")
include("ReadInversionSettingParametersRhoaPhaseTipper.jl")
include("ReadInversionSettingParametersTipper.jl")
include("ReadObservedDataRhoaPhase.jl")
include("ReadObservedDataRhoaPhaseTipper.jl")
include("ReadObservedDataTipper.jl")

end
