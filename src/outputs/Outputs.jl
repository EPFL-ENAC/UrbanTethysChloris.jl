module Outputs

using Dates
using ..UrbanTethysChloris: AbstractModel

include("output_manager.jl")
include("read_results.jl")

export OutputManager, initialize_outputs, save_outputs!, finalize_outputs!, unpack_results

end
