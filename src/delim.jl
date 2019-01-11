using Logging
using DataFrames, CSV
using JuMP

export write_results

function write_csv(data::DataFrame, filename::AbstractString)

    @info "Writing csv file with data: " data filename

    open(filename, "w") do io
        CSV.write(io, data)
    end

end

"""
    write_results(variable, name, model_data)
"""
function write_results(variable::AbstractArray, name::AbstractString, model_data::Dict)

    edges = model_data["edges"]
    nodes = model_data["nodes"]
    commodities = model_data["commodities"]
    years = model_data["years"]

    num_nodes = length(model_data["nodes"])
    num_edges = length(edges)
    num_comm = length(model_data["commodities"])

    df = DataFrame(source=String[],
                   sink=String[],
                   commodity=String[],
                   year=Int[],
                   value=Float64[])

    filename = "$(name).csv"

    function get_values(variable::Any, i::Integer, j::Integer, comm::Integer, year::Integer)
        source = nodes[i]
        sink = nodes[j]
        comm_name = commodities[comm]
        year_name = years[year]
        value = JuMP.value(variable[i, j, comm, year])
        value != 0 ? push!(df, (source, sink, comm_name, year_name, value)) : 0
    end

    function iterate_variable(variable::Array{JuMP.VariableRef,4})
        c_indices = keys(variable)
        for index in eachindex(variable)
            (i, j, comm, year) = c_indices[index].I
            get_values(variable, i, j, comm, year)
        end
    end

    function iterate_variable(variable::JuMP.Containers.SparseAxisArray{JuMP.VariableRef,4,NTuple{4,Any}})
        for index in eachindex(variable)
            (i, j, comm, year) = index
            get_values(variable, i, j, comm, year)
        end
    end

    iterate_variable(variable)
    write_csv(df, filename)
end
