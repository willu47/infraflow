module InfraFlow

include("delim.jl")

using JuMP, GLPK, Test, YAML, Logging
const MOI = JuMP.MathOptInterface

export get_data, formulate_gmcnf!, run
"""
Author: Will Usher
Date: 3rd January 2018

A linear programming new_capacity expansion generalized multi-commodity network 
flow implemented using JuMP/Julia based on generalized multi-commodity network 
flow as found in:

Ishimatsu, Takuto. “Generalized Multi-Commodity Network Flows : Case Studies 
in Space Logistics and Complex Infrastructure Systems.” 
Masssachusett Institute of Technology, 2013.
"""

"""
    get_data(file_path)

Load in model data from a YAML file
"""
function get_data(file_path::AbstractString)

    data = YAML.load(open(file_path))

    model_data = Dict{String,Any}()

    model_data["years"] = data["years"]
    years = Dict{Int,Int}([data["years"][i]=>i for i in eachindex(data["years"])])
    model_data["commodities"] = data["commodities"]
    commodities = Dict{String,Int}([data["commodities"][i]=>i for i in eachindex(data["commodities"])])

    model_data["discount_rate"] = data["discount_rate"]

    nodes = Dict{String,Int}()
    node_names = []
    self_loops = []
    for index in eachindex(data["nodes"])
        node = data["nodes"][index]
        nodes[node["name"]] = index
        push!(node_names, node["name"])
        if haskey(node, "requirements")
            push!(self_loops, index)
        end

    end
    model_data["nodes"] = node_names

    num_nodes = length(node_names)
    num_comm = length(commodities)
    num_years = length(years)
    
    cap2act = zeros(Float64, (num_nodes, num_nodes, num_comm))
    outflow_cost = zeros(Float64, (num_nodes, num_nodes, num_comm))
    inflow_cost = zeros(Float64,(num_nodes, num_nodes, num_comm))
    residual_capacity = zeros(Float64, (num_nodes, num_nodes, num_comm, num_years))
    maximum_new_capacity = fill(99999999999999999999, (num_nodes, num_nodes, num_comm, num_years))
    transformation = zeros(Float64, (num_nodes, num_nodes, num_comm, num_comm))

    requirements_outflow = zeros(Float64, (num_nodes, num_nodes, num_comm, num_comm))
    requirements_inflow = zeros(Float64, (num_nodes, num_nodes, num_comm, num_comm))

    edges = []
    for edge in data["edges"]
        source = nodes[edge["source"]]
        sink = nodes[edge["sink"]]
        push!(edges, source, sink)

        for flow in edge["flow"]
            comm_idx = commodities[flow["name"]]
            requirements_inflow[source, sink, comm_idx, comm_idx] = 1
            requirements_outflow[source, sink, comm_idx, comm_idx] = 1
        end

        if haskey(edge, "cap2act")
            for (comm, value) in edge["cap2act"]
                comm_idx = commodities[comm]
                cap2act[source, sink, comm_idx] = value
            end
        end

        if haskey(edge, "operational_cost")
            for (comm, value) in edge["operational_cost"]
                comm_idx = commodities[comm]
                outflow_cost[source, sink, comm_idx] = value
            end
        end

        if haskey(edge, "losses")
            for (comm, value) in edge["losses"]
                comm_idx = commodities[comm]
                transformation[source, sink, comm_idx, comm_idx] = (1.0 - value)
            end
        end

    end

    for loop in self_loops
        push!(edges, loop, loop)
    end

    edges = transpose(reshape(edges, 2, :))
    model_data["edges"] = edges

    demand = zeros((num_nodes, num_comm, num_years))
    capacity_cost = zeros((num_nodes, num_nodes, num_comm))

    for index in eachindex(data["nodes"])
        node = data["nodes"][index]
        source = node["name"]
        node_idx = nodes[source]
        if haskey(node, "demand")
            for (commodity, yearly_data) in node["demand"]
                for (year, value) in yearly_data
                    node_idx = nodes[source]
                    comm_index = commodities[commodity]
                    year_idx = years[year]
                    demand[node_idx, comm_index, year_idx] = value
                end
            end
        end

        if haskey(node, "investment_cost")
            output_commodity = node["output"]
            output_comm_idx = commodities[output_commodity]
            capacity_cost[node_idx, node_idx, output_comm_idx] = node["investment_cost"]
        end
    
        if haskey(node, "cap2act")
            output_commodity = node["output"]
            output_comm_idx = commodities[output_commodity]
            cap2act[node_idx, node_idx, output_comm_idx] = node["cap2act"]
        end

        if haskey(node, "requirements")
            for requirement in node["requirements"]
                comm = node["output"]
                comm_idx = commodities[comm]
                r_comm = requirement["name"]
                r_comm_idx = commodities[r_comm]            
                cap2act[node_idx, node_idx, r_comm_idx] = 1
                transformation[node_idx, node_idx, comm_idx, r_comm_idx] = 1
                requirements_outflow[node_idx, node_idx, r_comm_idx, comm_idx] = requirement["value"]
                requirements_inflow[node_idx, node_idx, comm_idx, comm_idx] = 1
                requirements_outflow[node_idx, node_idx, comm_idx, comm_idx] = 1
                if haskey(node, "residual_capacity")
                    for (year, value) in node["residual_capacity"]
                        year_idx = years[year]
                        residual_capacity[node_idx, node_idx, comm_idx, year_idx] = value
                    end
                end
                if haskey(node, "maximum_new_capacity")
                    for (year, value) in node["maximum_new_capacity"]
                        year_idx = years[year]
                        maximum_new_capacity[node_idx, node_idx, comm_idx, year_idx] = value
                    end
                end                
            end

        end


    end

    model_data["demand"] = demand
    model_data["capacity_cost"] = capacity_cost
    model_data["cap2act"] = cap2act
    model_data["outflow_cost"] = outflow_cost
    model_data["inflow_cost"] = inflow_cost
    model_data["residual_capacity"] = residual_capacity
    model_data["maximum_new_capacity"] = maximum_new_capacity
    model_data["transformation"] = transformation
    model_data["requirements_inflow"] = requirements_inflow
    model_data["requirements_outflow"] = requirements_outflow

    source_nodes = view(edges, :, 1)
    sink_nodes = view(edges, :, 2)

    outflow_edges = make_edge_dict(source_nodes, sink_nodes)
    inflow_edges = make_edge_dict(sink_nodes, source_nodes)

    model_data["outflow_edges"] = outflow_edges
    model_data["inflow_edges"] = inflow_edges

    return model_data
end


"""
    make_edge_dict(edge_nodes, other_nodes)

"""
function make_edge_dict(edge_nodes::AbstractArray, other_nodes::AbstractArray)
    edges = Dict{Int8,Array{Int8}}()
    for node in Set(edge_nodes)
        for i in eachindex(edge_nodes)
            if edge_nodes[i] == node
                if haskey(edges, node)
                    push!(edges[node], other_nodes[i])
                else
                    push!(edges, node=>[other_nodes[i]])
                end
            end
        end
    end
    return edges
end

"""
    formulate_gmcnf()

``\\sum_{ijky} c_{ijky}x_{ijky}`` 

"""
function formulate_gmcnf!(model_data::Dict{String,Any}, model; verbose = true)
    
    nodes = model_data["nodes"]
    edges = model_data["edges"]
    commodities = model_data["commodities"]
    discount_rate = model_data["discount_rate"]
    years = model_data["years"]

    outflow_edges = model_data["outflow_edges"]
    inflow_edges = model_data["inflow_edges"]

    num_nodes = length(nodes)
    num_edges = length(edges)
    num_comm = length(commodities)
    num_years = length(years)
    
    # demand at node by commodity
    demand = model_data["demand"]

    capacity_cost = model_data["capacity_cost"]
    cap2act = model_data["cap2act"]

    inflow_cost = model_data["inflow_cost"]
    outflow_cost = model_data["outflow_cost"]
    
    # describe the commodity requirements for an in- or out-flow
    requirements_outflow = model_data["requirements_outflow"]
    requirements_inflow = model_data["requirements_inflow"]
    
    # describe the flow gain/loss or transformation between commodities
    transformation = model_data["transformation"]
    
    # upper bound on operational decision variables (new_capacity)
    residual_capacity = model_data["residual_capacity"]

    # upper bound on new capacity
    maximum_new_capacity = model_data["maximum_new_capacity"]

    new_capacity = @variable(model,
                         new_capacity[i=1:num_nodes, j=1:num_nodes, k=1:num_comm, y=1:num_years],
                         lower_bound = 0)

    outflow = @variable(model, 
                        outflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm, y=1:num_years], 
                        lower_bound = 0)
    inflow = @variable(model, 
                       inflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm, y=1:num_years], 
                       lower_bound = 0)
    
    @constraint(
        model,
        max_new_capacity[i=1:num_nodes, j=1:num_nodes, k=1:num_comm, y=1:num_years],
        new_capacity[i, j, k, y] <= maximum_new_capacity[i, j, k, y]
    )

    @variable(
        model,
        total_annual_capacity[i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        lower_bound = 0
    )

    @constraint(
        model,
        accumulate_capacity[i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        total_annual_capacity[i, j, k, y] == residual_capacity[i, j, k, y] + sum(new_capacity[i, j, k, z] for z in 1:y)
    )

    @constraint(
        model,
        capacity_exp_outflow[i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        outflow[i, j, k, y] <= total_annual_capacity[i, j, k, y] * cap2act[i, j, k])

    @constraint(
        model,
        capacity_exp_inflow[i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        inflow[i, j, k, y] <= total_annual_capacity[i, j, k, y] * cap2act[i, j, k])
        
    function discount_factor(year::Int)
        return ((1 + discount_rate) ^ (years[year] - years[1]))
    end

    discounted_capital_cost = @expression(
        model,
        [i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        new_capacity[i, j, k, y] * capacity_cost[i, j, k] / discount_factor(y)
    )

    discounted_operational_cost = @expression(
        model,
        [i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        outflow_cost[i, j, k] * outflow[i, j, k, y] 
        + inflow_cost[i, j, k] * inflow[i, j, k, y] 
        / discount_factor(y)
    )

    @objective(
        model, 
        Min, 
        sum(discounted_operational_cost[i, j, k, y]
            + discounted_capital_cost[i, j, k, y]
            for i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm, y=1:num_years)
        )


    requirements_outflow_const = @expression(
        model,
        requirements_outflow_const[i=1:num_nodes, k=1:num_comm, y=1:num_years],
        if haskey(outflow_edges, i)
            sum(
                sum(requirements_outflow[i, j, k, l] for l in 1:num_comm) * outflow[i, j, k, y]
                for j in outflow_edges[i]
            )
                
        else
            @info "No outflow edge found for ($(nodes[i]), $(commodities[k]), $(years[y]))"
            0
        end
    )

    requirements_inflow_const = @expression(
        model,
        requirements_inflow_const[i=1:num_nodes, k=1:num_comm, y=1:num_years],
        if haskey(inflow_edges, i)
            sum(
                sum(requirements_inflow[h, i, k, m] for m in 1:num_comm) * inflow[h, i, k, y]
                for h in inflow_edges[i]
            )

        else
            @info "No inflow edge found for ($(nodes[i]), $(commodities[k]), $(years[y]))"
            0
        end
    )

    mass_balance = @constraint(
        model,
        mass_balance[i=1:num_nodes, k=1:num_comm, y=1:num_years],
        requirements_outflow_const[i, k, y] - requirements_inflow_const[i, k, y]
        <= demand[i, k, y])

    flow_transformation = @constraint(
        model,
        flow_transformation[i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm, l in 1:num_comm, y=1:num_years; transformation[i, j, k, l] != 0],
        transformation[i, j, k, l] * outflow[i, j, l, y] == inflow[i, j, k, y]
        )

    return model
end

function print_vars(var_object)
    for var in var_object
        if JuMP.value(var) != 0.0
            println("$(var): $(JuMP.value(var))")
        end
    end
end

function print_duals(con_object)
    for con in con_object
        if JuMP.shadow_price(con) != 0.0
            println("$(con): $(JuMP.shadow_price(con))")
        end
    end
end

function print_constraint(con_object)
    for con in con_object
        println("$(con)")
    end
end

function run(file_path::String)

    model_data = get_data(file_path)

    model = Model(with_optimizer(GLPK.Optimizer))

    @time model = formulate_gmcnf!(model_data, model, verbose = true)

    println("Compiled model, now running")
    @time JuMP.optimize!(model)
    println("Finished running, objective: £$(JuMP.objective_value(model))")

    status = JuMP.termination_status(model) == MOI.OPTIMAL
    feasible = JuMP.primal_status(model) == MOI.FEASIBLE_POINT

    if status
        println("\nSolution optimal\n")
        outflow = model[:outflow]
        inflow = model[:inflow]
        total_capacity = model[:total_annual_capacity]
        new_capacity = model[:new_capacity]

        write_results(outflow, "outflow", model_data)
        write_results(inflow, "inflow", model_data)
        write_results(total_capacity, "total_capacity", model_data)
        write_results(new_capacity, "new_capacity", model_data)
        
    else
        println("\nSolution not optimal\n")
        println("i,j=> $(model_data["nodes"])")
        println("k=> $(model_data["commodities"])")
        println("y=> $(model_data["years"])")
        println("mass_balance[i, k, y]")
        print_constraint(model[:mass_balance])
        println("flow_transformation[i, j, k, y]")
        print_constraint(model[:flow_transformation])
        print_constraint(model[:max_new_capacity])
        print_constraint(model[:accumulate_capacity])
    end

end

if length(ARGS) > 0
    run(ARGS[1])
end

end
