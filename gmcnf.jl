using JuMP, GLPK, Test
const MOI = JuMP.MathOptInterface

"""
Author: Will Usher
Date: 3rd January 2018

JuMP implementation of generalized multi-commodity network flow as found in 
Ishimatsu, Takuto. “Generalized Multi-Commodity Network Flows : Case Studies in 
Space Logistics and Complex Infrastructure Systems.” 
Masssachusett Institute of Technology, 2013.
"""

function make_edge_dict(edge_nodes, other_nodes)
    edges = Dict()
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


function gmcnf(; verbose = true)
    nodes = ["household", "power_station", "water_plant", "diesel_resource"]
    edges = [2 3; 2 1; 3 1; 4 2;]
    commodities = ["water", "electricity", "diesel"]
     
    source_nodes = view(edges, :, 1)
    sink_nodes = view(edges, :, 2)

    outflow_edges = make_edge_dict(source_nodes, sink_nodes)
    @test outflow_edges == Dict{Any,Any}(4=>[2],2=>[3, 1],3=>[1])
    inflow_edges = make_edge_dict(sink_nodes, source_nodes)
    @test inflow_edges == Dict{Any,Any}(3=>[2],1=>[2, 3], 2=>[4])

    num_nodes = length(nodes)
    num_edges = length(edges)
    num_comm = length(commodities)
    
    # demand at node by commodity
    demand = zeros((num_nodes, num_comm))
    demand[1, 1] = 3
    demand[1, 2] = 5

    inflow_cost = zeros((num_nodes, num_nodes, num_comm))
    outflow_cost = zeros((num_nodes, num_nodes, num_comm))

    requirements_outflow = zeros((num_nodes, num_nodes, num_comm, num_comm))
    requirements_inflow = zeros((num_nodes, num_nodes, num_comm, num_comm))
    transformation = zeros((num_nodes, num_nodes, num_comm, num_comm))
    concurrent_outflow = zeros((num_nodes, num_nodes, num_comm))
    concurrent_inflow = zeros((num_nodes, num_nodes, num_comm))

    flow_bounds = zeros((num_nodes, num_nodes, num_comm))

    model = Model(with_optimizer(GLPK.Optimizer))

    outflow = @variable(model, 
                        [i=1:num_nodes, j=1:num_nodes, k=1:num_comm], 
                        lower_bound = 0, 
                        upper_bound = flow_bounds[i, j, k],
                        base_name="outflow")
    inflow = @variable(model, 
                       inflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm], 
                       lower_bound = 0,
                       upper_bound = flow_bounds[i, j, k],
                       base_name="inflow")
    
    @objective(
        model, 
        Min, 
        sum(outflow_cost[i, j, k] * outflow[i, j, k] 
            + inflow_cost[i, j, k] * inflow[i, j, k] 
            for i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm)
        )
        

    requirements_outflow_const = @expression(model,
                [i=1:num_nodes, k=1:num_comm],
                if haskey(outflow_edges, i)
                    sum(sum(requirements_outflow[i, j, k, l] for l in 1:num_comm) * outflow[i, j, k] for j in outflow_edges[i])
                else
                    0
                end)

    requirements_inflow_const = @expression(model,
                [i=1:num_nodes, k=1:num_comm],
                if haskey(inflow_edges, i)
                    sum(sum(requirements_inflow[h, i, m, k] for m in 1:num_comm) * inflow[h, i, k] for h in inflow_edges[i])
                else
                    0
                end)

    mass_balance = @constraint(model,
                               [i=1:num_nodes, k=1:num_comm],
                               requirements_outflow_const[i, k] - requirements_inflow_const[i, k]
                               <= demand[i, k])

    flow_transformation = @constraint(model,
                [i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm],
                sum(transformation[i, j, k, l] * outflow[i, j, l] for l in 1:num_comm) 
                == inflow[i, j, k])

    self_constraint_outflow = @constraint(model,
                [i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm],
                concurrent_outflow[i, j, k] * outflow[i, j, k] <= 0)

    self_constraint_inflow = @constraint(model,
                [i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm],
                concurrent_inflow[i, j, k] * inflow[i, j, k] <= 0)

    println("Compiled model, now running")
    JuMP.optimize!(model)
    println("Finished running, objective: $(JuMP.objective_value(model))")

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    # @test JuMP.objective_value(model) == 225700.0

    print(mass_balance)
    @show model

end

@time gmcnf(verbose = true)
# @time gmcnf(verbose = true)
# @time gmcnf(verbose = true)