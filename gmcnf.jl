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


function formulate_gmcnf(; verbose = true)
    nodes = ["household", "power_station", "water_plant", "diesel_resource", "water_resource"]
    edges = [2 3; 2 1; 3 1; 4 2; 5 3; 2 2; 3 3; 4 4; 5 5;]
    commodities = ["water", "electricity", "diesel"]
     
    source_nodes = view(edges, :, 1)
    sink_nodes = view(edges, :, 2)

    outflow_edges = make_edge_dict(source_nodes, sink_nodes)
    @test outflow_edges == Dict{Int8, Array{Int8}}(4=>[2, 4], 2=>[3, 1, 2], 3=>[1, 3], 5=>[3, 5])
    inflow_edges = make_edge_dict(sink_nodes, source_nodes)
    @test inflow_edges == Dict{Int8, Array{Int8}}(3=>[2, 5, 3], 1=>[2, 3], 2=>[4, 2], 5=>[5], 4=>[4])

    num_nodes = length(nodes)
    num_edges = length(edges)
    num_comm = length(commodities)
    
    # demand at node by commodity
    demand = zeros((num_nodes, num_comm))

    inflow_cost = zeros((num_nodes, num_nodes, num_comm))
    outflow_cost = zeros((num_nodes, num_nodes, num_comm))
    
    # describe the commodity requirements for an in- or out-flow
    requirements_outflow = zeros((num_nodes, num_nodes, num_comm, num_comm))
    requirements_inflow = zeros((num_nodes, num_nodes, num_comm, num_comm))
    
    # describe the flow gain/loss or transformation between commodities
    transformation = zeros((num_nodes, num_nodes, num_comm, num_comm))
    concurrent_outflow = zeros((num_nodes, num_nodes, num_comm))
    concurrent_inflow = zeros((num_nodes, num_nodes, num_comm))
    flow_bounds = fill(999, (num_nodes, num_nodes, num_comm))

    demand[4, 3] = 9999  # unlimited diesel resources
    demand[5, 1] = 9999  # unlimited water resources
    demand[1, 1] = -3.0  # household requires 3 water
    demand[1, 2] = -5.0  # household requires 5 electricity

    outflow_cost[4, 2, 3] = 0.20  # diesel costs 0.20 £/kWh
    outflow_cost[2, 1, 2] = 0.01  # electricity distribution costs 0.01 £/kWh
    outflow_cost[2, 3, 2] = 0.01  # electricity distribution costs 0.01 £/kWh
    outflow_cost[3, 1, 1] = 0.03  # water production costs 0.03 £/l
    outflow_cost[5, 3, 1] = 0.0001  # water extraction costs little
    
    transformation[2, 2, 3, 2] = 1.0  # power plant requires diesel to produce electricity
    # transformation[2, 2, 2, 2] = 1.0  # power plant produces electricity  
    transformation[3, 3, 1, 1] = 1.0  # water plant produces water
    transformation[5, 5, 1, 1] = 1.0  # water resource produces water

    # water plant
    requirements_outflow[3, 3, 2, 1] = 0.2
    requirements_inflow[3, 3, 2, 1] = 1.0
    requirements_outflow[3, 3, 1, 1] = 1.0
    requirements_inflow[3, 3, 1, 1] = 1.0
    
    # power plant
    requirements_outflow[2, 2, 3, 2] = 3.0  # power plant requires 3 kWh diesel per 1 kWh electricity
    requirements_inflow[2, 2, 3, 2] = 1.0  # power plant takes diesel as an input
    # requirements_outflow[2, 2, 2, 2] = 1.0  # power plant produces electricity as an output
    requirements_inflow[2, 2, 2, 2] = 1.0  # power plant produces electricity as an output
    
    requirements_outflow[2, 1, 2, 2] = 1.0  # power plant sends electricity to household
    requirements_outflow[2, 3, 2, 2] = 1.0  # power plant sends electricity to water plant

    requirements_outflow[4, 2, 3, 3] = 1.0  # diesel is produced from diesel resource
    requirements_inflow[4, 2, 3, 3] = 1.0  # electricity is produced from diesel
    requirements_inflow[2, 3, 2, 1] = 1.0  # water plant requires 0.2 electricity to produce water
    requirements_outflow[3, 1, 1, 1] = 1.0  # water plant produces water
    requirements_inflow[3, 1, 1, 1] = 1.0
    requirements_inflow[2, 1, 2, 2] = 1.0
    requirements_outflow[5, 3, 1, 1] = 1.0  # water is produced from the water resource
    requirements_inflow[5, 3, 1, 1] = 1.0  # water is consumed by the water plant



    model = Model(with_optimizer(GLPK.Optimizer))

    outflow = @variable(model, 
                        outflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm], 
                        lower_bound = 0, 
                        upper_bound = flow_bounds[i, j, k])
    inflow = @variable(model, 
                       inflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm], 
                       lower_bound = 0,
                       upper_bound = flow_bounds[i, j, k])
    
    @objective(
        model, 
        Min, 
        sum(outflow_cost[i, j, k] * outflow[i, j, k] 
            + inflow_cost[i, j, k] * inflow[i, j, k] 
            for i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm)
        )
        

    requirements_outflow_const = @expression(
        model,
        requirements_outflow_const[i=1:num_nodes, k=1:num_comm],
        if haskey(outflow_edges, i)
            sum(
                sum(requirements_outflow[i, j, k, l] for l in 1:num_comm) * outflow[i, j, k]
                for j in outflow_edges[i]
            )
                
        else
            println("No outflow edge found for $(i)")
            0
        end
    )

    requirements_inflow_const = @expression(
        model,
        requirements_inflow_const[i=1:num_nodes, k=1:num_comm],
        if haskey(inflow_edges, i)
            sum(
                sum(requirements_inflow[h, i, k, m] for m in 1:num_comm) * inflow[h, i, k]
                for h in inflow_edges[i]
            )

        else
            println("No inflow edge found for $(i)")
            0
        end
    )

    mass_balance = @constraint(
        model,
        mass_balance[i=1:num_nodes, k=1:num_comm],
        requirements_outflow_const[i, k] - requirements_inflow_const[i, k]
        <= demand[i, k])

    flow_transformation = @constraint(
        model,
        flow_transformation[i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm],
        sum(transformation[i, j, l, k] * outflow[i, j, l] for l in 1:num_comm) 
        == inflow[i, j, k])

    self_constraint_outflow = @constraint(
        model,
        self_constraint_outflow[i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm],
        concurrent_outflow[i, j, k] * outflow[i, j, k] <= 0)

    self_constraint_inflow = @constraint(
        model,
        self_constraint_inflow[i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm],
        concurrent_inflow[i, j, k] * inflow[i, j, k] <= 0)

    return model
end

function populate_gmncf(model::JuMP.Model)

    println(model[:mass_balance])
    println(model[:flow_transformation])

end

@time model = formulate_gmcnf(verbose = true)
@time populate_gmncf(model)

println("Compiled model, now running")
@time JuMP.optimize!(model)
println("Finished running, objective: $(JuMP.objective_value(model))")

@test JuMP.termination_status(model) == MOI.OPTIMAL
@test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
# @test JuMP.objective_value(model) == 225700.0

outflow = model[:outflow]

for var in outflow
    if JuMP.value(var) != 0.0
        println("$(var): $(JuMP.value(var))")
    end
end

@test JuMP.value(outflow[3, 1, 1]) == 3.0
@test JuMP.value(outflow[2, 1, 2]) == 5.0
@test JuMP.value(outflow[4, 2, 3]) == 15.0
@test JuMP.value(outflow[2, 3, 2]) == 3.0 * 0.2

# @time gmcnf(verbose = true)
# @time gmcnf(verbose = true)