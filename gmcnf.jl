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
    nodes = ["household", "power_station", "diesel_resource"]
    edges = [3 2; 2 1; 2 2; 3 3;]
    commodities = ["electricity", "diesel"]
     
    source_nodes = view(edges, :, 1)
    sink_nodes = view(edges, :, 2)

    outflow_edges = make_edge_dict(source_nodes, sink_nodes)
    @test outflow_edges == Dict{Int8, Array{Int8}}(3=>[2, 3], 2=>[1, 2])
    inflow_edges = make_edge_dict(sink_nodes, source_nodes)
    @test inflow_edges == Dict{Int8, Array{Int8}}(1=>[2], 2=>[3, 2], 3=>[3])

    num_nodes = length(nodes)
    num_edges = length(edges)
    num_comm = length(commodities)
    
    # demand at node by commodity
    demand = zeros((num_nodes, num_comm))

    capacity_cost = zeros((num_nodes, num_nodes, num_comm))

    inflow_cost = zeros((num_nodes, num_nodes, num_comm))
    outflow_cost = zeros((num_nodes, num_nodes, num_comm))
    
    # describe the commodity requirements for an in- or out-flow
    requirements_outflow = zeros((num_nodes, num_nodes, num_comm, num_comm))
    requirements_inflow = zeros((num_nodes, num_nodes, num_comm, num_comm))
    
    # describe the flow gain/loss or transformation between commodities
    transformation = zeros((num_nodes, num_nodes, num_comm, num_comm))
    
    # upper bound on operational decision variables (capacity)
    flow_bounds = fill(0, (num_nodes, num_nodes, num_comm))

    # Demands and Resources
    demand[3, 2] = 999999  # unlimited diesel resources
    demand[1, 1] = -5000.0  # household requires 5 MWh electricity

    # Operational costs
    outflow_cost[3, 2, 2] = 0.20  # diesel costs £0.20/kWh
    outflow_cost[2, 1, 1] = 0.01  # electricity distribution costs £0.01/kWh

    # Investment costs
    capacity_cost[2, 2, 1] = 800  # £/kW for the diesel plant
    capacity_cost[3, 3, 2] = 10  # £/kW for diesel imports
    
    # Transformation of commodities
    transformation[2, 2, 2, 1] = 1.0  # power plant requires diesel to produce electricity
    transformation[2, 1, 1, 1] = 0.93  # 7% losses in distribution of electricity
    transformation[3, 2, 2, 2] = 1.0  # no losses in distribution of diesel

    # Power plant
    requirements_outflow[2, 2, 2, 1] = 3.0  # power plant requires 3 kWh diesel per 1 kWh electricity
    requirements_outflow[2, 2, 1, 1] = 1.0  # power plant produces electricity as an output
    requirements_inflow[2, 2, 1, 1] = 1.0  # power plant produces electricity as an output
    
    # Distribution of electricity
    requirements_outflow[2, 1, 1, 1] = 1.0  # Send electricity to household from power plant
    requirements_inflow[2, 1, 1, 1] = 1.0

    # Distribution of diesel
    requirements_outflow[3, 2, 2, 2] = 1.0  # Send diesel to power plant from diesel resource
    requirements_inflow[3, 2, 2, 2] = 1.0

    model = Model(with_optimizer(GLPK.Optimizer))

    capacity = @variable(model,
                         capacity[i=1:num_nodes, j=1:num_nodes, k=1:num_comm],
                         lower_bound = 0)

    outflow = @variable(model, 
                        outflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm], 
                        lower_bound = 0)
    inflow = @variable(model, 
                       inflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm], 
                       lower_bound = 0)
    
    @constraint(
        model,
        capacity_exp_outflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm],
        outflow[i, j, k] <= flow_bounds[i, j, k] + capacity[i, j, k])

    @constraint(
        model,
        capacity_exp_inflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm],
        inflow[i, j, k] <= flow_bounds[i, j, k] + capacity[i, j, k])

    @objective(
        model, 
        Min, 
        sum(outflow_cost[i, j, k] * outflow[i, j, k] 
            + inflow_cost[i, j, k] * inflow[i, j, k] 
            + capacity_cost[i, j, k] * capacity[i, j, k]
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

    return model
end

function populate_gmncf(model::JuMP.Model)

    println(model[:mass_balance])
    println(model[:flow_transformation])

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

@time model = formulate_gmcnf(verbose = true)
@time populate_gmncf(model)

println("Compiled model, now running")
@time JuMP.optimize!(model)
println("Finished running, objective: £$(JuMP.objective_value(model))")

@test JuMP.termination_status(model) == MOI.OPTIMAL
@test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
# @test JuMP.objective_value(model) == 225700.0

outflow = model[:outflow]
inflow = model[:inflow]
capacity = model[:capacity]

cap_con = model[:capacity_exp_outflow]

print_vars(outflow)
print_vars(inflow)
print_vars(capacity)
print_duals(cap_con)

@test JuMP.value(inflow[2, 1, 1]) ≈ 5000
@test JuMP.value(outflow[2, 1, 1]) ≈ 5376.344086021505
@test JuMP.value(outflow[2, 2, 2]) ≈ 5376.344086021505
@test JuMP.value(outflow[3, 2, 2]) ≈ 16129.032258064515

@time JuMP.optimize!(model)
@time JuMP.optimize!(model)
