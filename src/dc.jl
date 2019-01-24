module Fixed

using JuMP, GLPK, Test, Logging
using InfraFlow
const MOI = JuMP.MathOptInterface

struct Edge
    source::AbstractString
    sink::AbstractString
    bandwidth::AbstractFloat
end

struct Node
    name::AbstractString
end

function formulate!(model)

    nodes = [1 2 3 4]
    edges = [1 4; 4 2; 4 3;]

    source_nodes = view(edges, :, 1)
    sink_nodes = view(edges, :, 2)

    outflow_edges = InfraFlow.make_edge_dict(source_nodes, sink_nodes)
    inflow_edges = InfraFlow.make_edge_dict(sink_nodes, source_nodes)

    @show outflow_edges inflow_edges

    num_nodes = length(nodes)
    num_edges = length(edges)

    flow_cost = zeros(Float64, num_nodes, num_nodes)
    flow_cost[1, 4] = 0.1
    flow_cost[4, 2] = 0.1
    flow_cost[4, 3] = 0.1

    capacity = zeros(Float64, num_nodes, num_nodes)
    capacity[1, 4] = 110
    capacity[4, 2] = 10
    capacity[4, 3] = 100

    demand = [200 -10 -100 0]

    bandwidth = [10 10 10]

    flow = @variable(
        model,
        flow[i=keys(outflow_edges), j=outflow_edges[i]],
        lower_bound = 0
    )

    new_capacity = @variable(
        model,
        new_capacity[i=keys(outflow_edges), j=outflow_edges[i]],
        lower_bound = 0
    )

    @objective(
        model, 
        Min, 
        sum(flow_cost[i, j] * flow[i, j] for i=keys(outflow_edges), j=outflow_edges[i])
    )

    @expression(
        model,
        mass_balance_expression[i=1:num_nodes],
        if haskey(outflow_edges, i) & haskey(inflow_edges, i)
            (sum(flow[i, j] for j in outflow_edges[i]) - sum(flow[k, i] for k in inflow_edges[i]))
        elseif haskey(outflow_edges, i) ~haskey(inflow_edges, i)
            (sum(flow[i, j] for j in outflow_edges[i]))   
        elseif haskey(inflow_edges, i) ~haskey(outflow_edges, i)
            (- sum(flow[j, i] for j in inflow_edges[i]))
        else
            0
        end
    )        

    mass_balance = @constraint(
        model,
        mass_balance[i=1:num_nodes],
        mass_balance_expression[i] <= demand[i]
    )

    edge_capacity = @constraint(
        model,
        edge_capacity[i=keys(outflow_edges), j in outflow_edges[i]],
        flow[i, j] <= capacity[i, j]
    )

end


model = Model(with_optimizer(GLPK.Optimizer))
formulate!(model)
println("Compiled model, now running")
@time JuMP.optimize!(model)
println("Finished running, objective: £$(JuMP.objective_value(model))")

status = JuMP.termination_status(model) == MOI.OPTIMAL

@show JuMP.termination_status(model)
@show JuMP.primal_status(model)

feasible = JuMP.primal_status(model) == MOI.FEASIBLE_POINT

InfraFlow.print_vars(model[:flow])
InfraFlow.print_constraint(model[:mass_balance])
InfraFlow.print_constraint(model[:edge_capacity])
end