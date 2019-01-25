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

    cap_cost = zeros(Float64, num_nodes, num_nodes)
    flow_cost[1, 4] = 0.05
    flow_cost[4, 2] = 0.05
    flow_cost[4, 3] = 0.05

    existing_bandwidth = zeros(Float64, num_nodes, num_nodes)
    existing_bandwidth[1, 4] = 110
    existing_bandwidth[4, 2] = 10
    existing_bandwidth[4, 3] = 100

    demand = [200 -100 -100 0]

    bandwidth = [10 10 10]

    flow = @variable(
        model,
        flow[i=keys(outflow_edges), j=outflow_edges[i]],
        lower_bound = 0
    )

    new_bandwidth = @variable(
        model,
        new_bandwidth[i=keys(outflow_edges), j=outflow_edges[i]],
        lower_bound = 0
    )

    total_bandwidth = @variable(
        model,
        total_bandwidth[i=keys(outflow_edges), j=outflow_edges[i]],
        lower_bound = 0
    )

    @objective(
        model, 
        Min, 
        sum(flow_cost[i, j] * flow[i, j] + cap_cost[i, j] + new_bandwidth[i, j] for i=keys(outflow_edges), j=outflow_edges[i])
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

    total_cap = @constraint(
        model,
        total_cap[i=keys(outflow_edges), j in outflow_edges[i]],
        total_bandwidth[i, j] == existing_bandwidth[i, j] + new_bandwidth[i, j]
    )

    edge_bandwidth = @constraint(
        model,
        edge_bandwidth[i=keys(outflow_edges), j in outflow_edges[i]],
        flow[i, j] <= total_bandwidth[i, j]
    )

    for i in 1:num_nodes
        if haskey(outflow_edges, i) & haskey(inflow_edges, i)
            for j in outflow_edges[i]
                for k in inflow_edges[i]
                    @constraint(model, flow[i, j] <= total_bandwidth[k, i])
                end
            end
        end
    end

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
InfraFlow.print_vars(model[:new_bandwidth])
InfraFlow.print_vars(model[:total_bandwidth])
InfraFlow.print_constraint(model[:mass_balance])
InfraFlow.print_constraint(model[:edge_bandwidth])
end