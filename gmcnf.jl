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
function gmcnf(; verbose = true)
    nodes = [1 2 3]
    edges = Dict(1=>2, 3=>2)
    commodities = ["energy", "water"]
    
    
    num_nodes = length(nodes)
    num_edges = length(edges)
    num_comm = length(commodities)
    
    # demand at node by commodity
    demand = [[0, 0]; 
              [10, 10]; 
              [0, 0]]

    inflow_cost = zeros((num_nodes, num_nodes, num_comm))
    outflow_cost = zeros((num_nodes, num_nodes, num_comm))

    model = Model(with_optimizer(GLPK.Optimizer))

    @variable(model, 
              outflow[i=keys(edges), j=edges[i], k=1:num_comm], 
              lower_bound = 0, 
              base_name="outflow")
    @variable(model, 
              inflow[i=keys(edges), j=edges[i], k=1:num_comm], 
              lower_bound = 0,
              base_name="inflow")
    
    @objective(model, Min, sum(outflow_cost[i, j, k] * outflow[i, j, k] 
                             + inflow_cost[i, j, k] * inflow[i, j, k] 
                             for i in keys(edges), j in edges[i], k in 1:num_comm))


                                
    println("Compiled model, now running")
    JuMP.optimize!(model)
    println("Finished running, objective: $(JuMP.objective_value(model))")

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    # @test JuMP.objective_value(model) == 225700.0

end

@time gmcnf(verbose = true)
@time gmcnf(verbose = true)
@time gmcnf(verbose = true)