module InfraFlow

using JuMP, GLPK, Test, YAML
const MOI = JuMP.MathOptInterface

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

function get_data(file_path::AbstractString)
    data = YAML.load(open(file_path))

    model_data = Dict()

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
    model_data["years"] = data["years"]
    model_data["commodities"] = data["commodities"]
    model_data["discount_rate"] = data["discount_rate"]

    edges = []
    for edge in data["edges"]
        source = nodes[edge["source"]]
        sink = nodes[edge["sink"]]
        push!(edges, source, sink)
    end

    for loop in self_loops
        push!(edges, loop, loop)
    end

    edges = transpose(reshape(edges, 2, :))
    model_data["edges"] = edges

    return model_data
end


"""
    make_edge_dict(edge_nodes, other_nodes)

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

"""
    formulate_gmcnf()

``\\sum_{ijky} c_{ijky}x_{ijky}`` 

"""
function formulate_gmcnf(; verbose = true)
    nodes = ["household", "power_station", "diesel_resource"]
    edges = [3 2; 2 1; 2 2;]
    commodities = ["electricity", "diesel"]
    discount_rate = 0.10

    years = [2018, 2019, 2020, 2025]
     
    source_nodes = view(edges, :, 1)
    sink_nodes = view(edges, :, 2)

    outflow_edges = make_edge_dict(source_nodes, sink_nodes)
    inflow_edges = make_edge_dict(sink_nodes, source_nodes)

    num_nodes = length(nodes)
    num_edges = length(edges)
    num_comm = length(commodities)
    num_years = length(years)
    
    # demand at node by commodity
    demand = zeros((num_nodes, num_comm, num_years))

    capacity_cost = zeros((num_nodes, num_nodes, num_comm))
    cap2act = zeros((num_nodes, num_nodes, num_comm))

    inflow_cost = zeros((num_nodes, num_nodes, num_comm))
    outflow_cost = zeros((num_nodes, num_nodes, num_comm))
    
    # describe the commodity requirements for an in- or out-flow
    requirements_outflow = zeros((num_nodes, num_nodes, num_comm, num_comm))
    requirements_inflow = zeros((num_nodes, num_nodes, num_comm, num_comm))
    
    # describe the flow gain/loss or transformation between commodities
    transformation = zeros((num_nodes, num_nodes, num_comm, num_comm))
    
    # upper bound on operational decision variables (new_capacity)
    flow_bounds = fill(0, (num_nodes, num_nodes, num_comm))

    # Demands and Resources
    demand[3, 2, :] .= 999999.0  # unlimited diesel resources
    demand[1, 1, 1] = -5000.0  # household requires 5 MWh electricity
    demand[1, 1, 2] = -6000.0  # household requires 6 MWh electricity
    demand[1, 1, 3] = -7000.0  # household requires 7 MWh electricity
    demand[1, 1, 4] = -8000.0  # household requires 8 MWh electricity

    # Operational costs
    outflow_cost[3, 2, 2] = 0.20  # diesel costs £0.20/kWh
    outflow_cost[2, 1, 1] = 0.01  # electricity distribution costs £0.01/kWh

    # Investment costs
    capacity_cost[2, 2, 1] = 800  # £/kW for the diesel plant
    capacity_cost[3, 3, 2] = 10  # £/kW for diesel imports

    # Relate new_capacity of a node to its activity
    cap2act[2, 1, 1] = 1
    cap2act[2, 2, 2] = 1
    cap2act[2, 2, 1] = 1
    cap2act[2, 2, 1] = 8760  # Power plant produce 8760 kWh electricity per year per kW new_capacity
    cap2act[3, 2, 2] = 1
    cap2act[3, 3, 2] = 8760  # Diesel resource produce 8760 kWh diesel per year per kW new_capacity
    
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

    new_capacity = @variable(model,
                         new_capacity[i=1:num_nodes, j=1:num_nodes, k=1:num_comm, y=1:num_years],
                         lower_bound = 0)

    outflow = @variable(model, 
                        outflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm, y=1:num_years], 
                        lower_bound = 0)
    inflow = @variable(model, 
                       inflow[i=1:num_nodes, j=1:num_nodes, k=1:num_comm, y=1:num_years], 
                       lower_bound = 0)
    
    @variable(
        model,
        total_annual_capacity[i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        lower_bound = 0
    )

    @constraint(
        model,
        accumulate_capacity[i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        total_annual_capacity[i, j, k, y] == sum(new_capacity[i, j, k, z] for z in 1:y)
    )

    @constraint(
        model,
        capacity_exp_outflow[i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        outflow[i, j, k, y] <= flow_bounds[i, j, k] + total_annual_capacity[i, j, k, y] * cap2act[i, j, k])

    @constraint(
        model,
        capacity_exp_inflow[i in keys(outflow_edges), j in outflow_edges[i], k=1:num_comm, y=1:num_years],
        inflow[i, j, k, y] <= flow_bounds[i, j, k] + total_annual_capacity[i, j, k, y] * cap2act[i, j, k])
        
    function discount_factor(year) 
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
            println("No outflow edge found for $(i)")
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
            println("No inflow edge found for $(i)")
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
        flow_transformation[i in keys(outflow_edges), j in outflow_edges[i], k in 1:num_comm, y=1:num_years],
        sum(transformation[i, j, l, k] * outflow[i, j, l, y] for l in 1:num_comm) 
        == inflow[i, j, k, y])

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

function run()

    @time model = formulate_gmcnf(verbose = true)

    # println(model[:mass_balance])
    # println(model[:flow_transformation])

    println("Compiled model, now running")
    @time JuMP.optimize!(model)
    println("Finished running, objective: £$(JuMP.objective_value(model))")

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
    # @test JuMP.objective_value(model) == 225700.0

    outflow = model[:outflow]
    inflow = model[:inflow]
    new_capacity = model[:total_annual_capacity]

    print_vars(outflow)
    print_vars(inflow)
    print_vars(new_capacity)

    @test JuMP.value(inflow[2, 1, 1, 1]) ≈ 5000
    @test JuMP.value(outflow[2, 1, 1, 1]) ≈ 5376.344086021505
    @test JuMP.value(outflow[2, 2, 2, 1]) ≈ 5376.344086021505
    @test JuMP.value(outflow[3, 2, 2, 1]) ≈ 16129.032258064515

end

# run()

end
