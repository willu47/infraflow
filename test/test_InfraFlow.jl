
module Test_InfraFlow

using Test, InfraFlow, YAML

@testset "test_make_edges_dict" begin

    nodes = ["household", "power_station", "diesel_resource"]
    edges = [3 2; 2 1; 2 2; 3 3;]
    commodities = ["electricity", "diesel"]
    discount_rate = 0.10

    years = [2018, 2019, 2020, 2025]
        
    source_nodes = view(edges, :, 1)
    sink_nodes = view(edges, :, 2)

    outflow_edges = InfraFlow.make_edge_dict(source_nodes, sink_nodes)
    @test outflow_edges == Dict{Int8, Array{Int8}}(3=>[2, 3], 2=>[1, 2])
    inflow_edges = InfraFlow.make_edge_dict(sink_nodes, source_nodes)
    @test inflow_edges == Dict{Int8, Array{Int8}}(1=>[2], 2=>[3, 2], 3=>[3])

end

@testset "load_data" begin

    data = InfraFlow.get_data("./test_data.yml")

    @test data["nodes"] == ["household", "power_station", "diesel_resource"]
    @test data["edges"] == [3 2; 2 1; 2 2;]

    @test data["commodities"] == ["electricity", "diesel"]
    @test data["discount_rate"] == 0.10

    @test data["years"] == [2018, 2019, 2020, 2025]
    
    num_nodes = 3
    num_comm = 2
    num_years = 4

    # demand at node by commodity
    demand = zeros((num_nodes, num_comm, num_years))
    
    # Demands and Resources
    demand[3, 2, :] .= 999999.0  # unlimited diesel resources
    demand[1, 1, 1] = -5000.0  # household requires 5 MWh electricity
    demand[1, 1, 2] = -6000.0  # household requires 6 MWh electricity
    demand[1, 1, 3] = -7000.0  # household requires 7 MWh electricity
    demand[1, 1, 4] = -8000.0  # household requires 8 MWh electricity

    @test data["demand"] == demand

    capacity_cost = zeros((num_nodes, num_nodes, num_comm))

    # Investment costs
    capacity_cost[2, 2, 1] = 800  # £/kW for the diesel plant
    capacity_cost[3, 3, 2] = 10  # £/kW for diesel imports

    @test data["capacity_cost"] == capacity_cost

    cap2act = zeros((num_nodes, num_nodes, num_comm))

    # Relate new_capacity of a node to its activity
    cap2act[2, 1, 1] = 1
    cap2act[2, 2, 2] = 1
    cap2act[2, 2, 1] = 8760  # Power plant produce 8760 kWh electricity per year per kW new_capacity
    cap2act[3, 2, 2] = 1
    cap2act[3, 3, 2] = 8760  # Diesel resource produce 8760 kWh diesel per year per kW new_capacity

    @test data["cap2act"] == cap2act

    inflow_cost = zeros((num_nodes, num_nodes, num_comm))
    outflow_cost = zeros((num_nodes, num_nodes, num_comm))

    # Operational costs
    outflow_cost[3, 2, 2] = 0.20  # diesel costs £0.20/kWh
    outflow_cost[2, 1, 1] = 0.01  # electricity distribution costs £0.01/kWh

    @test data["inflow_cost"] == inflow_cost
    @test data["outflow_cost"] == outflow_cost
    
    # upper bound on operational decision variables (new_capacity)
    flow_bounds = zeros(Float64, (num_nodes, num_nodes, num_comm, num_years))

    flow_bounds[2, 2, 1, 1:3] .= 0.7
    flow_bounds[2, 2, 1, 4] = 0.0
    
    @test data["flow_bounds"] == flow_bounds

    # describe the flow gain/loss or transformation between commodities
    transformation = zeros(Float64, (num_nodes, num_nodes, num_comm, num_comm))
    
    # Transformation of commodities
    transformation[2, 2, 2, 1] = 1.0  # power plant requires diesel to produce electricity
    transformation[2, 1, 1, 1] = 0.93  # 7% losses in distribution of electricity
    transformation[3, 2, 2, 2] = 1.0  # no losses in distribution of diesel

    @test data["transformation"] ≈ transformation

    # describe the commodity requirements for an in- or out-flow
    requirements_outflow = zeros(Float64, (num_nodes, num_nodes, num_comm, num_comm))
    requirements_inflow = zeros(Float64, (num_nodes, num_nodes, num_comm, num_comm))
    
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

    @test data["requirements_inflow"] ≈ requirements_inflow
    @test data["requirements_outflow"] ≈ requirements_outflow


end

end
