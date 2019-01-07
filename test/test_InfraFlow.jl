
module Test_InfraFlow
    using Test, InfraFlow

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