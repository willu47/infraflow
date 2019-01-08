# infraflow: Julia/JuMP implementation of a generalised multi-commodity network flow

A linear programming new_capacity expansion generalized multi-commodity network
flow implemented using JuMP/Julia based on generalized multi-commodity network
flow as found in:

    Ishimatsu, Takuto. “Generalized Multi-Commodity Network Flows : Case Studies
    in Space Logistics and Complex Infrastructure Systems.”
    Masssachusett Institute of Technology, 2013
    http://hdl.handle.net/1721.1/82470

## Installing

Follow the [instructions](https://julialang.org/downloads/) for your platform
to install Julia.

You will also need to install the packages
[JuMP](http://www.juliaopt.org/JuMP.jl/dev/installation/) and
[GLPK](https://www.gnu.org/software/glpk/)

## Running for the first time

You can then run the model from the command line using

    julia InfraFlow.jl

## Building the documentation

    julia --project=. docs/make.jl