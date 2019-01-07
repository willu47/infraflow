using Documenter
using InfraFlow

makedocs(
    sitename = "InfraFlow",
    format = Documenter.HTML(),
    modules = [InfraFlow]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
