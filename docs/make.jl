using Documenter
using DocumenterMarkdown
using TorchPID

makedocs(
    sitename="TorchPID.jl",
    modules=[TorchPID],
    format=Markdown(),
    build="buildmd",
    pages=[
        "Home" => "index.md",
        "TORCH Detector" => "torch.md",
        "Guide" => "guide.md",
        "API" => [
            "Overview" => "api/api.md",
            "Physics Constants" => "api/physics_constants.md",
            "Hit Models" => "api/hit_models.md",
            "Patter Matcher" => "api/pattern_matcher.md",
        ],
    ],
    checkdocs=:none,
    linkcheck=false,
    doctest=false,
)

makedocs(
    sitename="TorchPID.jl",
    modules=[TorchPID],
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://rrabadan.github.io/TorchPID.jl",
        #assets=String[],
        assets=["assets/custom.css"],  # List specific files
        highlights=["yaml"]
    ),
    pages=[
        "Home" => "index.md",
        "TORCH Detector" => "torch.md",
        "Guide" => "guide.md",
        "API" => [
            "Overview" => "api/api.md",
            "Physics Constants" => "api/physics_constants.md",
            "Hit Models" => "api/hit_models.md",
            "Patter Matcher" => "api/pattern_matcher.md",
        ],
    ],
    checkdocs=:none,
    linkcheck=false,
    doctest=true,
)

deploydocs(
    repo="github.com/rrabadan/TorchPID.jl.git",
    target="gh-pages",
    push_preview=true,
    devbranch="main",
)