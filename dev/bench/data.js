window.BENCHMARK_DATA = {
  "lastUpdate": 1745752754728,
  "repoUrl": "https://github.com/lcauser/TeNe.jl",
  "entries": {
    "Julia benchmark result": [
      {
        "commit": {
          "author": {
            "email": "luke.causer@outlook.com",
            "name": "Luke Causer",
            "username": "lcauser"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "5f37317089f13b2f8a0e1fcce5091294e7f46a23",
          "message": "Merge pull request #80 from lcauser/formatter\n\nFormatting using JuliaFormatter.jl",
          "timestamp": "2025-04-27T12:14:30+01:00",
          "tree_id": "c85f43fd295ba9bdc2105264e3bca59ca4741589",
          "url": "https://github.com/lcauser/TeNe.jl/commit/5f37317089f13b2f8a0e1fcce5091294e7f46a23"
        },
        "date": 1745752754250,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 80477600,
            "unit": "ns",
            "extra": "gctime=2186590\nmemory=30810384\nallocs=331505\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 118786096.5,
            "unit": "ns",
            "extra": "gctime=2618947\nmemory=33964176\nallocs=476035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}