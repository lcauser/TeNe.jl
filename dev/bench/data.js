window.BENCHMARK_DATA = {
  "lastUpdate": 1735393950933,
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
          "id": "6ebaa16fbd779c9c59b23173b8203fc25cc048ac",
          "message": "Merge pull request #73 from lcauser/dependabot/github_actions/codecov/codecov-action-5.1.2\n\nBump codecov/codecov-action from 5.0.2 to 5.1.2",
          "timestamp": "2024-12-28T13:48:09Z",
          "tree_id": "3ad629282edfe7eafa3172842a23fa2262d0102f",
          "url": "https://github.com/lcauser/TeNe.jl/commit/6ebaa16fbd779c9c59b23173b8203fc25cc048ac"
        },
        "date": 1735393950532,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 61124404.5,
            "unit": "ns",
            "extra": "gctime=2243704\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 111174133,
            "unit": "ns",
            "extra": "gctime=2722377\nmemory=36224992\nallocs=512368\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}