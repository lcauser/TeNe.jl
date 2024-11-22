window.BENCHMARK_DATA = {
  "lastUpdate": 1732269118666,
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
          "id": "58a58b681712b93de187502507505063030e1b6f",
          "message": "Merge pull request #67 from lcauser/dependabot/github_actions/codecov/codecov-action-5.0.2\n\nBump codecov/codecov-action from 4.6.0 to 5.0.2",
          "timestamp": "2024-11-22T09:47:41Z",
          "tree_id": "c16bd7d753906510fea646ba79a22d49a6ec81ed",
          "url": "https://github.com/lcauser/TeNe.jl/commit/58a58b681712b93de187502507505063030e1b6f"
        },
        "date": 1732269118296,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 58964737,
            "unit": "ns",
            "extra": "gctime=2085085\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 97152133,
            "unit": "ns",
            "extra": "gctime=2640808\nmemory=36225808\nallocs=512358\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}