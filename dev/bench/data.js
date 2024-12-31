window.BENCHMARK_DATA = {
  "lastUpdate": 1735638844788,
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
          "id": "eedae1228aee3f5fa61fe8c8889448d8b31abacb",
          "message": "Merge pull request #75 from lcauser/compathelper/new_version/2024-12-31-01-14-31-057-01188373983\n\nCompatHelper: add new compat entry for KernelAbstractions at version …",
          "timestamp": "2024-12-31T09:51:02Z",
          "tree_id": "5cd6bd0b43ac04508023a55e616f27c3e227df01",
          "url": "https://github.com/lcauser/TeNe.jl/commit/eedae1228aee3f5fa61fe8c8889448d8b31abacb"
        },
        "date": 1735638844455,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 59826516,
            "unit": "ns",
            "extra": "gctime=1955522\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 101499187.5,
            "unit": "ns",
            "extra": "gctime=2271248\nmemory=32403984\nallocs=459035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}