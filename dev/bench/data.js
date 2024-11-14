window.BENCHMARK_DATA = {
  "lastUpdate": 1731604847774,
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
          "id": "55551b60732c19b555aac33fb2f49828ea2ca642",
          "message": "Merge pull request #66 from lcauser/fix/lc/oplists\n\nSmall changes to OpLists",
          "timestamp": "2024-11-14T17:17:30Z",
          "tree_id": "4b4e3f8c7a5eb0aa2d3a7253e7713de11cc1b624",
          "url": "https://github.com/lcauser/TeNe.jl/commit/55551b60732c19b555aac33fb2f49828ea2ca642"
        },
        "date": 1731604846591,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 59934078,
            "unit": "ns",
            "extra": "gctime=2349971\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 107342693,
            "unit": "ns",
            "extra": "gctime=2793558\nmemory=32398128\nallocs=459037\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}