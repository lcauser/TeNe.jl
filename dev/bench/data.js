window.BENCHMARK_DATA = {
  "lastUpdate": 1745754991784,
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
          "id": "90ae743a440dbea608bd0bfa04a9f47ba3a2c7c0",
          "message": "Merge pull request #81 from lcauser/formatter-workflow\n\nFormatting workflow",
          "timestamp": "2025-04-27T12:54:58+01:00",
          "tree_id": "47685675e0fa01f5fc2127050fc9326510cb05d0",
          "url": "https://github.com/lcauser/TeNe.jl/commit/90ae743a440dbea608bd0bfa04a9f47ba3a2c7c0"
        },
        "date": 1745754991432,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 78925576.5,
            "unit": "ns",
            "extra": "gctime=2245915.5\nmemory=30810384\nallocs=331505\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 115788652,
            "unit": "ns",
            "extra": "gctime=2557575.5\nmemory=33939344\nallocs=476042\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}