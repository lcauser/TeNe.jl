window.BENCHMARK_DATA = {
  "lastUpdate": 1732380620221,
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
          "id": "bf7fd96af4149ad7f39aae1c69c887e09ef724f0",
          "message": "Merge pull request #70 from lcauser/improvement/lc/consolidatiate_gmp_optimiser\n\nConsolidated support for general gmps arguments",
          "timestamp": "2024-11-23T16:48:29Z",
          "tree_id": "8d62359b47ae94eac37b1d8c0a1cbd2ab2600561",
          "url": "https://github.com/lcauser/TeNe.jl/commit/bf7fd96af4149ad7f39aae1c69c887e09ef724f0"
        },
        "date": 1732380619949,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 58716929,
            "unit": "ns",
            "extra": "gctime=2063841\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 103238727,
            "unit": "ns",
            "extra": "gctime=2415245.5\nmemory=32406256\nallocs=459041\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}