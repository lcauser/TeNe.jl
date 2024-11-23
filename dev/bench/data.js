window.BENCHMARK_DATA = {
  "lastUpdate": 1732375596175,
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
          "id": "d8321245420651f6f0c5a6bdc001cf5fe8b494ed",
          "message": "Merge pull request #69 from lcauser/improvement/lc/consolidatiate_gmp_optimiser\n\nMPSOptimiser allows GMPS, MPO and MPS",
          "timestamp": "2024-11-23T15:23:54Z",
          "tree_id": "0c942cf57e7acd1010bc219f5c12c17e4c8b9751",
          "url": "https://github.com/lcauser/TeNe.jl/commit/d8321245420651f6f0c5a6bdc001cf5fe8b494ed"
        },
        "date": 1732375595889,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 59386148,
            "unit": "ns",
            "extra": "gctime=1958474.5\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 98370271,
            "unit": "ns",
            "extra": "gctime=2333026\nmemory=36219376\nallocs=512366\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}