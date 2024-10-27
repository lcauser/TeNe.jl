window.BENCHMARK_DATA = {
  "lastUpdate": 1730056129665,
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
          "id": "20e70681a793364c73b55330193453709404aac9",
          "message": "Merge pull request #58 from lcauser/tests/benchmarking\n\nfixing workflows",
          "timestamp": "2024-10-27T19:04:28Z",
          "tree_id": "62dfc5dc7cccf95a8269f9ed6ac4cb6a8b2b1a47",
          "url": "https://github.com/lcauser/TeNe.jl/commit/20e70681a793364c73b55330193453709404aac9"
        },
        "date": 1730056128744,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 60437972,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28133904\nallocs=301924\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 60787655,
            "unit": "ns",
            "extra": "gctime=2204388\nmemory=28133904\nallocs=301924\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}