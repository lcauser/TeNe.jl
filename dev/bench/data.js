window.BENCHMARK_DATA = {
  "lastUpdate": 1730060629046,
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
          "id": "2af775062c3882f6f682c8dd167cd3d3d385b7e9",
          "message": "Merge pull request #60 from lcauser/fix/tensorproduct\n\nImproved test coverage on tensorproduct",
          "timestamp": "2024-10-27T20:22:14Z",
          "tree_id": "8f1a4651862972ba4e6cd5b2e8b6197b97c754c2",
          "url": "https://github.com/lcauser/TeNe.jl/commit/2af775062c3882f6f682c8dd167cd3d3d385b7e9"
        },
        "date": 1730060628710,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 60139095.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28133904\nallocs=301924\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 104804244,
            "unit": "ns",
            "extra": "gctime=2358523\nmemory=32403312\nallocs=459035\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}