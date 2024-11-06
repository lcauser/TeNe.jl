window.BENCHMARK_DATA = {
  "lastUpdate": 1730914878849,
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
          "id": "145ffa5f2031fc33da13fb35523eb9d9c1d26be5",
          "message": "Merge pull request #62 from lcauser/lc/bosonic_lattice_sites\n\nBosonic site types",
          "timestamp": "2024-11-06T17:36:51Z",
          "tree_id": "041ee3b09efcc5dfd0c7fd0e3aa4de09c42509c5",
          "url": "https://github.com/lcauser/TeNe.jl/commit/145ffa5f2031fc33da13fb35523eb9d9c1d26be5"
        },
        "date": 1730914878292,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 60599672,
            "unit": "ns",
            "extra": "gctime=2478344\nmemory=28133904\nallocs=301924\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 106390354,
            "unit": "ns",
            "extra": "gctime=2955954\nmemory=32401888\nallocs=459037\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}