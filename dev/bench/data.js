window.BENCHMARK_DATA = {
  "lastUpdate": 1731101986607,
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
          "id": "374c869ebdf80b3d2afb2201a7fd62ccda94e178",
          "message": "Merge pull request #64 from lcauser/lc/lioville_wrapper\n\nLiouville Wrappers",
          "timestamp": "2024-11-07T21:24:16Z",
          "tree_id": "fc28494fd453e89c1bef8a244fde46eade77b121",
          "url": "https://github.com/lcauser/TeNe.jl/commit/374c869ebdf80b3d2afb2201a7fd62ccda94e178"
        },
        "date": 1731014831549,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 61263953.5,
            "unit": "ns",
            "extra": "gctime=2431332.5\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 106150443,
            "unit": "ns",
            "extra": "gctime=2940829\nmemory=32399504\nallocs=459037\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
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
          "id": "4f7b4e1ba92aeaabd680d905788c5b9e597c675b",
          "message": "Merge pull request #65 from lcauser/lc/fix_wrapper\n\nfixed wrapper",
          "timestamp": "2024-11-08T21:38:00Z",
          "tree_id": "b46048f7db8de275ecccd53455823ebba1f320de",
          "url": "https://github.com/lcauser/TeNe.jl/commit/4f7b4e1ba92aeaabd680d905788c5b9e597c675b"
        },
        "date": 1731101985540,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 59873818.5,
            "unit": "ns",
            "extra": "gctime=1906319\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 100185454,
            "unit": "ns",
            "extra": "gctime=2514654\nmemory=36222112\nallocs=512368\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}