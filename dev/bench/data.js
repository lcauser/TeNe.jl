window.BENCHMARK_DATA = {
  "lastUpdate": 1730058571912,
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
          "id": "13ea5d8f03f76287a525f9b8062cf921b0af34ec",
          "message": "Merge pull request #59 from lcauser/temp/test\n\nfixed CI performance workflow",
          "timestamp": "2024-10-27T19:37:36Z",
          "tree_id": "06731164040524c19fb5c57cf7da3b6d538570a0",
          "url": "https://github.com/lcauser/TeNe.jl/commit/13ea5d8f03f76287a525f9b8062cf921b0af34ec"
        },
        "date": 1730057950155,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 60593149,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28133904\nallocs=301924\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 109223073,
            "unit": "ns",
            "extra": "gctime=2218303\nmemory=36234352\nallocs=512358\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
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
          "id": "13ea5d8f03f76287a525f9b8062cf921b0af34ec",
          "message": "Merge pull request #59 from lcauser/temp/test\n\nfixed CI performance workflow",
          "timestamp": "2024-10-27T19:37:36Z",
          "tree_id": "06731164040524c19fb5c57cf7da3b6d538570a0",
          "url": "https://github.com/lcauser/TeNe.jl/commit/13ea5d8f03f76287a525f9b8062cf921b0af34ec"
        },
        "date": 1730058571035,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 59571068,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28133904\nallocs=301924\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 106011006.5,
            "unit": "ns",
            "extra": "gctime=2416497\nmemory=36231712\nallocs=512354\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}