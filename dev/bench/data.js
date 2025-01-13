window.BENCHMARK_DATA = {
  "lastUpdate": 1736803020775,
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
          "id": "94b608be023e218c590f6b60118e8e524114711c",
          "message": "Merge pull request #76 from lcauser/compathelper/new_version/2025-01-13-01-19-09-353-00633633208\n\nCompatHelper: add new compat entry for Revise at version 3, (keep exi…",
          "timestamp": "2025-01-13T21:12:27Z",
          "tree_id": "197a640e5e3a4030e174fa95e05ccf2b118a4a01",
          "url": "https://github.com/lcauser/TeNe.jl/commit/94b608be023e218c590f6b60118e8e524114711c"
        },
        "date": 1736803020218,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 60841947,
            "unit": "ns",
            "extra": "gctime=2263779\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 106596497,
            "unit": "ns",
            "extra": "gctime=2669994\nmemory=36229680\nallocs=512342\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}