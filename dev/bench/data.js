window.BENCHMARK_DATA = {
  "lastUpdate": 1730063163703,
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
          "id": "3e2e5e8d9116eb03a3a64cec520a40da65864607",
          "message": "Merge pull request #61 from lcauser/test/workflow_test\n\nChanged doc make to copy benchmarks",
          "timestamp": "2024-10-27T20:50:29Z",
          "tree_id": "e583d2a8f75cc888179c767e2f1e0ae20d961fa3",
          "url": "https://github.com/lcauser/TeNe.jl/commit/3e2e5e8d9116eb03a3a64cec520a40da65864607"
        },
        "date": 1730063162785,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 60956317.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=28133904\nallocs=301924\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 106102502,
            "unit": "ns",
            "extra": "gctime=2137364.5\nmemory=36214592\nallocs=512358\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}