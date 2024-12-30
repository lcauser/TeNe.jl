window.BENCHMARK_DATA = {
  "lastUpdate": 1735547247452,
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
          "id": "f968a3aa88eb7604e0065706e1a32d2dbbc881c0",
          "message": "Merge pull request #74 from lcauser/compathelper/new_version/2024-12-30-01-17-13-306-02286174077\n\nCompatHelper: add new compat entry for Documenter at version 1, (keep…",
          "timestamp": "2024-12-30T08:25:35Z",
          "tree_id": "2a8896a6bb4794d2684aa44948bc26addeb32219",
          "url": "https://github.com/lcauser/TeNe.jl/commit/f968a3aa88eb7604e0065706e1a32d2dbbc881c0"
        },
        "date": 1735547246528,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 60608834,
            "unit": "ns",
            "extra": "gctime=2114394\nmemory=28135072\nallocs=301932\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 106003022.5,
            "unit": "ns",
            "extra": "gctime=2489214\nmemory=36235632\nallocs=512370\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}