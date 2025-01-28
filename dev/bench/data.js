window.BENCHMARK_DATA = {
  "lastUpdate": 1738061447509,
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
          "id": "7b750299e3c30ac2ea573e143d001f0d23802ef6",
          "message": "Merge pull request #77 from lcauser/dependabot/github_actions/codecov/codecov-action-5.3.1\n\nBump codecov/codecov-action from 5.1.2 to 5.3.1",
          "timestamp": "2025-01-28T10:46:21Z",
          "tree_id": "37ce16b0434afa48c3040a0dc9eaf5232330e723",
          "url": "https://github.com/lcauser/TeNe.jl/commit/7b750299e3c30ac2ea573e143d001f0d23802ef6"
        },
        "date": 1738061447159,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 80089882,
            "unit": "ns",
            "extra": "gctime=2380565\nmemory=30808272\nallocs=331505\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 111566378,
            "unit": "ns",
            "extra": "gctime=2629735\nmemory=37764752\nallocs=529359\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}