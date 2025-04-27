window.BENCHMARK_DATA = {
  "lastUpdate": 1745753114713,
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
          "id": "819d44d9d688d9b740760f940f56683dcd54fbfe",
          "message": "Merge pull request #79 from lcauser/dependabot/github_actions/codecov/codecov-action-5.4.2\n\nBump codecov/codecov-action from 5.3.1 to 5.4.2",
          "timestamp": "2025-04-27T12:23:37+01:00",
          "tree_id": "3aed866e5b6f3cba5f026fa40f06c9770a2a55dd",
          "url": "https://github.com/lcauser/TeNe.jl/commit/819d44d9d688d9b740760f940f56683dcd54fbfe"
        },
        "date": 1745753113785,
        "tool": "julia",
        "benches": [
          {
            "name": "dmrg",
            "value": 78649287.5,
            "unit": "ns",
            "extra": "gctime=2076549\nmemory=30810384\nallocs=331505\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "state_optimiser",
            "value": 114911425,
            "unit": "ns",
            "extra": "gctime=2275932\nmemory=37764480\nallocs=529362\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}