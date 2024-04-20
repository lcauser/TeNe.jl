var documenterSearchIndex = {"docs":
[{"location":"#TeNe.jl-Documentation","page":"TeNe.jl Documentation","title":"TeNe.jl Documentation","text":"","category":"section"},{"location":"","page":"TeNe.jl Documentation","title":"TeNe.jl Documentation","text":"contract(x, y, cix, ciy)","category":"page"},{"location":"#TeNe.contract-NTuple{4, Any}","page":"TeNe.jl Documentation","title":"TeNe.contract","text":"contract(x, y, cix, ciy, [conjx=false, conjy=false]; kwargs...)\n\nContract tensors x and y across dimensions cix and ciy, and returns it as z.\n\nArguments\n\n- `x`: first tensor to contract.\n- `y': second tensor to contract.\n- `cix`: the dimensions of the first tensor to contract.\n- `ciy`: the dimensions of the second tensor to contract.\n- `conjx::Bool=false`: Take the complex conjugate of argument x?\n- `conjy::Bool=false`: Take the complex conjugate of argument y?\n\nOptional Keyword Arguments\n\n- `tocache::Bool=true`: store the result in the second level of the cache?\n- `sublevel=:auto`: if stored in cache, at which sublevel? :auto finds non-aliased memory\n\nExamples\n\njulia> x = randn(ComplexF64, 2, 3, 4);\njulia> y = randn(ComplexF64, 3, 5, 6);\njulia> z = contract(x, y, 2, 1);\njulia> size(z)\n(2, 4, 5, 6)\n\njulia> x = randn(ComplexF64, 2, 3, 4, 5);\njulia> y = randn(ComplexF64, 6, 5, 2, 7);\njulia> z = contract(x, y, (1, 4), (3, 2));\njulia> size(z)\n(3, 4, 6, 7)\n\n\n\n\n\n","category":"method"}]
}
