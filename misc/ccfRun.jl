# ccfRun.jl -- called by ccfJl in Buffering.wl.
# julia ccfRun.jl Rpat.tsv C.tsv out.m
# Rpat: reactions x species 0/1 (rate of reaction rho depends on species m);
# C: reactions x d integer basis of ker S (d may be 0, file then empty).
# Writes a Mathematica list {blocks, orderPairs, covers, C0, Cinf, R0T, RinfT}
# with blocks = {{reactions}, {species}, {cycleRows}} in CCF order.
using LayeredMixedMatrices, DelimitedFiles

R = Bool.(readdlm(ARGS[1], '\t', Int))
nr = size(R, 1)
Craw = filesize(ARGS[2]) == 0 ? zeros(Int, nr, 0) : readdlm(ARGS[2], '\t', Int)
C = reshape(Craw, nr, :)
res = ccf(R, C; columnLM = true)

ml(v) = "{" * join(string.(collect(v)), ",") * "}"
mlp(ps) = "{" * join(["{" * string(p[1]) * "," * string(p[2]) * "}" for p in ps], ",") * "}"
blocks = "{" * join([ "{" * ml(res.Cblocks[i]) * "," * ml(res.RblocksT[i]) * "," * ml(res.RblocksQ[i]) * "}"
                      for i in eachindex(res.Cblocks)], ",") * "}"
open(ARGS[3], "w") do io
    print(io, "{", blocks, ",", mlp(column_block_order_pairs(res)), ",", mlp(block_poset(res)), ",",
          ml(res.C0), ",", ml(res.C∞), ",", ml(res.R0T), ",", ml(res.R∞T), "}")
end
