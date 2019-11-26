import Base.copy

mutable struct tabgraph
    nv::Int
    ne::Int
    adj::Array{Bool,2}
    pds::Array{Real,2}
    pred::Array{Int,2}
    tabgraph(nbs::Int32) =
        new(nbs,0,
            zeros(Bool,(nbs,nbs)),
            Array{Real}(undef,nbs,nbs),
            Array{Int}(undef,nbs,nbs))
end

function relie!(G::tabgraph,s::Int32,d::Int32,p::Float64)
    if s <= G.nv && d <= G.nv
        G.adj[s,d] = true
        G.pds[s,d] = p
        G.ne+=1
    end
end

function aff(G::tabgraph)
    G.pds
end

function alea!(G::tabgraph, p::Float64)
    for i in 1:G.nv
        for j in vcat(1:i-1,i+1:G.nv)
            if rand(Float64) > p
                relie!(G,i,j,rand(1:10))
            end
        end
    end
end

function voisins(G::tabgraph,s::Int32)
    [i for i in 1:G.nv if G.adj[s,i]]
end

function preds(G::tabgraph,s::Int32)
    [i for i in 1:G.nv if G.adj[i,s]]
end

function lire!(G::tabgraph,nomf)
    f = open(nomf)
    for (i, line) in enumerate(eachline(f))
        line = split(line,",")
        relie!(G,parse(Int,line[1]),parse(Int,line[2]),parse(Int,line[3]))
    end
    close(f)
end

function connexitéforte(G::tabgraph)
    F = tabgraph(G.nv)
    F.adj = copy(G.adj)
    F.pds = copy(G.pds)
    for i in 1:F.nv
        voisin = voisins(F,i)
        pred = preds(F,i)
        for j in pred
             for k in voisin
                  relie!(F,j,k,1)
              end
          end
      end
      return F
end
