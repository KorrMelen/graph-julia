import Base.copy

mutable struct tabgraph
    nv::Int
    ne::Int
    adj::Array{Bool,2}
    pds::Array{Real,2}
    pred::Array{Int,2}
    tabgraph(nbs::Int) =
        new(nbs,0,
            zeros(Bool,(nbs,nbs)),
            Array{Real}(undef,nbs,nbs),
            Array{Int}(undef,nbs,nbs))
end

function relie!(G::tabgraph,s::Int,d::Int,p::Int)
    if s <= G.nv && d <= G.nv && s!=d
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

function voisins(G::tabgraph,s::Int)
    [i for i in 1:G.nv if G.adj[s,i]]
end

function preds(G::tabgraph,s::Int)
    [i for i in 1:G.nv if G.adj[i,s]]
end

function lire!(G::tabgraph,nomf)
    f = open(nomf)
    for line in eachline(f)
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
                F.adj[j,k] = true
                #relie!(F,j,k,-1)
            end
        end
    end
    return F
end


function connexitéforte2(G::tabgraph)
    F = tabgraph(G.nv)
    F.adj = copy(G.adj)
    F.pds = copy(G.pds)
    m = (F.adj.|(F.adj*F.adj)).> 0
    while m != F.adj
        F.adj = m
        m = (F.adj.|(F.adj*F.adj)).> 0
    end
    return F
end

function comparTime(taille,proba)
        G = tabgraph(taille)
        alea!(G,proba)
        println(mesure(connexitéforte,G))
        println(mesure(connexitéforte2,G))
end

function mesure(connexite,G::tabgraph)
    connexite(G)
    t=time()
    for i in 1:1000
        connexite(G)
    end
    return (time()-t)/1000
end

function courtchemin(G::tabgraph)
    F = tabgraph(G.nv)
    F.adj = copy(G.adj)
    F.pds = copy(G.pds)
    for i in 1:F.nv
        voisin = voisins(F,i)
        pred = preds(F,i)
        for j in pred
            for k in voisin
                if F.adj[j,k]
                    F.pds[j,k] = min(F.pds[j,k],F.pds[j,i]+F.pds[i,k])
                else
                    relie!(F,j,k,F.pds[j,i]+F.pds[i,k])
                end
            end
        end
    end
    return F
end




g = tabgraph(5)
f = tabgraph(5)

lire!(g,"test")
lire!(f,"test2")
