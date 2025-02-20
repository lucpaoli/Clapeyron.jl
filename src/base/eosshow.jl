

function eosshow(io::IO, ::MIME"text/plain", Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    length(model) == 1 && println(io, " with 1 component:")
    length(model) > 1 && println(io, " with ", length(model), " components:")
    for i in 1:length(model)
        print(io, " \"", model.components[i], "\"")
        i != length(model) && println(io)
    end
    if hasfield(typeof(model),:params)
        println(io)
        paramnames = fieldnames(typeof(model.params))
        len_params = length(paramnames)
        !iszero(len_params) && print(io,"Contains parameters: ")
        firstloop = true
        for fieldname in paramnames
            firstloop == false && print(io, ", ")
            print(io, fieldname)
            firstloop = false
        end
    end
    show_references(io,model)
end

function eosshow(io::IO, Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    firstloop = true
    print(io, "(")
    for i in 1:length(model.components)
        firstloop == false && print(io, ", ")
        print(io, "\"", model.components[i], "\"")
        firstloop = false
    end
    print(io, ")")
end

function gc_eosshow(io::IO, ::MIME"text/plain", Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    length(model) == 1 && println(io, " with 1 component:")
    length(model) > 1 && println(io, " with ", length(model), " components:")
    param = model.groups
    for i in 1:length(param.components)
        
        print(io, " \"", param.components[i], "\": ")
        firstloop = true
        for j in 1:length(param.n_groups[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.groups[i][j], "\" => ", param.n_groups[i][j])
            firstloop = false
        end
        i != length(param.components) && println(io)
    end
    println(io)
    print(io,"Group Type: ",param.grouptype)
    fields = fieldnames(typeof(model.params))
    if length(fields) != 0
        println(io)
        print(io, "Contains parameters: ")
        firstloop = true
        for fieldname in fields
            firstloop == false && print(io, ", ")
            print(io, fieldname)
            firstloop = false
        end
    end
    show_references(io,model)
end

function gc_eosshow(io::IO, Base.@nospecialize(model::EoSModel))
    return eosshow(io,model)
end

function show_references(io::IO,model)
    if get(ENV,"CLAPEYRON_SHOW_REFERENCES","FALSE") == "TRUE"
        citations = cite(model)
        iszero(length(citations)) && return nothing #do not do anything if there isnt any citations
        println(io)
        print(io,"References: ") 
        for (i,doi) in enumerate(cite(model))
            i != 1 && print(io,", ")
            print(io,doi)   
        end
    end
    return nothing
end

export eosshow