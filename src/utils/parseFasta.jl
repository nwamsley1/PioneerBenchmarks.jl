function parseFasta(fasta_path::String)

    function getReader(fasta_path::String)
        if endswith(fasta_path, ".fasta.gz")
            return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
        elseif endswith(fasta_path, ".fasta")
            return FASTA.Reader(open(fasta_path))
        else
            throw(ErrorException("fasta_path \"$fasta_path\" did not end with `.fasta` or `.fasta.gz`"))
        end
    end

    function parse_identifier(x)
        return split(x, "|")
        if length(x_split) >=2
            return x_split[2]
        else
            return x_split[1]
        end
    end
    #I/O for Fasta
    reader = getReader(fasta_path)
    fasta = Vector{Tuple{String, String}}()
    #@time begin
        for record in reader
                try 
                    protein_name = split(FASTA.identifier(record), "|")
                    #a = split(parse_identifier(FASTA.identifier(record)), "|")
                    if length(protein_name)>=2
                        protein_name = protein_name[2]
                        species_name = split(split(split(FASTA.description(record), ' ')[1], '|')[end], '_')[end]
                        if species_name == "MYCHR"
                            species_name = "CONTAMINANT"
                        end
                    else
                        protein_name = protein_name[1]
                        species_name = "CONTAMINANT"
                    end

                        push!(fasta, 
                                (
                                    protein_name, 
                                    species_name
                                    #FASTA.sequence(record),
                                    #false
                                )
                        )
                catch e
                        println("Error parsing fasta record \n : $record \n ")
                        rethrow(e)
                end
        end
    #end

    return fasta
end

function getAccessionToSpecies(
    fasta_dir::String
)
    acc_spec_pairs = Vector{
        Vector{Tuple{String, String}}
    }()

    fasta_paths = readdir(fasta_dir)
    n_fasta_paths = 0
    for fasta_path in fasta_paths
        if endswith(fasta_path, "fasta.gz") | endswith(fasta_path,".fasta")
            push!(
                acc_spec_pairs,
                parseFasta(joinpath(fasta_dir, fasta_path))
            )
            n_fasta_paths += 1
        end
    end
    if n_fasta_paths == 0
        throw("The fasta directory supplied did not contain any files ending in '.fasta.gz' or '.fasta'!. \n Directory was \n $fasta_dir")
    end
    return Dict(vcat(acc_spec_pairs...))
end