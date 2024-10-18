function parseFasta(fasta_path::String, parse_identifier::Function = x -> split(x,"|")[2])

    function getReader(fasta_path::String)
        if endswith(fasta_path, ".fasta.gz")
            return FASTA.Reader(GzipDecompressorStream(open(fasta_path)))
        elseif endswith(fasta_path, ".fasta")
            return FASTA.Reader(open(fasta_path))
        else
            throw(ErrorException("fasta_path \"$fasta_path\" did not end with `.fasta` or `.fasta.gz`"))
        end
    end

    #I/O for Fasta
    reader = getReader(fasta_path)

    #In memory representation of FastaFile
    #fasta = Vector{FastaEntry}()
    fasta = Vector{Tuple{String, String}}()
    @time begin
        for record in reader
                push!(fasta, 
                        (parse_identifier(FASTA.identifier(record)),
                         split(split(split(FASTA.description(record), ' ')[1], '|')[end], '_')[end],
                                #FASTA.sequence(record),
                                #false
                        )
                )
        end
    end

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