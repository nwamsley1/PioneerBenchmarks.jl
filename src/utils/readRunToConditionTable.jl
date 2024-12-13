function readRunToConditionTable(filepath::String)
    # Initialize empty dictionary
    run_condition_dict = Dict{String, String}()
    
    # Read the file line by line
    open(filepath) do file
        for line in eachline(file)
            # Split each line by tab character
            fields = split(line, ',')
            # Check if we have both key and value
            if length(fields) == 2
                key = strip(fields[1])    # First column
                value = strip(fields[2])   # Second column
                run_condition_dict[key] = value
            end
        end
    end
    
    return run_condition_dict
end