function parseMods(mods_string::AbstractString)::Base.RegexMatchIterator
    #Example: "1(7,M,Oxidation)(13,K,AnExampleMod)"
    mods_regex = r"(?<=\().*?(?=\))"
    return eachmatch(mods_regex, mods_string)
end

function parseMods(mods_string::Missing)::Base.RegexMatchIterator
    #Example: "1(7,M,Oxidation)(13,K,AnExampleMod)"
    mods_regex = r"(?<=\().*?(?=\))"
    return eachmatch(mods_regex, "")
end

function getModName(mod_string::AbstractString)::String
    match(r"[^,]+(?=$)", mod_string).match
end

function getModIndex(mod_string::AbstractString)::UInt8
    parse(UInt8, match(r"^[0-9]+(?=,)", mod_string).match)
end

function insert_at_indices(original::String, insertions::Vector{Tuple{String, UInt8}})
    # Convert the original string into an array of characters for easier manipulation
    char_array = collect(original)

    # Sort the insertions by index in ascending order
    sorted_insertions = sort(insertions, by = x -> x[2])

    # Adjust the index for each insertion
    offset = 0
    for (substr, idx) in sorted_insertions
        # Adjust the index with the current offset
        insertion_index = idx + offset
        # Insert each character of the substring at the specified index
        for (i, char) in enumerate(substr)
            insert!(char_array, insertion_index + i, char)
        end
        # Update the offset by the length of the inserted substring
        offset += length(substr)
    end

    # Join the array of characters back into a single string
    return join(char_array)
end

function getModifiedSequence(
    sequence::String,
    isotope_mods::String,
    structural_mods::String)

    mods = structural_mods*isotope_mods
    mods = [("["*uppercase(getModName(mod.match))*"]", getModIndex(mod.match)) for mod in parseMods(mods)]
    return insert_at_indices(sequence, mods)
end

function getModifiedSequence(
    sequence::String,
    isotope_mods::Missing,
    structural_mods::String)

    mods = structural_mods
    mods = [("["*uppercase(getModName(mod.match))*"]", getModIndex(mod.match)) for mod in parseMods(mods)]
    return insert_at_indices(sequence, mods)
end

function getModifiedSequence(
    sequence::String,
    isotope_mods::String,
    structural_mods::Missing)

    mods = isotope_mods
    mods = [("["*uppercase(getModName(mod.match))*"]", getModIndex(mod.match)) for mod in parseMods(mods)]
    return insert_at_indices(sequence, mods)
end