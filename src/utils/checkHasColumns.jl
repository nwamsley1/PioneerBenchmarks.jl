function checkHasColumns(
    expected_columns::Vector{Symbol},
    column_names::Set{Symbol}
)
    for column_name in expected_columns
        if column_name∉column_names
            throw("The column: "*column_name*" was not found.")
        end
    end
    return nothing
end
