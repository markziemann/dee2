module MainHelpers exposing (..)

import Array
import Dict
import Http
import MainTypes exposing (..)
import Maybe.Extra as MExtra
import SearchBarTypes


extractSelectedRows : Array.Array SearchBarTypes.SearchResult -> String
extractSelectedRows rows =
    Array.map
        (\row ->
            if row.selected then
                List.map (\key -> Dict.get key row.data) ["SRR_accession"]

            else
                Nothing
        )
        rows
        |> Array.toList
        |> MExtra.values
        |> String.join ","


generateDownloadUrl selectedResults =

