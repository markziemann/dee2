module Helpers exposing (..)

import Array
import Dict exposing (Dict)
import Dict.Extra as DExtra
import Maybe.Extra as MExtra
import SearchPage.Types
import Url.Builder


getKeyValue dict =
    case ( Dict.get "species" dict, Dict.get "SRR_accession" dict ) of
        ( Just key, Just value ) ->
            Just ( key, value )

        ( _, _ ) ->
            Nothing


extractSelectedRows : Array.Array SearchPage.Types.SearchResult -> List Url.Builder.QueryParameter
extractSelectedRows rows =
    Array.map
        (\row ->
            if row.selected then
                getKeyValue row.data

            else
                Nothing
        )
        rows
        |> Array.toList
        |> MExtra.values
        |> DExtra.fromListDedupe (\a b -> a ++ "," ++ b)
        |> Dict.toList
        |> List.map (\( a, b ) -> Url.Builder.string a b)


queryString : Maybe (Array.Array SearchPage.Types.SearchResult) -> String
queryString maybeRows =
    MExtra.unwrap "" (extractSelectedRows >> Url.Builder.toQuery) maybeRows
