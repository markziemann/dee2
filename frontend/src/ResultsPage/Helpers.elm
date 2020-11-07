module ResultsPage.Helpers exposing (..)

import Dict exposing (Dict)
import Dict.Extra as DExtra
import ResultsPage.Types exposing (..)
import SearchPage.Types
import Url.Builder
import Bool.Extra as BExtra
import SearchPage.Types exposing (SearchResult)
import Table

hideWhenTrue : String -> Bool -> String
hideWhenTrue classString true =
    BExtra.ifElse
        (String.join " " [ classString, "invisible" ])
        classString
        true

disableWhenTrue : String -> Bool -> String
disableWhenTrue classString true =
    BExtra.ifElse
        (String.join " " [ classString, "disabled" ])
        classString
        true

get : String -> SearchResult -> String
get key searchResult =
    Maybe.withDefault "" (Dict.get key searchResult.data)


getId =
    .id >> String.fromInt


defaultTable =
    Table.defaultCustomizations

highlightRowIfTrue style true =
    BExtra.ifElse style "" true

updateSearchData : Model -> SearchPage.Types.OutMsg -> Model
updateSearchData model outMsg =
    { model
        | searchResults = outMsg.searchResults
        , paginationOffset = outMsg.paginationOffset
    }



stageResultForDownload : SearchPage.Types.SearchResult -> Maybe ( String, String )
stageResultForDownload searchResult =
    case ( Dict.get "species" searchResult.data, Dict.get "SRR_accession" searchResult.data ) of
        ( Just key, Just value ) ->
            Just ( key, value )

        ( _, _ ) ->
            Nothing


extractSelectedRows : SelectedResults -> List Url.Builder.QueryParameter
extractSelectedRows selectedResults =
    Dict.toList selectedResults
        |> List.map Tuple.second
        |> DExtra.fromListDedupe (\a b -> a ++ "," ++ b)
        |> Dict.toList
        |> List.map (\( a, b ) -> Url.Builder.string a b)


queryString : SelectedResults -> String
queryString selectedResults =
    (extractSelectedRows >> Url.Builder.toQuery) selectedResults
