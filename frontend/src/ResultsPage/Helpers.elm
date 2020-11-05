module ResultsPage.Helpers exposing (..)

import Array
import Dict exposing (Dict)
import Dict.Extra as DExtra
import Maybe.Extra as MExtra
import ResultsPage.Types exposing (..)
import SearchPage.Types
import Url.Builder


updateSearchData : Model -> SearchPage.Types.SearchOutMsg -> Model
updateSearchData model outMsg =
    { model
        | searchHits = Just outMsg.hits
        , searchResultRows = Just outMsg.rows
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
