module ResultsPage.Helpers exposing (..)

import Bool.Extra as BExtra
import Dict exposing (Dict)
import Dict.Extra as DExtra
import Html exposing (text)
import Html.Attributes exposing (class, style, title)
import Html.Events exposing (on)
import Json.Decode as Decode exposing (Decoder)
import ResultsPage.Types exposing (..)
import SearchPage.Types exposing (SearchResult)
import Table
import Url.Builder


noOverflowColumn : String -> (data -> ( Int, String )) -> Table.Column data Msg
noOverflowColumn name toStr =
    Table.veryCustomColumn
        { name = name
        , viewData = truncatedToolTip toStr
        , sorter = Table.unsortable
        }


truncatedToolTip toValues data =
    let
        ( id, string ) =
            toValues data
    in
    Table.HtmlDetails
        [ --class "text-truncate"
          style "overflow" "hidden"
        , style "text-overflow" "ellipsis"
        , on "mouseenter" (maybeShowToggleTip id)
        ]
        [ text string ]


maybeShowToggleTip : Int -> Decoder Msg
maybeShowToggleTip id =
    Decode.map (ShowToggleTip id) (Decode.field "target" Decode.value)



--Decode.map2 (\offsetWidth scrollWidth -> offsetWidth < scrollWidth)
--    (Decode.at ["target", "offsetWidth"] Decode.int)
--    (Decode.at ["target", "scrollWidth"] Decode.int)
--    |> Decode.andThen
--        (BExtra.ifElse
--            (Decode.succeed ShowToggleTip )
--            (Decode.fail "Text not truncated")
--        )


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
