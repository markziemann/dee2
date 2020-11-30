module ResultsPage.Helpers exposing (..)

import Bool.Extra as BExtra
import Dict exposing (Dict)
import Dict.Extra as DExtra
import Html exposing (text)
import Html.Attributes exposing (class, style)
import Html.Events exposing (on)
import Json.Decode as Decode exposing (Decoder)
import ResultsPage.Types exposing (..)
import SearchPage.Types exposing (SearchParameters, SearchResult, SearchResults)
import Table
import Url.Builder

maybeExpiredData: MaybeExpired a -> a
maybeExpiredData maybeExpired =
    case maybeExpired of
        Current a -> a
        Expired a -> a


expired : MaybeExpired a -> MaybeExpired a
expired maybeExpired =
    case maybeExpired of
        Current a ->
            Expired a
        _ ->
            maybeExpired


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
        [ class "text-truncate"
        , on "mouseenter" (maybeShowToggleTip id)
        ]
        [ text string ]


maybeShowToggleTip : Int -> Decoder Msg
maybeShowToggleTip id =
    Decode.map2 (\offsetWidth scrollWidth -> offsetWidth < scrollWidth)
        (Decode.at [ "target", "offsetWidth" ] Decode.float)
        (Decode.at [ "target", "scrollWidth" ] Decode.float)
        |> Decode.andThen
            (BExtra.ifElse
                (Decode.succeed (ShowToggleTip id))
                (Decode.fail "Text not truncated")
            )


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


stageResultForDownload : ColumnMapping -> SearchPage.Types.SearchResult -> Maybe ( String, String )
stageResultForDownload mapping searchResult =
    case ( Dict.get mapping.species searchResult.data, Dict.get mapping.accession searchResult.data ) of
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
