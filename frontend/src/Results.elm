module Results exposing (..)

import Html exposing (..)
import Html.Attributes exposing (..)
import Json.Decode as Decode



-- Represents a single row in a .tsv file


type alias SearchData =
    List ( String, String )


type alias SearchResult =
    { data : SearchData
    , selected : Bool
    }


type alias SearchResults = List SearchResult

listWrap a = (::) a []

searchResultDecoder =
    Decode.list
        (Decode.map (\data -> SearchResult data False)
            (Decode.keyValuePairs Decode.string)
        )


highlightSelectedRow : Bool -> Html.Attribute msg
highlightSelectedRow selected =
    if selected then
        class "active"

    else
        class "table-primary"


viewSearchResults : SearchResults -> Html msg
viewSearchResults searchResults =
    searchResults
        |> List.map
            (\{data, selected}->
                tr [highlightSelectedRow selected]
                    (List.map (\( key, value ) -> td [] [ text value ]) data)
            )
        |> tbody []
        |> listWrap
        |> table [ class "table table-hover table-sm table-bordered table-responsive" ]
