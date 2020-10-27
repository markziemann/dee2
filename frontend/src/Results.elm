module Results exposing (..)

import Array exposing (Array, set)
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (onClick)
import Json.Decode as Decode exposing (Decoder)
import Http



-- Represents a single row in a .tsv file


type alias SearchData =
    List ( String, String )


type alias SearchResult =
    { id : Int
    , data : SearchData
    , selected : Bool
    }


type alias Model =
    Array SearchResult


listWrapped a =
    (::) a []


init =
    Array.empty


searchResultDecoder : Decoder Model
searchResultDecoder =
    Decode.list (Decode.keyValuePairs Decode.string)
        -- Decode.map doesn't iterate (confusing!) its more like function application
        -- it should be called apply
        |> Decode.map Array.fromList
        |> Decode.map (Array.indexedMap (\idx data -> SearchResult idx data False))


type Msg
    = GotHttpResults (Result Http.Error Model)
    | ResultClicked (Model -> Model)


update : Msg -> Model -> Model
update msg model =
    case msg of
        ResultClicked function ->
            function model

        GotHttpResults result ->
            Result.withDefault init result


selectClickedResult : SearchResult -> List (Html.Attribute Msg)
selectClickedResult ({ id, selected } as result) =
    [ class
        (if result.selected then
            "table-primary"

         else
            ""
        )
    , ResultClicked (\results -> Array.set id { result | selected = not result.selected } results)
        |> onClick
    ]


viewSearchResults : Model -> Html Msg
viewSearchResults searchResults =
    searchResults
        |> Array.map
            (\result ->
                tr (selectClickedResult result)
                    (List.map (\( key, value ) -> td [] [ text value ]) result.data)
            )
        |> Array.toList
        |> tbody []
        |> listWrapped
        |> table [ class "table table-hover table-sm table-bordered table-responsive" ]
