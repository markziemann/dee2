module MainViews exposing (..)
import SearchBarTypes exposing (SearchResult, SearchResults)
import Html exposing (..)
import Html.Attributes as Attr
import Html.Events as Events
import SearchBarHelpers exposing (listWrapped)
import MainTypes exposing (..)
import Array

selectClickedResult : SearchResult -> List (Html.Attribute Msg)
selectClickedResult result =
    [ Attr.class
        (if result.selected then
            "table-primary"

         else
            ""
        )

    --    Probably shouldn't be passing functions to update
    , ResultClicked result
        |> Events.onClick
    ]


viewSearchResults : SearchResults -> Html Msg
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
        |> table [ Attr.class "table table-hover table-sm table-bordered table-responsive" ]