module SearchBarViews exposing (..)

import Array exposing (isEmpty)
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (..)
import SearchBarHelpers exposing (highlightMatchingText, listWrapped)
import SearchBarTypes exposing (..)


viewLargeSearchBar : Model -> Html Msg
viewLargeSearchBar model =
    div []
        [ input
            [ onInput SearchUpdate
            , attribute "aria-label" "Search"
            , class "form-control form-control-lg"
            , placeholder "Human epilepisy | SRP070529"
            , type_ "search"
            , value model.searchString
            ]
            []
        , viewSuggestions model

        -- Alternate ^^^'viewSuggestions' func with no highlighting useful to debug
        --, ul [] (Array.toList (Array.map (\str -> li [][text str]) model.searchSuggestions))
        , div [ class "d-flex justify-content-center" ]
            [ button
                [ onClick Search
                , class "btn btn-lg btn-outline-success my-5"
                , type_ "button"
                ]
                [ text "Search" ]
            ]
        ]


viewSuggestions : Model -> Html Msg
viewSuggestions { suggestionsVisible, searchString, searchSuggestions, activeSuggestion } =
    let
        selector =
            case activeSuggestion of
                Nothing ->
                    \_ -> ""

                Just value ->
                    \idx ->
                        if idx == value then
                            "active"

                        else
                            ""

        show =
            if isEmpty searchSuggestions then
                identity

            else
                (\value ->
                    value
                |> (\str -> [ str, "show" ])
                |> (\strings -> if suggestionsVisible then strings else (::) "invisible" strings)
                |> String.join " " )
    in
    div [ class (show "dropdown") ]
        [ div [ class (show "dropdown-menu") ]
            (List.indexedMap
                (\idx suggestion ->
                    li
                        [ String.join " "
                            [ "dropdown-item"
                            , selector idx
                            ]
                            |> class
                        , onClick (SuggestionSelected idx)
                        ]
                        (highlightMatchingText searchString suggestion)
                )
                (Array.toList searchSuggestions)
            )
        ]


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
        |> table [ class "table table-hover table-sm table-bordered table-responsive" ]
