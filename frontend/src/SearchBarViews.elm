module SearchBarViews exposing (..)

import Array exposing (isEmpty)
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (..)
import SearchBarHelpers exposing (highlightMatchingText)
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
            [text model.searchString]
        , viewSuggestions model

        -- Alternate ^^^'viewSuggestions' func with no highlighting useful to debug
        --, ul [] (Array.toList (Array.map (\str -> li [][text str]) model.searchSuggestions))
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
                \value ->
                    value
                        |> (\str -> [ str, "show" ])
                        |> (\strings ->
                                if suggestionsVisible then
                                    strings

                                else
                                    (::) "invisible" strings
                           )
                        |> String.join " "
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



