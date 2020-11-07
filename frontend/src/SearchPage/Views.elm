module SearchPage.Views exposing (..)

import Array exposing (isEmpty)
import Bool.Extra as BExtra
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (..)
import SearchPage.Helpers exposing (highlightMatchingText)
import SearchPage.Types exposing (..)
import SharedTypes


viewLargeSearchBar : Model -> Html Msg
viewLargeSearchBar model =
    div []
        [ input
            [ onInput SearchUpdate
            , attribute "aria-label" "Search"
            , class "form-control form-control-lg"
            , placeholder "e.g. Human epilepisy | SRP070529"
            , type_ "search"
            , value model.searchString
            , id "search-bar"
            ]
            [ text model.searchString ]
        , viewSuggestions model
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


viewSearchButton : Model -> SharedTypes.PaginationOffset -> Html Msg
viewSearchButton model paginationOffset =
    div [ class "d-flex justify-content-center" ]
        [ div [ class "btn-group dropright my-5" ]
            [ button
                [ onClick <| Search model.searchMode model.searchString paginationOffset
                , class "btn btn-lg btn-outline-success"
                , type_ "button"
                ]
                [ text "Search" ]
            ]
        ]


viewSearchModeSelector : SearchMode -> Html Msg
viewSearchModeSelector searchMode =
    div [ class "d-sm-flex justify-content-center" ]
        [ div [ class "form-check mx-2 my-4" ]
            [ input
                [ onInput StrictSelected
                , class "form-check-input"
                , type_ "radio"
                , name "search-mode"
                , checked <| BExtra.ifElse True False (searchMode == Strict)
                , id "simple-query-string"
                ]
                []
            , label [ class "form-check-label", for "simple-query-string" ] [ text "Simple Query String" ]
            ]
        , div [ class "form-check mx-2 my-4" ]
            [ input
                [ onInput FuzzySelected
                , class "form-check-input"
                , type_ "radio"
                , name "search-mode"
                , checked <| BExtra.ifElse True False (searchMode == Fuzzy)
                , id "fuzzy-text"
                ]
                []
            , label [ class "form-check-label", for "fuzzy-text" ] [ text "Fuzzy Text Search" ]
            ]
        ]
