module SearchBarViews exposing (..)

import Array exposing (isEmpty)
import Bool.Extra as BExtra
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


viewSearchButton : SearchMode -> Html Msg
viewSearchButton searchMode =
    let
        checkedWhen = (\mode -> mode == searchMode)
    in
    div [ class "d-flex justify-content-center" ]
        [ div [ class "btn-group dropright my-5" ]
            [ button
                [ onClick Search
                , class "btn btn-lg btn-outline-success"
                , type_ "button"
                ]
                [ text "Search" ]
            , button
                [ class "btn btn-success dropdown-toggle dropdown-toggle-split"
                , attribute "data-toggle" "dropdown"
                , attribute "aria-haspopup" "true"

                --, attribute "aria-expanded" "true" -- changes
                ]
                []
            , div
                [ class "dropdown-menu"

                --, attribute "x-placement" "bottom-start"
                ]
                [ div [ class "mx-3" ]
                    [ h6 [ class "dropdown-header" ] [ text "Search Mode" ]
                    , radioButton "strict-mode" "search-mode" "Strict" (checkedWhen Strict)
                    , radioButton "strict-mode" "search-mode" "Strict" (checkedWhen Fuzzy)
                    ]
                ]
            ]
        ]


radioButton : String -> String -> String -> Bool -> Msg -> Html Msg
radioButton idValue nameValue textLabel checked msg =
    let
        checkedString =
            BExtra.ifElse "checked" "" checked
    in
    div [ class "form-check" ]
        [ input
            [ onClick msg
            , attribute "checked" checkedString
            , id idValue
            , class "form-check-input"
            , type_ "radio"
            , name nameValue
            ]
            []
        , label [ class "form-check-label", for nameValue ] [ text textLabel ]
        ]
