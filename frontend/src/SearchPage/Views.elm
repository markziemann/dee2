module SearchPage.Views exposing (..)

import Array exposing (isEmpty)
import Bool.Extra as BExtra
import Helpers exposing (errorToString)
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (..)
import SearchPage.Helpers exposing (highlightMatchingText, suggestionHighlightFunc)
import SearchPage.Types exposing (..)
import SharedTypes exposing (RemoteData(..))


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
        highlight =
            suggestionHighlightFunc activeSuggestion

        classString =
            if suggestionsVisible then
                "show"

            else
                "invisible"

        dropdown =
            \items ->
                div [ class <| "dropdown " ++ classString ]
                    [ div [ class <| "dropdown-menu " ++ classString ] items ]
    in
    case searchSuggestions of
        Success suggestions ->
            List.indexedMap
                (\idx suggestion ->
                    li
                        [ String.join " "
                            [ "dropdown-item"
                            , highlight idx
                            ]
                            |> class
                        , onClick (SuggestionSelected idx)
                        ]
                        (highlightMatchingText searchString suggestion)
                )
                (Array.toList suggestions)
                |> dropdown

        Failure err ->
            [ div [ class "text-warning" ] [ text <| errorToString err ] ]
                |> dropdown

        NotAsked ->
            dropdown []

        Loading ->
            [ span
                [ class "spinner-border spinner-border-sm mr-2"
                , attribute "role" "status"
                , attribute "aria-hidden" "true"
                ]
                []
            , text "Loading..."
            ]
                |> dropdown


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
