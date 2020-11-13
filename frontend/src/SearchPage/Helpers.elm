module SearchPage.Helpers exposing (..)

import Array exposing (Array)
import Html exposing (..)
import Html.Attributes exposing (class)
import Json.Decode as Decode exposing (Decoder, array, field, string)
import Process exposing (sleep)
import SearchPage.Types exposing (..)
import SharedTypes exposing (PaginationOffset, RemoteData(..), WebData)
import Task


getWebDataIndex : Int -> WebData (Array.Array a) -> Maybe a
getWebDataIndex index webData =
    case webData of
        Success array ->
            Array.get index array

        _ ->
            Nothing


lengthWebData : WebData (Array.Array a) -> Maybe Int
lengthWebData webData =
    case webData of
        Success array ->
            Just (Array.length array)

        _ ->
            Nothing


emptyWebData : WebData (Array.Array a) -> Bool
emptyWebData webData =
    case webData of
        Success array ->
            Array.isEmpty array

        _ ->
            False


updateActiveSuggestion : Model -> Int -> Model
updateActiveSuggestion model value =
    { model | activeSuggestion = Just value }


clearSearchSuggestions : Model -> Model
clearSearchSuggestions model =
    { model | searchSuggestions = NotAsked }


clearActiveSuggestion : Model -> Model
clearActiveSuggestion model =
    { model | activeSuggestion = Nothing }


showSuggestions : Model -> Model
showSuggestions model =
    { model | suggestionsVisible = True }


hideSuggestions : Model -> Model
hideSuggestions model =
    { model | suggestionsVisible = False }


delay : Float -> msg -> Cmd msg
delay time msg =
    sleep time
        |> Task.andThen (always <| Task.succeed msg)
        |> Task.perform identity


decodeSearchSuggestions : Decoder SearchSuggestions
decodeSearchSuggestions =
    field "suggestions" (array string)


suggestionHighlightFunc : Maybe Int -> (Int -> String)
suggestionHighlightFunc maybeActiveSuggestion =
    case maybeActiveSuggestion of
        Nothing ->
            \_ -> ""

        Just value ->
            \idx ->
                if idx == value then
                    "active"

                else
                    ""


highlightMatchingText : String -> String -> List (Html msg)
highlightMatchingText searchString suggestion =
    -- Note splitting a string removes the split string eg.> split "a" "James" = ["J", "mes"]
    if searchString /= "" then
        -- Determine location of matches (Case insensitive!)
        String.indexes (String.toLower searchString) (String.toLower suggestion)
            -- Get List of matching characters
            |> List.map (\idx -> String.slice idx (idx + String.length searchString) suggestion)
            -- Consecutively split the suggestion on the matching strings
            |> List.foldl (\str -> List.concatMap (\innerStr -> String.split str innerStr)) [ suggestion ]
            -- Convert each split string to NON-bold text
            |> List.map (\t -> p [ class "d-inline" ] [ text t ])
            -- Insert BOLD text of matching strings
            |> List.intersperse (p [ class "d-inline font-weight-bold" ] [ text searchString ])

    else
        []




withPagination : PaginationOffset -> SearchParameters -> SearchParameters
withPagination paginationOffset (SearchParameters searchMode searchString _) =
    SearchParameters searchMode searchString paginationOffset


withSearchString : String -> SearchParameters -> SearchParameters
withSearchString searchString (SearchParameters searchMode _ paginationOffset) =
    SearchParameters searchMode searchString paginationOffset


withSearchMode : SearchMode -> SearchParameters -> SearchParameters
withSearchMode searchMode (SearchParameters _ searchString paginationOffset) =
    SearchParameters searchMode searchString paginationOffset


getSearchString : SearchParameters -> String
getSearchString (SearchParameters _ searchString _) =
    searchString


getSearchMode : SearchParameters -> SearchMode
getSearchMode (SearchParameters searchMode _ _) =
    searchMode


sameModeAndString : SearchParameters -> SearchParameters -> Bool
sameModeAndString (SearchParameters modeA stringA _) (SearchParameters modeB stringB _) =
    modeA == modeB && stringA == stringB
