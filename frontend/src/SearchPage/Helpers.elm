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
highlightMatchingText query suggestion =
    -- Note splitting a string removes the split string eg.> split "a" "James" = ["J", "mes"]
    if query /= "" then
        -- Determine location of matches (Case insensitive!)
        String.indexes (String.toLower query) (String.toLower suggestion)
            -- Get List of matching characters
            |> List.map (\idx -> String.slice idx (idx + String.length query) suggestion)
            -- Consecutively split the suggestion on the matching strings
            |> List.foldl (\str -> List.concatMap (\innerStr -> String.split str innerStr)) [ suggestion ]
            -- Convert each split string to NON-bold text
            |> List.map (\t -> p [ class "d-inline" ] [ text t ])
            -- Insert BOLD text of matching strings
            |> List.intersperse (p [ class "d-inline font-weight-bold" ] [ text query ])

    else
        []


withPagination : PaginationOffset -> SearchParameters -> SearchParameters
withPagination paginationOffset (SearchParameters level mode query _) =
    SearchParameters level mode query paginationOffset


withquery : String -> SearchParameters -> SearchParameters
withquery query (SearchParameters level mode _ paginationOffset) =
    SearchParameters level mode query paginationOffset


toSearchParameters : Model -> Maybe SearchParameters
toSearchParameters model =
    case model.level of
        Just level ->
            SearchParameters
                level
                model.mode
                model.query
                model.defaultPaginationOffset
                |> Just

        Nothing ->
            Nothing


withMode : Mode -> SearchParameters -> SearchParameters
withMode mode (SearchParameters level _ query paginationOffset) =
    SearchParameters level mode query paginationOffset


getQuery : Maybe SearchParameters -> Maybe String
getQuery parameters =
    case parameters of
        Just (SearchParameters _ _ query _) ->
            Just query

        Nothing ->
            Nothing


getMode : SearchParameters -> Mode
getMode (SearchParameters _ mode _ _) =
    mode

getLevel : SearchParameters -> Level
getLevel (SearchParameters level _ _ _) =
    level

differentSearch : Maybe SearchParameters -> SearchParameters -> Bool
differentSearch maybeSearchParameters (SearchParameters levelB modeB queryB _) =
    case maybeSearchParameters of
        Just (SearchParameters levelA modeA queryA _) ->
            levelA /= levelB || modeA /= modeB || queryA /= queryB

        Nothing ->
            False


updateLevel : Level -> Model -> Model
updateLevel level model =
    { model | level = Just level }
