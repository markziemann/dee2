module SearchBarHelpers exposing (..)

import Array exposing (Array)
import Html exposing (..)
import Html.Attributes exposing (class)
import Json.Decode as Decode exposing (Decoder, array, field, string)
import Keyboard.Event exposing (KeyboardEvent)
import Process exposing (sleep)
import SearchBarTypes exposing (..)
import Task


setActiveSuggestion : Model -> Int -> Model
setActiveSuggestion model value =
    { model | activeSuggestion = Just value }


notWaiting : Model -> Model
notWaiting model =
    { model | waitingForResponse = False }


isWaiting : Model -> Model
isWaiting model =
    { model | waitingForResponse = True }


clearSearchSuggestions : Model -> Model
clearSearchSuggestions model =
    { model | searchSuggestions = Array.empty }


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


broadCast : msg -> Cmd msg
broadCast x =
    Task.perform identity (Task.succeed x)


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


listWrapped a =
    (::) a []


decodeSearchResults : Decoder SearchResults
decodeSearchResults =
    Decode.list (Decode.keyValuePairs Decode.string)
        -- Decode.map doesn't iterate (confusing!) its more like function application
        -- it should be called apply
        |> Decode.map Array.fromList
        |> Decode.map (Array.indexedMap (\idx data -> SearchResult idx data False))

--toKey : KeyboardEvent -> Maybe Msg
--toKey keyboardEvent =
--    -- If Nothing is returned the update
--    -- function is never called which simplifies
--    -- downstream logic
--    -- Visit: clause 5.6.2 for other keycodes
--    -- https://w3c.github.io/uievents/#dom-keyboardevent-key
--    case keyboardEvent.key of
--        Just value ->
--            if value == "ArrowUp" then
--                Just (KeyPressed ArrowUp)
--
--            else if value == "ArrowDown" then
--                Just (KeyPressed ArrowDown)
--
--            else if value == "Enter" then
--                Just (KeyPressed Enter)
--
--            else
--                Nothing
--
--        _ ->
--            Nothing
