module SearchBar exposing (..)

import Array exposing (Array)
import Browser.Events exposing (onKeyDown)
import Elastic as Elastic exposing (parse, serialize)
import Http as Http exposing (get)
import Json.Decode as Decode
import KeyBoardHelpers exposing (arrowDown, arrowUp, orTry)
import Keyboard.Event exposing (considerKeyboardEvent)
import Result
import SearchBarHelpers exposing (..)
import SearchBarTypes as SearchBarTypes exposing (..)



-- Allows exporting of imported types


type alias Model =
    SearchBarTypes.Model


type alias Msg =
    SearchBarTypes.Msg



-- Expose this Msg constructor to other modules
searchMsg =
    SearchBarTypes.Search

suggestionSelectedMsg =
    SearchBarTypes.SuggestionSelected

init : Model
init =
    { searchString = ""
    , searchSuggestions = Array.empty
    , activeSuggestion = Nothing
    , suggestionsVisible = True
    , searchResults = Array.empty
    , waitingForResponse = False
    }


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        updateActiveSuggestion =
            \value -> setActiveSuggestion model value |> showSuggestions

        wrapAround =
            \idx ->
                idx
                    |> modBy (Array.length model.searchSuggestions)
                    |> updateActiveSuggestion
    in
    case msg of
        SearchUpdate value ->
            if value == "" then
                -- Don't wait to clear the suggestions if there is nothing
                -- in the search box. Having search suggestions with an empty
                -- searchString causes issues for the highlightMatchingText function
                ( { model | searchString = value }
                    |> clearActiveSuggestion
                    |> clearSearchSuggestions
                , Cmd.none
                )

            else
                ( { model | searchString = value }
                    |> showSuggestions
                    |> clearActiveSuggestion
                , delay 1000 (GetSearchSuggestions value)
                )

        Search ->
            let
                -- Elastic.Word makes no modification to the search string
                -- https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-simple-query-string-query.html#_simple_query_string_syntax
                search_string =
                    Result.withDefault (Elastic.Word model.searchString) (parse model.searchString)
                        |> serialize

                serverQuery =
                    get
                        { url = "/search/" ++ search_string
                        , expect = Http.expectJson GotHttpSearchResponse decodeSearchResults
                        }
            in
            ( model |> clearSearchSuggestions |> isWaiting, serverQuery )

        GetSearchSuggestions value ->
            -- This is some debounce on the search string to prevent spamming ElasticSearch with queries
            if model.searchString == value then
                ( model
                , get
                    { url = "/search_as_you_type/" ++ value
                    , expect = Http.expectJson GotSearchSuggestions decodeSearchSuggestions
                    }
                )

            else
                ( model |> showSuggestions, Cmd.none )

        GotSearchSuggestions result ->
            ( { model | searchSuggestions = Result.withDefault Array.empty result }
                |> showSuggestions
            , Cmd.none
            )

        SuggestionSelected value ->
            let
                searchString =
                    case Array.get value model.searchSuggestions of
                        Just string ->
                            string

                        Nothing ->
                            model.searchString
            in
            ( { model | searchString = searchString }
                |> clearActiveSuggestion
                |> hideSuggestions
            , Cmd.none
            )

        ArrowUp ->
            case model.activeSuggestion of
                Nothing ->
                    ( Array.length model.searchSuggestions
                        |> updateActiveSuggestion
                    , Cmd.none
                    )

                Just value ->
                    ( wrapAround (value - 1), Cmd.none )

        ArrowDown ->
            case model.activeSuggestion of
                Nothing ->
                    ( updateActiveSuggestion 0, Cmd.none )

                Just value ->
                    ( wrapAround (value + 1), Cmd.none )


        GotHttpSearchResponse result ->
            ( { model | searchResults = Result.withDefault Array.empty result }
                |> notWaiting
            , Cmd.none
            )

        ResultClicked function ->
            ( { model | searchResults = function model.searchResults }, Cmd.none )

        ClickOutOfSuggestions ->
            ( model |> hideSuggestions, Cmd.none )


subscriptions : Model -> Sub Msg
subscriptions model =
    Sub.batch
        [ onKeyDown (considerKeyboardEvent (arrowUp ArrowUp |> orTry (arrowDown ArrowDown)))
        , Browser.Events.onClick (Decode.succeed ClickOutOfSuggestions)
        ]
